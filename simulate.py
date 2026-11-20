"""
Grid Inertia, RoCoF & Frequency Stability RK4 Simulator (Devlog #10)
Single-bus swing-equation RK4 numerical integration (dt = 0.001 s) with:
- Post-trip synchronous unit inertia loss
- Mechanical governor deadband & reserve saturation
- Fair mechanism separation: Virtual Machine Mode (VMM) vs Fast Frequency Response (FFR <200 ms) vs Synchronous Condenser
- TEIAS Grid Code Art. 67 (49.0 Hz first UFLS stage) & EirGrid 75% SNSP benchmark
"""

import numpy as np
import matplotlib.pyplot as plt


def simulate(
    p_base: float = 32000.0,
    delta_p: float = -1200.0,
    ibr_share: float = 0.75,
    h_sync: float = 4.8,
    h_override: float | None = None,
    trip_loses_inertia: bool = True,
    h_extra: float = 0.0,
    h_vmm_mws: float = 0.0,
    tau_gov: float = 4.2,
    gov_delay: float = 0.0,
    deadband_hz: float = 0.02,
    d_frac: float = 0.01,
    k_gov: float = 3200.0,
    r_gov_max: float = 1800.0,
    bess_mw: float = 0.0,
    tau_gfm: float = 0.15,
    dt: float = 0.001,
    t_max: float = 15.0,
) -> dict:
    f0 = 50.0
    h_pre = h_override if h_override is not None else h_sync * (1.0 - ibr_share)
    h_trip_loss = (abs(delta_p) * h_sync / p_base) if trip_loses_inertia else 0.0
    h_sys = max(0.05, h_pre - h_trip_loss + h_extra + (h_vmm_mws / p_base))

    d_mw_per_hz = d_frac * p_base
    rocof_init = (f0 * delta_p) / (2.0 * h_sys * p_base)
    m_sys = (2.0 * h_sys * p_base) / f0  # MW*s / Hz

    def rhs(t_val: float, f_val: float, p_g: float, p_b: float):
        df_err = max(0.0, f0 - f_val)

        if t_val >= gov_delay and df_err > deadband_hz:
            target_gov = min(r_gov_max, k_gov * (df_err - deadband_hz))
        else:
            target_gov = 0.0
        dp_gov = (target_gov - p_g) / tau_gov

        if bess_mw > 0.0 and df_err > 0.03:
            target_bess = min(bess_mw, (df_err / 0.25) * bess_mw)
        else:
            target_bess = 0.0
        dp_bess = (target_bess - p_b) / tau_gfm

        p_imbalance = delta_p + p_g + p_b + d_mw_per_hz * df_err
        df_dt = p_imbalance / m_sys
        return df_dt, dp_gov, dp_bess

    f = f0
    p_gov = 0.0
    p_bess = 0.0

    f_nadir = f0
    t_nadir = 0.0
    t_to_49 = None

    n_steps = int(round(t_max / dt))
    stride = 20
    t_pts, f_pts = [], []

    for step in range(n_steps + 1):
        t = step * dt
        if f < f_nadir:
            f_nadir = f
            t_nadir = t
        if t_to_49 is None and f <= 49.0:
            t_to_49 = t
        if step % stride == 0 or step == n_steps:
            t_pts.append(t)
            f_pts.append(f)
        if step == n_steps:
            break

        k1_f, k1_g, k1_b = rhs(t, f, p_gov, p_bess)
        k2_f, k2_g, k2_b = rhs(
            t + 0.5 * dt, f + 0.5 * dt * k1_f, p_gov + 0.5 * dt * k1_g, p_bess + 0.5 * dt * k1_b
        )
        k3_f, k3_g, k3_b = rhs(
            t + 0.5 * dt, f + 0.5 * dt * k2_f, p_gov + 0.5 * dt * k2_g, p_bess + 0.5 * dt * k2_b
        )
        k4_f, k4_g, k4_b = rhs(
            t + dt, f + dt * k3_f, p_gov + dt * k3_g, p_bess + dt * k3_b
        )

        f += (dt / 6.0) * (k1_f + 2.0 * k2_f + 2.0 * k3_f + k4_f)
        p_gov += (dt / 6.0) * (k1_g + 2.0 * k2_g + 2.0 * k3_g + k4_g)
        p_bess += (dt / 6.0) * (k1_b + 2.0 * k2_b + 2.0 * k3_b + k4_b)

    return {
        "h_sys": h_sys,
        "rocof": rocof_init,
        "f_nadir": f_nadir,
        "t_nadir": t_nadir,
        "f_steady": f,
        "t_to_49": t_to_49,
        "ufls_49": f_nadir <= 49.0,
        "t": np.asarray(t_pts),
        "f": np.asarray(f_pts),
    }


def fmt_row(label: str, r: dict) -> str:
    t49 = f"{r['t_to_49']:.2f}s" if r["t_to_49"] is not None else "None"
    return (
        f"{label:<44} | H={r['h_sys']:4.2f}s | RoCoF={r['rocof']:+6.3f} Hz/s | "
        f"Nadir={r['f_nadir']:6.3f} Hz (@{r['t_nadir']:4.2f}s) | 49.0Hz={t49:<6} | UFLS={str(r['ufls_49']):<5} | Steady={r['f_steady']:5.2f} Hz"
    )


def plot_charts(results: dict) -> None:
    plt.style.use("dark_background")

    # Chart 1: ENTSO-E CE vs Island Mode IBR Sweep
    fig, ax = plt.subplots(figsize=(10, 5.6), dpi=160)
    fig.patch.set_facecolor("#0a101c")
    ax.set_facecolor("#0f172a")

    ax.plot(results["ce"]["t"], results["ce"]["f"], color="#38bdf8", lw=2.4, ls="--", label="ENTSO-E CE Connected (350 GW) -> 49.897 Hz")
    ax.plot(results["ibr20"]["t"], results["ibr20"]["f"], color="#34d399", lw=2.2, label="Island 20% IBR (H=3.66 s) -> 49.348 Hz")
    ax.plot(results["ibr50"]["t"], results["ibr50"]["f"], color="#facc15", lw=2.2, label="Island 50% IBR (H=2.22 s) -> 49.184 Hz")
    ax.plot(results["ibr75"]["t"], results["ibr75"]["f"], color="#f43f5e", lw=2.8, label="Island 75% IBR (H=1.02 s) -> 48.807 Hz (UFLS @ 1.59 s)")

    ax.axhline(49.0, color="#ef4444", ls=":", lw=1.8, label="TEIAS Art. 67 First UFLS Stage (49.0 Hz)")
    ax.set_title("Grid Frequency Response to -1,200 MW Contingency (RK4, dt = 0.001 s)", fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel("Time after contingency (s)")
    ax.set_ylabel("Grid Frequency (Hz)")
    ax.set_xlim(0, 15)
    ax.set_ylim(48.65, 50.08)
    ax.grid(True, alpha=0.2)
    ax.legend(loc="lower right", frameon=True, facecolor="#0f172a", edgecolor="#334155", fontsize=9)
    plt.tight_layout()
    plt.savefig("frequency_nadir_comparison.png")
    plt.close(fig)

    # Chart 2: Fair Mechanism Separation (VMM vs Sync Condenser vs FFR vs GFM-BESS)
    fig, ax = plt.subplots(figsize=(10, 5.6), dpi=160)
    fig.patch.set_facecolor("#0a101c")
    ax.set_facecolor("#0f172a")

    ax.plot(results["ibr75"]["t"], results["ibr75"]["f"], color="#f43f5e", lw=2.0, ls="--", alpha=0.85, label="Baseline 75% IBR (No Mitigation) -> 48.807 Hz")
    ax.plot(results["vmm"]["t"], results["vmm"]["f"], color="#f97316", lw=2.2, label="VMM Virtual Inertia Only (+6,000 MWs) -> 48.901 Hz (UFLS)")
    ax.plot(results["cond"]["t"], results["cond"]["f"], color="#38bdf8", lw=2.4, label="Synchronous Condenser (+25,600 MWs) -> 49.102 Hz")
    ax.plot(results["ffr"]["t"], results["ffr"]["f"], color="#facc15", lw=2.4, label="FFR Only (450 MW, tau=150 ms) -> 49.336 Hz")
    ax.plot(results["gfm"]["t"], results["gfm"]["f"], color="#34d399", lw=2.8, label="Combined GFM-BESS (450 MW FFR + 6,000 MWs VMM) -> 49.379 Hz")

    ax.axhline(49.0, color="#ef4444", ls=":", lw=1.8, label="TEIAS Art. 67 First UFLS Stage (49.0 Hz)")
    ax.set_title("75% IBR Island Mode Mitigation: VMM vs Sync Condenser vs FFR vs GFM-BESS", fontsize=12.5, fontweight="bold", pad=12)
    ax.set_xlabel("Time after contingency (s)")
    ax.set_ylabel("Grid Frequency (Hz)")
    ax.set_xlim(0, 15)
    ax.set_ylim(48.65, 50.08)
    ax.grid(True, alpha=0.2)
    ax.legend(loc="lower right", frameon=True, facecolor="#0f172a", edgecolor="#334155", fontsize=9)
    plt.tight_layout()
    plt.savefig("ffr_vs_vmm_solution_comparison.png")
    plt.close(fig)


def main() -> None:
    results = {
        "ce": simulate(p_base=350000.0, delta_p=-1200.0, h_override=3.5, trip_loses_inertia=True, tau_gov=4.2, k_gov=15000.0, r_gov_max=3000.0),
        "ibr20": simulate(ibr_share=0.20, trip_loses_inertia=True),
        "ibr50": simulate(ibr_share=0.50, trip_loses_inertia=True),
        "ibr75": simulate(ibr_share=0.75, trip_loses_inertia=True),
        "vmm": simulate(ibr_share=0.75, trip_loses_inertia=True, bess_mw=0.0, h_vmm_mws=6000.0),
        "cond": simulate(ibr_share=0.75, trip_loses_inertia=True, bess_mw=0.0, h_vmm_mws=25600.0),
        "ffr": simulate(ibr_share=0.75, trip_loses_inertia=True, bess_mw=450.0, tau_gfm=0.15, h_vmm_mws=0.0),
        "gfm": simulate(ibr_share=0.75, trip_loses_inertia=True, bess_mw=450.0, tau_gfm=0.15, h_vmm_mws=6000.0),
    }

    print("=== TABLE 1: ENTSO-E CE CONNECTED vs ISLAND MODE (32 GW, -1,200 MW CONTINGENCY) ===")
    print(fmt_row("ENTSO-E CE Connected (350 GW)", results["ce"]))
    print(fmt_row("Island Mode (20% IBR)", results["ibr20"]))
    print(fmt_row("Island Mode (50% IBR)", results["ibr50"]))
    print(fmt_row("Island Mode (75% IBR, No Mitigation)", results["ibr75"]))

    print("\n=== TABLE 2: MITIGATION COMPARISON (75% IBR ISLAND MODE) ===")
    print(fmt_row("VMM Virtual Inertia Only (+6,000 MWs)", results["vmm"]))
    print(fmt_row("Synchronous Condenser (+25,600 MWs)", results["cond"]))
    print(fmt_row("FFR Only (450 MW, tau=150 ms < 200 ms)", results["ffr"]))
    print(fmt_row("Combined GFM-BESS (450 MW FFR + 6,000 MWs)", results["gfm"]))

    plot_charts(results)
    print("\n[OK] Generated frequency_nadir_comparison.png and ffr_vs_vmm_solution_comparison.png")


if __name__ == "__main__":
    main()
