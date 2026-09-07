"""M5 figures: payload-vs-inclination headline curve and rotational-assistance
explanation, built from the authoritative CSV (scripts/m5_finalize_results.py).

NOTE on scope (DESIGN.md M5 S13): a separate, denser "fine sweep" across many more
inclinations was originally planned for the headline curve, but was dropped after the
near-polar guidance search demonstrated severe, genuine numerical cost (a single
90 deg trajectory costs 15-30x a 28.5 deg one even after the V_FLOOR/deadlock fixes).
The headline figure below uses the 6 independently re-verified authoritative points
only; DESIGN.md documents this as an explicit, deliberate scope reduction, not a
silently smaller result set.

Usage:
    python scripts/m5_figures.py
Writes: figures/m5_payload_vs_inclination.png, figures/m5_rotational_assistance.png
"""

import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(__file__)


def load_csv(name):
    with open(os.path.join(HERE, name)) as f:
        return list(csv.DictReader(f))


def main():
    rows = load_csv("m5_inclination_sweep_results.csv")
    i_arr = np.array([float(r["inclination_deg"]) for r in rows])
    p_arr = np.array([float(r["max_payload_kg"]) for r in rows])
    boost_arr = np.array([float(r["useful_rotational_boost_ms"]) for r in rows])

    baseline = p_arr[i_arr == 28.5][0]
    polar = p_arr[i_arr == 90.0][0]

    # --- Figure 1: headline payload vs inclination ---------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(i_arr, p_arr, "o-", color="tab:blue", markersize=8,
            label="authoritative points (independently re-verified)")
    ax.axhline(baseline, color="gray", linestyle="--", linewidth=1,
               label=f"28.5 deg baseline ({baseline:.0f} kg)")
    ax.scatter([28.5], [baseline], color="tab:green", zorder=6, s=140, marker="*",
               label="28.5 deg (M4 baseline)")
    ax.scatter([90.0], [polar], color="tab:purple", zorder=6, s=140, marker="*",
               label=f"90 deg polar endpoint ({polar:.0f} kg)")
    for x, y in zip(i_arr, p_arr):
        ax.annotate(f"{y:.0f}", (x, y), textcoords="offset points", xytext=(0, 10),
                    fontsize=8, ha="center")
    ax.set_xlabel("target inclination [deg]  (direct-ascent domain: 28.5-90 deg)")
    ax.set_ylabel("maximum payload to 400 km [kg]")
    ax.set_ylim(min(p_arr) - 150, max(p_arr) + 200)
    ax.set_title("M5 headline: maximum payload to 400 km vs. target inclination\n"
                  "(direct-ascent, M4 orbit-capable study vehicle, launch site 28.5 deg N)\n"
                  "NOT strictly monotonic -- see DESIGN.md M5 S9 for why", fontsize=10)
    ax.legend(fontsize=8, loc="lower left")

    # Secondary top axis: launch azimuth is a strictly DEcreasing function of
    # inclination (90 deg at i=28.5 down to 0 deg at i=90), so match that direction
    # explicitly rather than relying on an inverse-function auto-mapping (which
    # produced overlapping/garbled tick labels near the axis ends when tried).
    azimuth_deg = np.degrees(np.arcsin(np.clip(
        np.cos(np.radians(i_arr)) / np.cos(np.radians(28.5)), -1, 1)))
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(i_arr)
    ax2.set_xticklabels([f"{a:.0f}" for a in azimuth_deg])
    ax2.set_xlabel("direct-ascent launch azimuth [deg from north]")

    fig.tight_layout()
    out1 = os.path.join(HERE, "..", "figures", "m5_payload_vs_inclination.png")
    fig.savefig(out1, dpi=150)
    print(f"figure saved to {os.path.abspath(out1)}")

    # --- Figure 2: rotational-assistance explanation --------------------------------
    fig2, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    ax.plot(i_arr, boost_arr, "o-", color="tab:orange")
    ax.set_xlabel("target inclination [deg]")
    ax.set_ylabel("useful (in-plane) rotational boost [m/s]")
    ax.set_title("Earth-rotation assistance vs. inclination\n(strictly monotonic -- pure geometry)")

    ax = axes[1]
    boost_loss = boost_arr[0] - boost_arr
    payload_loss = p_arr[0] - p_arr
    ax.scatter(boost_loss, payload_loss, color="tab:red", zorder=5, s=60)
    for x, y, i in zip(boost_loss, payload_loss, i_arr):
        ax.annotate(f"{i:.1f} deg", (x, y), textcoords="offset points", xytext=(6, 6),
                    fontsize=8)
    ax.axhline(0, color="gray", linewidth=0.8)
    ax.set_xlabel("lost rotational assistance vs. 28.5 deg [m/s]")
    ax.set_ylabel("payload loss vs. 28.5 deg [kg]")
    ax.set_title("Direct-ascent inclination penalty vs. lost rotational assistance\n"
                  "(weak/no correlation here -- guidance-search landscape dominates;\n"
                  "NOT a plane-change delta-v -- this is a direct-ascent azimuth trade)",
                  fontsize=9)

    fig2.tight_layout()
    out2 = os.path.join(HERE, "..", "figures", "m5_rotational_assistance.png")
    fig2.savefig(out2, dpi=150)
    print(f"figure saved to {os.path.abspath(out2)}")


if __name__ == "__main__":
    main()
