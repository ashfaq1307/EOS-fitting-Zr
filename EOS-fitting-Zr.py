import numpy as np
from pathlib import Path
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")

SCRIPT_DIR = Path(__file__).resolve().parent
DAT_FILE   = SCRIPT_DIR / "Zr_EV.dat"
OUT_DIR    = SCRIPT_DIR

# Unit conversion
BOHR3_TO_ANG3  = 0.148185       # 1 bohr^3 = 0.148185 Å³
RY_PER_BOHR3_TO_GPA = 14710.5   # 1 Ry/bohr^3 = 14710.5 GPa

# Load data
data = np.loadtxt(DAT_FILE, comments="#")
V_raw = data[:, 0]   # bohr^3/atom  (descending – sort ascending)
E_raw = data[:, 1]   # Ry/atom

idx = np.argsort(V_raw)
V = V_raw[idx]
E = E_raw[idx]

V_ang = V * BOHR3_TO_ANG3   # for display only

# Initial guess from data
i0      = np.argmin(E)
E0_init = E[i0]
V0_init = V[i0]
B0_init = 0.015   # Ry/bohr^3  ≈ 100 GPa
BP_init = 4.0


#  EOS DEFINITIONS  (all in bohr^3 / Ry)

#  E0: Energy at equilibrium (Ry/atom)
#  V0: Equilibrium volume (bohr^3/atom)
#  B0: Bulk modulus at equilibrium (Ry/bohr^3)
#  BP: Pressure derivative of bulk modulus at equilibrium (dimensionless)

def murnaghan(V, E0, V0, B0, BP):
    eta = V0 / V
    return E0 + B0 * V0 / BP * (eta**BP / (BP - 1) + 1) - B0 * V0 / (BP - 1)

def birch_murnaghan(V, E0, V0, B0, BP):
    f = 0.5 * ((V0 / V)**(2/3) - 1)
    xi = 3/4 * (BP - 4)
    return E0 + 9 * V0 * B0 / 16 * (f**2 * (6 - 4 * f + xi * 9 * f**2))

def vinet_rose(V, E0, V0, B0, BP):
    x = (V / V0)**(1/3)
    eta = 1.5 * (BP - 1)
    return E0 + 2 * B0 * V0 / (BP - 1)**2 * (
        2 - (5 + 3 * BP * (x - 1) - 3 * x) * np.exp(-eta * (x - 1))
    )

def poirier_tarantola(V, E0, V0, B0, BP):
    x = np.log(V / V0)
    h = B0 * V0
    xi = 1 + 0.5 * (BP - 2) * x
    return E0 + h / 2 * x**2 * xi

# Bounds / parameter order
p0 = [E0_init, V0_init, B0_init, BP_init]
bounds = ([-np.inf, 0, 0, 0], [np.inf, np.inf, np.inf, 20])

MODELS = {
    "Murnaghan": murnaghan,
    "Birch-Murnaghan": birch_murnaghan,
    "Vinet-Rose": vinet_rose,
    "Poirier-Tarantola": poirier_tarantola,
}


#  FIT ALL MODELS

results = {}
print("\n" + "-"*72)
print("   COMPARATIVE EOS FITTING  -  Zirconium (Zr)")
print("-"*72)
header = f"{'Model':<22} {'E0 (Ry/at)':>12} {'V0 (A3/at)':>11} {'B0 (GPa)':>10} {'B0p':>7}  {'R2 (%)':>8}  {'RMSE (mRy)':>10}"
print(header)
print("-"*72)

for name, func in MODELS.items():
    try:
        popt, pcov = curve_fit(func, V, E, p0=p0, bounds=bounds, maxfev=50000)
        E_fit   = func(V, *popt)
        residuals = E - E_fit
        ss_res  = np.sum(residuals**2)
        ss_tot  = np.sum((E - E.mean())**2)
        r2      = 1 - ss_res / ss_tot
        rmse    = np.sqrt(np.mean(residuals**2)) * 1000   # mRy

        E0, V0b, B0b, BP = popt
        V0_ang3 = V0b * BOHR3_TO_ANG3
        B0_GPa  = B0b * RY_PER_BOHR3_TO_GPA

        results[name] = dict(
            popt=popt, pcov=pcov, E_fit=E_fit,
            E0=E0, V0_bohr=V0b, V0_ang=V0_ang3,
            B0_GPa=B0_GPa, BP=BP,
            r2=r2, rmse=rmse
        )
        print(f"{name:<22} {E0:>12.6f} {V0_ang3:>11.4f} {B0_GPa:>10.3f} {BP:>7.3f}  {r2*100:>8.5f}  {rmse:>10.5f}")

    except Exception as exc:
        print(f"{name:<22}  FAILED: {exc}")

print("-"*72)

# Best model by R^2
best_name = max(results, key=lambda k: results[k]["r2"])
print(f"\n    Best fit: {best_name}  (R² = {results[best_name]['r2']*100:.6f} %)\n")


#  FIGURE 1 – E vs V with all four fitted curves

V_fine     = np.linspace(V.min() * 0.98, V.max() * 1.02, 500)
V_fine_ang = V_fine * BOHR3_TO_ANG3

COLORS = {
    "Murnaghan": "#E63946",
    "Birch-Murnaghan": "#457B9D",
    "Vinet-Rose": "#2A9D8F",
    "Poirier-Tarantola": "#E9C46A",
}
LSTYLE = {
    "Murnaghan": "-",
    "Birch-Murnaghan": "--",
    "Vinet-Rose": "-.",
    "Poirier-Tarantola": ":",
}

fig1, ax1 = plt.subplots(figsize=(9, 5.5))
ax1.scatter(V_ang, E, s=45, zorder=5, color="k", label="DFT data", marker="o")

for name, res in results.items():
    E_fine = MODELS[name](V_fine, *res["popt"])
    ax1.plot(V_fine_ang, E_fine, color=COLORS[name], ls=LSTYLE[name], lw=2.0, label=f"{name}  (R²={res['r2']*100:.4f}%)")

ax1.set_xlabel("Volume  (Å^3/atom)", fontsize=13)
ax1.set_ylabel("Energy  (Ry/atom)", fontsize=13)
ax1.set_title("Equation-of-State Fitting — Zirconium", fontsize=14, fontweight="bold")
ax1.legend(fontsize=9.5, loc="upper right")
ax1.grid(alpha=0.35, ls=":")
fig1.tight_layout()
fig1.savefig(OUT_DIR / "Zr_EOS_fit.png", dpi=180)
plt.close(fig1)


#  FIGURE 2 – Residuals panel (4 sub-plots)

fig2, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
axes = axes.ravel()

for ax, (name, res) in zip(axes, results.items()):
    res_mRy = (E - res["E_fit"]) * 1000
    ax.bar(V_ang, res_mRy, width=1.5, color=COLORS[name], alpha=0.8, edgecolor="k", lw=0.5)
    ax.axhline(0, color="k", lw=1.2)
    ax.set_title(name, fontsize=11, fontweight="bold")
    ax.set_ylabel("Residual (mRy)", fontsize=9)
    ax.set_xlabel("Volume (Å^3/atom)", fontsize=9)
    ax.annotate(f"RMSE = {res['rmse']:.4f} mRy\nR² = {res['r2']*100:.5f} %", xy=(0.04, 0.92), xycoords="axes fraction", fontsize=8.5, va="top", bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
    ax.grid(alpha=0.3, ls=":")

fig2.suptitle("Residuals of EOS Fits — Zirconium", fontsize=14, fontweight="bold")
fig2.tight_layout()
fig2.savefig(OUT_DIR / "Zr_EOS_residuals.png", dpi=180)
plt.close(fig2)


#  FIGURE 3 – Parameter comparison bar chart

names  = list(results.keys())
V0s    = [results[n]["V0_ang"] for n in names]
B0s    = [results[n]["B0_GPa"] for n in names]
BPs    = [results[n]["BP"]     for n in names]

fig3, axes3 = plt.subplots(1, 3, figsize=(13, 4.5))
params = [(V0s, "V₀  (Å^3/atom)", "#4A90D9"),
          (B0s, "B₀  (GPa)", "#E05C5C"),
          (BPs, "B₀'", "#50C878")]

for ax, (vals, ylabel, color) in zip(axes3, params):
    bars = ax.bar(names, vals, color=color, edgecolor="k", lw=0.8, alpha=0.9)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(ylabel, fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", labelsize=8.5, rotation=15)
    ax.grid(axis="y", alpha=0.35, ls=":")
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.005, f"{val:.3f}", ha="center", va="bottom", fontsize=8)

fig3.suptitle("Fitted EOS Parameters — Zirconium", fontsize=13, fontweight="bold")
fig3.tight_layout()
fig3.savefig(OUT_DIR / "Zr_EOS_parameters.png", dpi=180)
plt.close(fig3)


#  FIGURE 4 – R^2 and RMSE comparison

r2s   = [results[n]["r2"]*100    for n in names]
rmses = [results[n]["rmse"]      for n in names]

fig4, (ax_r2, ax_rmse) = plt.subplots(1, 2, figsize=(10, 4))

bars_r2 = ax_r2.bar(names, r2s, color=["#E63946","#457B9D","#2A9D8F","#E9C46A"], edgecolor="k", lw=0.8)
ax_r2.set_ylabel("R^2  (%)", fontsize=11)
ax_r2.set_title("Goodness of Fit — R^2", fontsize=11, fontweight="bold")
ax_r2.set_ylim([min(r2s)*0.9999995, 100.0000005])
ax_r2.tick_params(axis="x", labelsize=9, rotation=12)
ax_r2.grid(axis="y", alpha=0.3, ls=":")
for bar, val in zip(bars_r2, r2s):
    ax_r2.text(bar.get_x() + bar.get_width()/2, bar.get_height(), f"{val:.5f}", ha="center", va="bottom", fontsize=8)

bars_rm = ax_rmse.bar(names, rmses, color=["#E63946","#457B9D","#2A9D8F","#E9C46A"], edgecolor="k", lw=0.8)
ax_rmse.set_ylabel("RMSE  (mRy)", fontsize=11)
ax_rmse.set_title("Root-Mean-Square Error", fontsize=11, fontweight="bold")
ax_rmse.tick_params(axis="x", labelsize=9, rotation=12)
ax_rmse.grid(axis="y", alpha=0.3, ls=":")
for bar, val in zip(bars_rm, rmses):
    ax_rmse.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.02, f"{val:.4f}", ha="center", va="bottom", fontsize=8)

fig4.suptitle("Fit Quality Comparison — Zirconium EOS", fontsize=13, fontweight="bold")
fig4.tight_layout()
fig4.savefig(OUT_DIR / "Zr_EOS_quality.png", dpi=180)
plt.close(fig4)

print("  Figures saved:")
print("   • Zr_EOS_fit.png")
print("   • Zr_EOS_residuals.png")
print("   • Zr_EOS_parameters.png")
print("   • Zr_EOS_quality.png\n")
