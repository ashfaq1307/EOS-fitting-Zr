# Comparative EOS Fitting — Zirconium (Zr)

A Python-based curve-fitting tool for analyzing **Energy–Volume (E–V) data** of Zirconium using four classical **Equations of State (EOS)**. Built as part of a computational materials science exercise.

---

## Overview

In computational materials science, fitting DFT energy–volume data to an EOS is the standard method for extracting equilibrium bulk properties of a material. This project implements and compares four EOS models:

| Model | Type |
|---|---|
| Murnaghan | 1st-order, simple |
| Birch–Murnaghan | 3rd-order finite strain |
| Vinet–Rose | Universal EOS |
| Poirier–Tarantola | Logarithmic strain |

---

## Results Summary

| Model | E₀ (Ry/atom) | V₀ (Å^3/atom) | B₀ (GPa) | B₀′ | R^2 (%) | RMSE (mRy) |
|---|---|---|---|---|---|---|
| Murnaghan | −0.730978 | 18.796 | 165.797 | 11.721 | 85.9989 | 5.5426 |
| Birch–Murnaghan | −0.747576 | 22.800 | 117.343 | 3.439 | 99.9990 | 0.0477 |
| Vinet–Rose | −0.747583 | 22.814 | 87.916 | 3.212 | 99.9998 | 0.0231 |
| **Poirier–Tarantola** | **−0.747584** | **22.815** | **87.925** | **1.197** | **99.9998** | **0.0222** |

**Best Model: Poirier–Tarantola** (lowest RMSE, highest R²)

The fitted V₀ ≈ 22.8 Å^3/atom and B₀ ≈ 88 GPa are consistent with experimental values for hcp Zr (B₀ ≈ 83–92 GPa).

---

## Output Plots

| File | Description |
|---|---|
| `Zr_EOS_fit.png` | E–V curves for all four EOS overlaid on DFT data |
| `Zr_EOS_residuals.png` | Residual bar plots for each model |
| `Zr_EOS_parameters.png` | Bar chart comparison of V₀, B₀, B₀′ |
| `Zr_EOS_quality.png` | R² and RMSE comparison |

---

## Requirements
numpy
scipy
matplotlib

Install with:

```bash
pip install numpy scipy matplotlib
```

---

## Usage

1. Clone the repository:
```bash
   git clone https://github.com/YOUR_USERNAME/EOS-fitting-Zr.git
   cd EOS-fitting-Zr
```

2. Run the script:
```bash
   python EOS-fitting-Zr.py
```

3. Output: four `.png` plots are saved in the same directory, and fitted parameters are printed to the terminal.

---

# Project Structure

```
EOS-fitting-Zr/
├── EOS-fitting-Zr.py       # Main fitting script
├── Zr_EV.dat               # Raw DFT Energy–Volume data (bohr³/Ry units)
├── Summary.txt             # Human-readable results summary
├── Zr_EOS_fit.png          # E–V fit plot
├── Zr_EOS_residuals.png    # Residuals plot
├── Zr_EOS_parameters.png   # Parameter comparison
├── Zr_EOS_quality.png      # Fit quality comparison
|── .gitignore
└── README.md
```


## Physical Notes

- Input data is in atomic units: **bohr^3/atom** (volume) and **Ry/atom** (energy)
- Unit conversions applied internally: 1 bohr^3 = 0.148185 Å^3, 1 Ry/bohr^3 = 14710.5 GPa
- Murnaghan EOS fails here because the volume range spans ~70% compression (111–190 bohr^3/atom), beyond its valid 1st-order approximation
- Vinet–Rose and Poirier–Tarantola are recommended for large-compression datasets

---


