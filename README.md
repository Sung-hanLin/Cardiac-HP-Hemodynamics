# Hemodynamic Modeling and Phase-Resolved Reconstruction of Cardiac HP 13C MRS

This repository provides the MATLAB implementation for processing hyperpolarized [1-13C]pyruvate cardiac data using 1H cine and ECG timing.

## Workflow
### 1. Digital Phantom Simulation
Located in `scripts_simulation/`, these scripts reproduce **Figure 2** and **Supplementary Figure 1** from the manuscript, demonstrating the impact of heart rate variability on signal fidelity.

### 2. In-Vivo Data Fitting & Reconstruction
The core analysis framework is located in `scripts_analysis/`:
- `RUN_reconstruction.m`: The primary entry point. It performs multi-stage fitting to RV, LV, and recirculation components and reconstructs signals at targeted cardiac phases (e.g., end-systole).
- `reproduce_figures_combined.m`: Visualizes the comparison between end-diastolic and end-systolic envelopes (similar to **Figure 7**).

## Key Features
- [cite_start]**3-Component Modeling**: Separates first-pass and recirculation kinetics[cite: 1189].
- [cite_start]**Volume-Guided Correction**: Uses subject-specific 1H cine as a template to correct phase misalignment[cite: 1214].
- [cite_start]**Robust Optimization**: Employs Trust-Region Reflective algorithms for reliable parameter estimation[cite: 1145].

## Citation
Lin SH, et al. "Hemodynamic modeling and phase-resolved reconstruction of hyperpolarized cardiac 13C MRS using 1H cine and ECG timing." (2025).