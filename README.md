# Emergence of Criticality MATLAB Code

This repository contains MATLAB scripts used to analyze Curie-Weiss maximum-entropy models fit to neural population statistics. The draft paper currently in this folder is:

`Criticality_Paper__Nature_ (1).pdf`

The code is organized around helper functions in the repository root, dataset-specific preprocessing in the data folders, and figure scripts in `Plots/`.

## Repository Layout

- `Plots/` contains scripts used to generate paper figures and diagnostics.
- `Allen/`, `Hippo/`, `Salamander/`, and `Stringer/` contain dataset-specific parsing, cleaning, and saved data files.
- Root-level `.m` files contain shared Curie-Weiss calculations, exact statistics, Jacobians, and inverse mappings.

## Data Included In Git

The repository tracks the processed data needed by the current plotting scripts:

- `Allen/Allenhldata.mat`
- `Hippo/hippomuchidata.mat`
- `Stringer/stringerhldata.mat`
- `Stringer/stringerhldata2.mat`

Large binary data are stored with Git LFS. Install Git LFS before cloning if you need the `.mat` data files locally. Raw Allen, hippocampus, and Stringer source recordings are not tracked because the plotting scripts use the processed files above and the raw recordings are much larger than the processed plotting data.

## Main Dependencies

Most plotting scripts rely on these root-level helper functions:

- `hlambda.m`: infers Curie-Weiss parameters from `m`, `chi`, and `N`.
- `muChiExact2Spin.m`: evaluates exact finite-size statistics and moments.
- `Jacobian.m`: evaluates the local mapping volume between parameter space and statistics space.

## Figure Script Map

The mapping below follows the current draft paper captions.

- Fig. 1: `Plots/Muchiandhlambdaplots.m`
  Generates empirical points in `(mu, chi)` and inferred model points in `(h, lambda)`. Use the full-data setting (`figureSets = {[-1]}`) for the Fig. 1 panels.

- Fig. 2: `Plots/JacobianCartoon.m`
  Generates the parameter-space circle grid, the mapped statistics-space shapes, susceptibility-vs-field curves, and Jacobian-vs-field curves for `N = 20`.

- Fig. 3a: `Plots/CountourPlots.m`
  Generates Jeffreys-prior contour plots for `N = 20`, `100`, and `1000`. The same script also includes a Jacobian contour diagnostic.

- Fig. 3b: `Plots/RGvectorFieldFixedN.m`
  Generates fixed-system-size vector fields showing how model parameters flow as `N` changes.

- Fig. 3c: `Plots/dotplots.m`
  Maps a rectangular grid in `(mu, chi)` into `(h, lambda)` for increasing system sizes and draws the mapped box boundaries.

- Fig. 4: `Plots/Muchiandhlambdaplots.m`
  Uses the subsampled setting (`figureSets = {[20 100 1000 -1]}`) to show how real-data points move in statistics space and model space as `N` changes.

- Fig. 5b: `Plots/doublewellplot.m`
  Plots the inferred free-energy landscape for a fixed pair of statistics while varying system size.

- Fig. 5c-d: `Plots/lambdamu.m`
  Compares exact finite-size inference against the double-well approximation for selected populations from each dataset.

- Fig. 5 assembled figure: `Plots/Figure5Panels.m`
  Generates the current Figure 5 layout: neural data relative to the critical region, inferred free energy for fixed statistics, separate 3-by-2 parameter panels for `h` and `lambda` versus `N`, and the full-statistic double-well trajectories in `(h, lambda)` space. The script also prints separate comparison/inset figures used for Illustrator assembly.

## Supplementary Figure Scripts

- Curie-Weiss signatures supplement: `Plots/SupplementaryCurieWeissFigures.m`
  Generates the six-panel supplementary Curie-Weiss figure: specific heat versus temperature, specific heat versus interaction scale, Jacobian magnitude versus temperature, Jacobian magnitude versus interaction scale, response to field perturbations with the `~|\delta h^{-1}|` guide, and entropy versus energy.

- Original Curie-Weiss signatures exploration: `Plots/signaturespapercurieweiss.m`
  Retained as the source exploratory script for the Curie-Weiss signature panels and dataset-specific diagnostics.

- Shuffled full-data control: `Plots/ShuffledFullDataFigure1.m`
  Recomputes full-data statistics after independently shuffling each neuron's time series and plots the Figure 1-style statistics and parameter panels.

## Running Notes

Run plotting scripts from the `Plots/` folder so their relative path setup resolves the root helper functions and data folders. Some scripts perform dense exact finite-size calculations and may take time to finish.

The data files loaded by the plotting scripts are expected to already exist in the data folders or on the MATLAB path, for example `Allenhldata.mat`, `hippomuchidata.mat`, and `stringerhldata.mat`.
