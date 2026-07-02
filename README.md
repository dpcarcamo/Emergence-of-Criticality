# Emergence of Criticality MATLAB Code

This repository contains MATLAB scripts used to analyze Curie-Weiss maximum-entropy models fit to neural population statistics. The draft paper currently in this folder is:

`Criticality_Paper__Nature_ (1).pdf`

The code is organized around helper functions in the repository root, dataset-specific preprocessing in the data folders, and figure scripts in `Plots/`.

## Repository Layout

- `Plots/` contains scripts used to generate paper figures and diagnostics.
- `Allen/`, `Hippo/`, `Salamander/`, and `Stringer/` contain dataset-specific parsing, cleaning, and saved data files.
- Root-level `.m` files contain shared Curie-Weiss calculations, exact statistics, Jacobians, and inverse mappings.

## Main Dependencies

Most plotting scripts rely on these root-level helper functions:

- `hlambda.m`: infers Curie-Weiss parameters from `mu`, `chi`, and `N`.
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

## Running Notes

Run plotting scripts from the `Plots/` folder so their relative path setup resolves the root helper functions and data folders. Some scripts perform dense exact finite-size calculations and may take time to finish.

The data files loaded by the plotting scripts are expected to already exist in the data folders or on the MATLAB path, for example `Allenhldata.mat`, `hippomuchidata.mat`, and `stringerhldata.mat`.
