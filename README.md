# Emergence of criticality in models of real neurons

MATLAB code accompanying **David P. Carcamo and Christopher W. Lynn, _Emergence of criticality in models of real neurons_ (2026)**. The manuscript is available at [arXiv:2609.09438](https://arxiv.org/abs/2609.09438).

This repository implements the forward and inverse mappings between observable neural-population statistics and Curie–Weiss maximum-entropy models. It includes exact finite-size statistics, parameter inference, Jacobian calculations, parameter flow with increasing population size, and comparisons with the double-well approximation. The figure scripts in [`Plots/`](Plots/) provide the main entry points for reproducing the analyses.

## Repository contents

| Location | Contents |
| --- | --- |
| [`Plots/`](Plots/) | Main and supplementary figure scripts, with plotting and numerical settings defined in each script. |
| [`Allen/`](Allen/) | Processed summaries and preprocessing code for Allen Institute thalamus and visual-cortex recordings. |
| [`Hippo/`](Hippo/) | Processed summaries and preprocessing code for hippocampal recordings. |
| [`Stringer/`](Stringer/) | Processed summaries and preprocessing code for visual-cortex recordings from Stringer et al. |
| Root-level `.m` files | Shared functions for exact statistics, inverse inference, Jacobians, and continuous-size parameter flow. |

## Requirements

The analyses for the accompanying manuscript were run using **MATLAB R2025a**.

| Software | Use |
| --- | --- |
| MATLAB | Numerical calculations and plotting. The scripts use `xscale` and `yscale`, which were [introduced in R2023b](https://www.mathworks.com/help/matlab/ref/xscale.html); R2025a is the reference release used for the paper. |
| Signal Processing Toolbox | Required by `Figure5Panels.m` and `doublewellplot.m` for [`findpeaks`](https://www.mathworks.com/help/signal/ref/findpeaks.html). |
| Parallel Computing Toolbox | Optional for `ShuffledFullDataFigure1.m`. The script attempts parallel processing and falls back to serial execution if a pool cannot be started. Set `useParallelProcessing = false` in the script to request serial execution explicitly. |
| Statistics and Machine Learning Toolbox | Required by spatial preprocessing (`pdist`, `squareform`) and by the Stringer raw-data step in `ShuffledFullDataFigure1.m` ([`zscore`](https://www.mathworks.com/help/stats/zscore.html)). |
| Git LFS | Required to retrieve the processed `.mat` files from GitHub. |

Run scripts with `Plots/` as MATLAB's current folder and the repository root on the MATLAB path. Some standalone scripts construct paths using Windows directory separators; the current path setup in those scripts assumes Windows.

## Installation

Install MATLAB and the toolboxes needed for the analyses of interest. With Git and [Git LFS](https://git-lfs.com/) available, run the following in a terminal:

```sh
git lfs install
git clone https://github.com/dpcarcamo/Emergence-of-Criticality.git
cd Emergence-of-Criticality
git lfs pull
```

For an existing clone, run `git lfs pull` from its root to retrieve any missing LFS objects. The MATLAB source requires no compilation. Installation time depends primarily on the MATLAB installation and data download; the four tracked processed files occupy approximately 0.94 GB in total on disk.

Open MATLAB and set the current folder to the cloned repository root. Then run:

```matlab
addpath(pwd)
cd Plots
```

Keep the directory structure shown above so the scripts can locate the shared functions and data. Most figure scripts clear workspace variables when they start; edit their settings in the script itself before running.

## Processed data

The following three files are the empirical inputs used by the main figure scripts and the supplementary Curie–Weiss signature analysis:

| File | Dataset | Approximate file size |
| --- | --- | ---: |
| [`Allen/Allenhldata.mat`](Allen/Allenhldata.mat) | Allen Institute thalamus and visual cortex | 8.3 MB |
| [`Hippo/hippomuchidata.mat`](Hippo/hippomuchidata.mat) | Hippocampus | 888 MB |
| [`Stringer/stringerhldata.mat`](Stringer/stringerhldata.mat) | Stringer et al. visual cortex | 18.7 MB |

The repository also tracks `Stringer/stringerhldata2.mat` (23.5 MB), the output filename used by `Stringer/dataprocessing.m`. The current main figure scripts, including `Figure5Panels.m`, load **`stringerhldata.mat`**.

The processed files store population statistics and fitted model parameters. The plotting scripts use the following variables:

| Variable | Meaning |
| --- | --- |
| `Data(1,:,:,:)` | Mean population activity, denoted `m` in the manuscript and often `mu` in the code. |
| `Data(2,:,:,:)` | Correlation strength or susceptibility: `N` times the variance of instantaneous population activity. |
| `Data(3,:,:,:)` | Number of neurons, `N`. |
| `hs`, `ls` | Stored inferred external fields `h` and interaction strengths `lambda`. Both `Muchiandhlambdaplots.m` and the data loader in `Figure5Panels.m` require these fields; Figure 5 refits its parameter-comparison curves. |

The remaining `Data` dimensions index spatial subpopulations, population sizes, and recordings; a single-recording dataset may omit the final singleton dimension. The plotting scripts average across the subpopulation dimension and omit missing values. Existing loaders also interpret zero entries as missing placeholders. Activity is represented by spins taking values `-1` and `+1`.

Figures 2 and 3 use model calculations and do not require empirical data. The main empirical figures can be generated from the processed files without rerunning preprocessing. Recomputing the shuffled controls requires the additional inputs described below.

## Minimal example

After the MATLAB path setup above, the following example evaluates a system of 20 independent spins. It requires no empirical data or additional toolboxes:

```matlab
h = 0.2;
lambda = 0;
N = 20;

[m, chi] = muChiExact2Spin(h, lambda, N);
fprintf('m = %.10f, chi = %.10f\n', m, chi)
```

The analytic result at `lambda = 0` is `m = tanh(h)` and `chi = 1 - tanh(h)^2`, giving:

```text
m = 0.1973753202, chi = 0.9610429830
```

This example demonstrates the forward calculation. To generate empirical figure panels, run:

```matlab
Muchiandhlambdaplots
```

With its current defaults, this script opens the full-population statistics and inferred-parameter panels for Figure 1, the population-size grid for Figure 4, and auxiliary inset and marker-size legend figures. It reads the saved summaries and fitted parameters from the three core data files.

## Reproducing the main figures

Run each entry point below from `Plots/` after completing the path setup. Scripts generate MATLAB figures; some also open separate inset, legend, or colorbar figures for manuscript assembly.

| Figure | Entry point | Analysis and outputs |
| --- | --- | --- |
| **1** | [`Muchiandhlambdaplots.m`](Plots/Muchiandhlambdaplots.m) | Full-population mean activity and correlation strength; corresponding stored model parameters; linear-field inset. |
| **2** | [`JacobianCartoon.m`](Plots/JacobianCartoon.m) | Colored regions in parameter space and their images in statistics space for `N = 20, 100, 1000`; Jacobian versus field and population size. |
| **3** | [`Figure3CombinedPanels.m`](Plots/Figure3CombinedPanels.m) | Combined three-row figure: Jacobian contours, parameter-flow vector fields, and inverse images of a fixed statistics grid for `N = 20, 100, 1000`. |
| **4** | [`Muchiandhlambdaplots.m`](Plots/Muchiandhlambdaplots.m) | Statistics and stored inferred parameters for subpopulations near `N = 20, 100, 1000`, displayed in a two-row grid. |
| **5** | [`Figure5Panels.m`](Plots/Figure5Panels.m) | Empirical statistics relative to the mean-field bound; inferred free-energy profiles; exact and double-well parameter comparisons; double-well trajectories toward the thermodynamic limit. |

### Figure settings

In `Muchiandhlambdaplots.m`, the current settings enable both the full-data panels and the population-size grid:

```matlab
figureSets = {[-1]};          % -1 selects the final/full-population point
plotFullDataPanelsTogether = true;
plotNGridFigure = true;
gridNValues = [20 100 1000];
```

`plotNGridLinearParameterFigure` controls the additional linear-field grid. The population-size grid is generated independently of `figureSets`; switching `figureSets` is not necessary to obtain the current Figure 4 layout.

In `Figure5Panels.m`, `fixedStatisticsMu`, `fixedStatisticsChi`, and `freeEnergyNs` define the free-energy example. The settings `parameterPopulationMode` and `selectedPopulationByDataset` control the recordings used for panels c and d. In the default `"single"` mode, a `nan` entry selects a recording randomly from those with the largest number of valid population sizes. The order of entries is **Allen, Stringer, Hippocampus**. The selected indices are printed and stored in the figure cache. To reproduce a particular selection, set the corresponding recording indices explicitly and recompute the cache. Setting `parameterPopulationMode = "all"` includes all available recordings.

### Standalone panels and related analyses

| Script | Purpose |
| --- | --- |
| [`CountourPlots.m`](Plots/CountourPlots.m), Jacobian section | Standalone Jacobian contours. Use `Figure3CombinedPanels.m` for the current manuscript layout. |
| [`RGvectorFieldFixedN.m`](Plots/RGvectorFieldFixedN.m) | Standalone fixed-size parameter-flow panels for Figure 3b. |
| [`dotplots.m`](Plots/dotplots.m) | Standalone mapping of a rectangular statistics grid and its boundary into parameter space for Figure 3c. |
| [`doublewellplot.m`](Plots/doublewellplot.m) | Standalone free-energy illustration with its own fixed-statistics settings. |
| [`lambdamu.m`](Plots/lambdamu.m) | Standalone exact-versus-double-well parameter comparisons for a selected recording from each dataset. |
| [`FullStatsDoubleWellTrajectories.m`](Plots/FullStatsDoubleWellTrajectories.m) | Double-well trajectories at fixed full-recording statistics, with logarithmic and linear field views. |
| [`RGvectorFieldScaledMagnitude.m`](Plots/RGvectorFieldScaledMagnitude.m) | Alternative vector-field visualization with arrow lengths scaled by the magnitude in plotted coordinates. |

The standalone scripts retain their own settings; use the main entry points above for the combined figure layouts.

## Supplementary analyses

| Script | Analysis |
| --- | --- |
| [`SupplementaryCurieWeissFigures.m`](Plots/SupplementaryCurieWeissFigures.m) | Six-panel analysis of specific heat and Jacobian versus temperature or interaction scale, field-perturbation response, and entropy versus energy. Uses the core processed data. |
| [`RGtrajectories.m`](Plots/RGtrajectories.m) | Integrates continuous-size parameter-flow trajectories from a grid of initial statistics. Initial conditions, population-size range, ODE tolerances, and optional invariant checks are set in the script. Requires no empirical data. |
| [`ShuffledFullDataFigure1.m`](Plots/ShuffledFullDataFigure1.m) | Full-population statistics and inferred parameters after independently shuffling each neuron's time series, together with independent-neuron reference panels. Requires raw inputs when its cache is absent. |
| [`signaturespapercurieweiss.m`](Plots/signaturespapercurieweiss.m) | Extended Curie–Weiss signature analysis with empirical or manually specified targets. |
| [`plot_CW_criticality_signatures.m`](Plots/plot_CW_criticality_signatures.m) | Wrapper that runs `signaturespapercurieweiss.m`. |

### Inputs for shuffled controls

On a fresh run, `ShuffledFullDataFigure1.m` reads:

| Location | Required input |
| --- | --- |
| `Allen/Allen_data_share/` | Raw `.mat` recordings containing the binary neuron-by-time matrix `X`. |
| `Hippo/Binary_1485neurons_30Hz_share.csv` | Sparse hippocampal activity events, with neuron and time indices in the first two columns. |
| `Stringer/data/` | Recording `.mat` files containing the `stim.resp` or `stim.spont` response arrays. |
| `Stringer/data_info.csv` | Recording metadata with `file_name` and `stimulus` columns selecting the file and response type. |

These raw recordings and metadata are not included in the GitHub checkout. The original experiments are described and cited in the manuscript. The three processed summaries used by the main figures do not contain the individual time series needed for reshuffling.

The script caches shuffled statistics in `Plots/ShuffledFullDataFigure1_cache.mat` and creates separate shuffled and independent-reference figures. Time-series shuffling is stochastic. Preserve the generated statistics and random-stream configuration when comparing a particular realization; serial and parallel execution require appropriate random-stream control for matching realizations.

## Computation, caching, and outputs

Exact finite-size calculations sum over the `N + 1` possible population activities. Full-resolution calculations can therefore be substantially more expensive than the minimal example. In particular, Figure 3 uses a `1000 × 1000` contour grid at each of three population sizes and also computes flow fields and inverse mappings. Runtime depends on grid resolution, population size, inference settings, and whether cached values are available. Measured runtime and peak-memory benchmarks are not currently provided.

Several scripts save intermediate `.mat` caches in `Plots/` and reuse compatible cached results on subsequent runs. `Figure5Panels.m` also saves partial results for its empirical parameter calculations. Use the script's `forceRecompute` setting to regenerate values; `Figure3CombinedPanels.m` has separate switches for its contour, vector-field, and dot-plot calculations. After changing source data or shared numerical routines, recompute affected caches rather than relying on unchanged filenames or grid settings. Cache files are excluded from Git.

The scripts primarily create figures in the MATLAB session. Save or export the required figure windows using MATLAB's figure tools; a full set of manuscript PDFs is not exported automatically. Input summaries can expand substantially in memory when loaded, so memory requirements depend on the arrays and analysis, as well as the number of parallel workers.

## Shared numerical functions

| Function | Role |
| --- | --- |
| [`muChiExact2Spin.m`](muChiExact2Spin.m) | Exact finite-size moments, mean activity, susceptibility, and energy variance. |
| [`muchi.m`](muchi.m) | Forward-statistics interface used during parameter inference. |
| [`hlambda.m`](hlambda.m) | Numerical inversion from mean activity and susceptibility to `h` and `lambda` at fixed `N`. |
| [`Jacobian.m`](Jacobian.m) | Jacobian of the mapping from `(h, lambda)` to `(m, chi)` and its determinant. |
| [`DhDNtot.m`](DhDNtot.m), [`DlamDNtot.m`](DlamDNtot.m) | Parameter flow with changing population size, using the root-level integral-statistics and derivative functions. |

## Citation

Please cite the accompanying manuscript when using this code:

```bibtex
@article{carcamo2026emergence,
  title         = {Emergence of criticality in models of real neurons},
  author        = {Carcamo, David P. and Lynn, Christopher W.},
  year          = {2026},
  eprint        = {2609.09438},
  archivePrefix = {arXiv},
  primaryClass  = {physics.bio-ph},
  url           = {https://arxiv.org/abs/2609.09438}
}
```

Please also cite the original experimental datasets when using their data, as detailed in the manuscript. For reproducibility, record the repository commit and any modified script settings used for an analysis. Correspondence information is provided in the accompanying paper.
