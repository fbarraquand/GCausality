# File description

This folder contains the data, scripts, results and figures for 2 species models (i.e., based on Veilleux data sets, deterministic chaos model from Sugihara and stochastic model with no driver). Some of the files (scripts, figures or results) are in explo folders, corresponding to previous analyses which are not used in the final paper. 

An example of the time series for all models is plotted in `time_series_small.r`.

### Veilleux data sets

*  `Veilleux_lag.R` computes and plots the lag order for the two datasets.
*  `predatorPrey_Veilleux.R` computes and writes Granger-causality testing on the Veilleux datasets and plot the CCM curves corresponding to the Veilleux dataset.
*  `Veilleux_CCM_varp.R` computes and plots the results for CCM using both Veilleux data sets and simulations based on the dynamics observed in the datasets.
*  `Veilleux_CCM_appendices.R` computes and plots the CCM results using both raw and logged-values

### Deterministic chaos model

*  `theoreticalSugiharaModels.R` looks for the best lag for GC analyses on the chaotic model
*  `chaosSugihara.R` plots the results for GC and CCM
*  `chaosSugihara_diff_pval.R` writes down the pvalues for both GC and CCM analyses

### Stochastic model

* `StochasticCompet.R` does the GC and CC analyses on the stochastic model

All results from `chaosSugihara_diff_pval.R` and `StochasticCompet.R` are analyzed and figures are plotted with `results_StochasticCompet.R` 
