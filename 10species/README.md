# File description

This folder contains the data, scripts, results and figures for high dimension models. Some of the files (scripts, figures or results) are in explo folders, corresponding to previous analyses which are not used in the final paper. Some of the scripts are used to plot figures for both 10 and 20 species models for the sake of clarity in the article. 

In the `script` folder,  

* `10Species_dataGeneration_ref.R` and `10Species_dataGeneration_random.R` produce the simulated dynamics of the 10 species model that are then stored in the `data` folder

* `lag_order_large_simulation.R` computes and represents the best lag for the 10 and 20 species models.

* `analysis_CCM_10species.R` computes the probability associated to the interactions between each species using the CCM method in the 10 species model. 

* `large_simulation_output_per_interaction.R` computes the probability associated to the interactions between each species using the GC method in the 10 and 20 species models.

* `results_large_per_inter_panels.R` produces a figure showing the probability to detect each interaction using the GC or CCM method in the 10 and 20 species models.

* `results_large_ROC_panels.R` produces a figure showing the ROC scatter plots using the GC or CCM method in the 10 and 20 species models.

