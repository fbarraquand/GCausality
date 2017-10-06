# GCausality

Granger and other causal approaches like CCM for inferring interactions from time series. Here, we stick to time-domain approaches implemented in ``R``. 
The analyses rely on packages ``vars`` and ``lmtest`` for VAR, as well as ``rEDM`` for convergent cross-mapping. The analyses are organised in folders depending on the number of species and drivers considered. 


* Two-species models include analyses on the Veilleux predator-prey data and some two-species theoretical (deterministic) models of Sugihara et al., as well as stochastic variants

* The Veilleux predator-prey dataset is used to fit a VAR model (which fits well and predicts Granger causality), as well as to demonstrate how CCM performs on VAR1-simulated data (surprinsingly well, given it should not work with what they call linear models...)

* I consider also 2 species model with an external driver. The idea being that we show how we can easily use conditional GC causality in this context (Sugihara et al. 2012 suggested synchronous drivers may be very problematic for GC). 

* To see bigger, a modular 10 species model with noise is included, to see how pairwise GC fares in that case (note: the full VAR model is trickier to fit, at least without a priori knowledge of the causal structure or a regularization technique - and even with that, probably). 

* The 5 species model of Sugihara will be considered later - most relevant to a paper rather than a presentation. 
