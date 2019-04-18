### List of things to do when re-working on the article

* Compare the evaluation by p-values, effect sizes and combination for the two species models. Only after that, work on the results for the 10- and 20-species models

* Include - in Appendix or additional figure - (stochastic) Lyapunov exponents computed for the models

* Something to document additionnally: the construction of surrogates for CCM (and perhaps to get better p-values for GC?)







### Stuff that might be better in another, follow-up article

* New p-values based on good surrogates, including for GC-type models

* A detailed investigation of confounding with different data structures (currently we only have X<-U->Y):
  1. U->X->Y
  2. X->Y<-U
  3. U->X-|Y
This could be combined with a non-dynamic causal evaluation à la Pearl. 
Note: our many-species case actually include confounding to some degree, we can use it to evaluate the causal inference on some of these causal structures

