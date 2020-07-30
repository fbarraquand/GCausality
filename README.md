[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3967592.svg)](https://doi.org/10.5281/zenodo.3967592)


# GCausality

Companion code to *Inferring species interactions using Granger causality and convergent cross mapping* by F. Barraquand, C. Picoche, M. Detto and F. Hartig. [http://arxiv.org/abs/1909.00731](http://arxiv.org/abs/1909.00731)

Linear Granger causality and convergent cross-mapping are implemented using ``R``. Here, we stick to time-domain approaches from packages ``vars`` and ``lmtest`` for Granger causality, as well as [SIMoNe](https://github.com/cran/simone) for regularized models, and call [rEDM](https://ha0ye.github.io/rEDM/) for convergent cross-mapping (see Hao Ye et al.'s [tutorial](https://ha0ye.github.io/rEDM/articles/rEDM.html) for more information). 

The analyses are organised in folders corresponding to our case studies, with 2 species interacting systems, 2 species and some added abiotic forcing, 10 and finally 20 species networks;  see the Methods of the [arXiv preprint](http://arxiv.org/abs/1909.00731) for a description. 

``Lyapunov`` computes the Lyapunov exponents for all case studies. The READMEs within each folder add more information about their structure and content. 

