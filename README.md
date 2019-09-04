# GCausality

Companion code to *Inferring species interactions using Granger causality and convergent cross mapping* by F. Barraquand, C. Picoche, M. Detto and F. Hartig. 

Linear Granger causality and convergent cross-mapping are implemented using standard ``R`` packages. Here, we stick to time-domain approaches from packages ``vars`` and ``lmtest`` for Granger causality, and call ``rEDM`` for convergent cross-mapping. 

The analyses are organised in folders corresponding to our case studies, with 2 species interacting systems, 2 species and some added abiotic forcing, 10 and finally 20 species networks;  see the Methods of the [arXiv preprint](http://arxiv.org/abs/1909.00731) for a description. 

``Lyapunov`` computes the Lyapunov exponents for all case studies. The READMEs within each folder add more information about their structure and content. 

