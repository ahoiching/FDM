# FDM
[![DOI](https://zenodo.org/badge/494231969.svg)](https://zenodo.org/badge/latestdoi/494231969)

This repository contains the code to perform the Functional Decomposition of Metabolism (FDM) on a given [cobrapy](https://opencobra.github.io/cobrapy/) compatible metabolic model (in json format). The method and its implementation is described in: 

Mori M., Cheng C. et al., _"Functional Decomposition of Metabolism allows a system-level quantification of fluxes and protein allocation towards specific metabolic functions"_, Nature Communications (2023).

Please cite this publication when using this software. The code is maintained by Chuankai Cheng at https://github.com/ahoiching/FDM.

For the flux balance analysis involved in our pipeline, we used GUROBIpy for the convex optimization. As we ran all of our code on [Google Colab environment](https://colab.research.google.com/), it is necessary to obtain a [Gurobi Web License](https://www.gurobi.com/features/web-license-service/) to perform FDM. After obtaining the license, you can run our [example Colab Notebook](https://github.com/ahoiching/FAM/blob/main/FDM_example.ipynb) to run FDM on a test model.
