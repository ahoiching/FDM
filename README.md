# FDM
This repository contains the code to perform the Functional Decomposition of Metabolism (FDM) on a given [cobrapy](https://opencobra.github.io/cobrapy/) compatible metabolic model (in json format). The method is described in: 

Mori M., Cheng C. et al., _"Functional Decomposition of Metabolism allows a system-level quantification of fluxes and protein allocation towards specific metabolic functions"_, Nature Communications (2023)

For the flux balance analysis involved in our pipeline, we used GUROBIpy for the convex optimization. As we ran all of our code on [Google Colab environment](https://colab.research.google.com/). I used the [Web License Service](https://www.gurobi.com/features/web-license-service/) for applying my GUROBI license.

To get started, see the [Colab Notebook](https://github.com/ahoiching/FAM/blob/main/FDM_example.ipynb).
