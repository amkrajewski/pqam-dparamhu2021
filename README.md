# PyQAlloy-compatible Model for D Parameter prediction based on Hu2021

Release: ![PyPI](https://img.shields.io/pypi/v/pqam-dparamhu2021)

Tests: [![small runtime test](https://github.com/amkrajewski/pqam-dparamhu2021/actions/workflows/runtimeTest.yml/badge.svg)](https://github.com/amkrajewski/pqam-dparamhu2021/actions/workflows/runtimeTest.yml)

License: [![License: LGPL v3](https://img.shields.io/badge/License-LGPL_v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

This repository contains a PyQAlloy-compatible compositionalModel, compatible with [ULTERA Database (ultera.org)](https://ultera.org) infrastructure, for D Parameter Prediction Based on [Yong-Jie Hu 2021 (10.1016/j.actamat.2021.116800)](https://doi.org/10.1016/j.actamat.2021.116800) that 
accepts a chemical formula string of an alloy or a `pymatgen.Composition` object. It return:

Output Order: [`gfse`, `surf`, `dparam`]

Output Meaning (all based on [10.1016/j.actamat.2021.116800](https://doi.org/10.1016/j.actamat.2021.116800)):
- `gfse` - Genralized Stacking Fault Energy (GSF) [J/m^2]
- `surf` - Surface Energy (Surf) [J/m^2]
- `dparam` - D Parameter [unitless] calculated as `surf/gfse`



## Install and use

To run this model you will need **Python 3.9+** and **R 4.1.0+** installed on your system, ideally **before** you install
this software. For Python, we recommend you use a virtual Conda environment, which chan be created with minimal effort 
(see [Miniconda install instructions](https://docs.conda.io/en/latest/miniconda.html)). For R, it can be downloaded 
pre-compiled from a _Comprehensive R Archive Network_ repository (e.g. [Case CRAN](https://cran.case.edu)) and should 
work on most systems, including ARM-based (e.g. Apple M1).

If you have Python and R, you can simply install this model with:
    
    pip install pqam_dparamhu2021

Then, use should be as simple as:

    import pqam_dparamhu2021
    
    print(pqam_dparamhu2021.predict("W30 Mo25 Ta45"))

***

In some cases, required `locfit` R package may not be installed automatically. If you get an error message about it,
try to go to your R console, typically by typing `R` in your terminal, and install it manually with:

    install.packages("locfit")

Or, if you are automating things and need a single-liner, on Mac OS and Linux, the following should work:

    Rscript -e "install.packages('locfit', repos='http://cran.us.r-project.org')"

## Attribution

This repository has been created by Adam M. Krajewski (https://orcid.org/0000-0002-2266-0099) and is licensed under the MIT License. 
**Please cite this repository if you use it in your work.**

The featurizer and predictive model (HEA_pred.R and dependencies) have been optimized across and re-styled by Adam M.
Krajewski based on code originally developed by Young-Jie Hu (https://orcid.org/0000-0003-1500-4015) et al. for their
journal publication and published in Materials Commons at https://doi.org/10.13011/m3-rkg0-zh65, where original code
can be accessed,  distributed under ODC Open Database License (ODbL) v1.0. **Please cite this publication as well:** 
- Yong-Jie Hu, Aditya Sundar, Shigenobu Ogata, Liang Qi, Screening of generalized stacking fault energies, 
surface energies and intrinsic ductile potency of refractory multicomponent alloys, Acta Materialia, 
Volume 210, 2021, 116800, https://doi.org/10.1016/j.actamat.2021.116800

The gbm-locfit package (Gradient Boosting Machine-Locfit: A GBM framework using local regresssion via Locfit) has been 
developed by Materials Project in 2016 and is distributed under the terms of the MIT License. Details can be found in
its code.


## Hu's README File

>Hello, thank you for your interest in our work!
>Here we provide a script written in R language to take an alloy composition of interest and correspondingly predict the GSF energy, surface energy, and the >ductility parameter based on the SL models in our manuscript ( https://doi.org/10.1016/j.actamat.2021.116800)
>To run the script and make predictions, you need to:
>1)	Download the RStudio platform. (https://www.rstudio.com/) ## No worry, it is open access 😊
>2)	Put all the files you downloaded from Materials Commons (basically all our files) into one local folder
>3)	Open the “predict.R” file in RStudio, input the alloy composition there, execute every line there, and the prediction will jump out in the console window >below. 
>
>Please contact qiliang@umich.edu or yh593@drexel.edu if you have any questions. 
>
>P.S.,
>“predict_screen_4nary_all.csv” is the original data for plotting Figure 7&8 in the manuscript. Other figures in the manuscript can be reproduced by the data >listed in the tables in the manuscript.
