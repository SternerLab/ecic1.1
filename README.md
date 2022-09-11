[DRAFT IN PROGRESS] 

INTRO 

The ECIC package provides tools for Neyman-Pearson style statistical tests in a generalized setting with three or more models and no designated null model. Non-nested models are also included. See Cullan, Lidgard, and Sterner 2020 for more information: https://www.tandfonline.com/doi/abs/10.1080/02664763.2019.1701636

Instead of a likelihood ratio, ECIC works with a generalized model score function such as AIC, BIC, or an Akaike weights ratio. The procedure implements control over the false positive rate, alpha, by simulating distributions of model score differences and using these to set a numerically-derived threshold for accepting the observed best model or declining to choose a best model due to insufficient evidence. ECIC avoids the use of null models by simulating score differences from each alternative model to the observed best one and using the most stringent threshold derived from these alternatives. 

IMPLEMENTATION

ECIC is intended to be modular in design in order to accommodate commonly used data types, model types, model scores, and model selection decision procedures. Currently implemented features include:

Data and model types: ECIC is able to handle "lm" models from the "lm" R package and "PaleoTS" models and data objects from the PaleoTS R package.

Model scores: AIC, Akaike weights, and BIC scores are implemented.

Decision procedures: Beyond the ECIC procedure itself, the package reports the standard "just pick the best model" procedure by default, along with the Burnham-Anderson rule of thumb using Akaike weights. 

EXAMPLES
See scripts in spanos_test.R, rwalktest.R, paleofit_test.r.

GUIDE TO FILES:

ecic.R: defines the main ECIC function that takes a dataset and list of models and applies a specified score function and decision procedure.

ecic_model.R: provides wrappers for ECIC model objects (and lists of ECIC model objects) of various types (e.g. lm, random walk, and normal distributions). 

ic_scores.R: defines functions for model scores, including the AIC, BIC, and Akaike weight

bias_correct.R: implements a heuristic correction to maximum likelihood estimates of model parameters before simulating from alternative models. The correction accounts for conditioning the simulations on a particular model being observed best.

ecic_control.R: defines the ecicControl function to calculate parametric bootstrap distributions of score differences assuming alternative models to the observed best are true

model_frequences.R: calculate frequency the observed best model also scores best assuming an alternate model is true


