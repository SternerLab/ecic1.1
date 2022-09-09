[DRAFT IN PROGRESS] 

INTRO 

The ECIC package provides tools for Neyman-Pearson style statistical tests in a generalized setting with three or more models and no designated null model. Non-nested models are also included. See Cullan, Lidgard, and Sterner 2020 for more information: https://www.tandfonline.com/doi/abs/10.1080/02664763.2019.1701636

Instead of a likelihood ratio, ECIC works with a generalized model score function such as AIC, BIC, or an Akaike weights ratio. The procedure implements control over the false positive rate, alpha, by simulating distributions of model score differences and using these to set a numerically-derived threshold for accepting the observed best model or declining to choose a best model due to insufficient evidence. ECIC avoids the use of null models by simulating score differences from each alternative model to the bserved best one and using the most stringent threshold derived from these alternatives.  
