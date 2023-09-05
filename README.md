# Random Forest of Q Trees
The data-adaptive effect heterogeneity analysis in Mendelian randomization (MR) with high-dimensional covariates and collider-robust techniques. 

Random forest of Q tree (RFQT) can estimate the individual or subgroups heterogenous causal effects and find the most important effect modifications. The package also helps to draw effect heterogenity testing and build the decision structure forindicidual effect prediction. 

RFQT is a data-adaptive method, which means one do not need to pre-specify the effect modification, and allows the covariate to be high-dimensional. It also allows the covariates to be the colliders (downstream effect of the exposure) and to most extent avoid the collider bias.

This manuscript gives the simple step-by-step guide to draw the effect heterogenity anlysis in MR with RFQT package.

## Load the Package
Install the RFQT package from Github with:
```R
devtools::install_github("HDTian/RFQT")
```
```R
library(RFQT)
```


## Data Preparation
Before fitting, it is better to do some simple data check. Assume you have the data named as `Dat`. Make sure: 

(i) you have the data as a data frame
```R
is.data.frame(Dat) 
```
```R
#TRUE
```

(ii) the column variables in order are: ID, the (one-dimensional[^1]) instrument, the exposure of interest, the outcome of interest, and (possibly high-dimensional) covariates.

(iii) rename the column variables of the data frame in order: `I`, `Z`, `X`, `Y`, `...`. The name of the covariate can be arbitrary.
```R
colnames(Dat)
```
```R
#"I"     "Z"     "X"     "Y"     "ages"  "sbp"   "dbp"  ... 
```

(iv) check all variables are numeric variables[^2]
```R
apply(Dat, 2, is.numeric)
```
```R
#   I        Z        X        Y     ages      sbp      dbp  ...
#TRUE     TRUE     TRUE     TRUE     TRUE     TRUE     TRUE  ...
```



[^1]: When the instruments are high-dimensional like the multiple genetic variants in MR, first convet them in a single weighted gene score `Z` and them store into our data `Dat`. 
[^2]: Through our method support discrete or coarsened variables, we still encourage continuous variable in RFQT fitting for more reliable results. One can do the selection for binary variables like gender and ancestry in begining (i.e. prior to `Dat`)  

## RFQT Fitting
The most important function used for RFQT fitting is `RFQTfit()`, where multiple tuning parameters exist (try `?RFQTfit` for details) but the default settings are ok for most general anlysis. Here assume you split the original data `Dat` into the traning set `odat` and testing set `vdat` (the testing set is usually unneccessary for real application, where simply let `odat<-Dat` and ignore `vdat`). Now run
```R
ALLRES_real<-RFQTfit(odat,vdat,Nb=200)
```
where `Nb=200` indicates the forest contains 200 Q-trees. `ALLRES_real` contains the fitting results.

## Results Analysis
You may wish to obtain or present different kinds of results. Some of them can be easily achieved by calling`ALLRES_real` directly; for example

*Get the individual predicted heterogenous effects

*Get the variable importance

*Get any single tree end node information

## Results Analysis - Extensions



