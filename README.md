# Random Forest of Q Trees
The data-adaptive effect heterogeneity analysis in Mendelian randomization (MR) with high-dimensional covariates and collider-robust techniques. 

Random forest of Q tree (RFQT) can estimate the individual or subgroups heterogenous causal effects and find the most important effect modifications. The package also helps to draw effect heterogenity testing and build the decision structure for individual effect prediction. 

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
[^2]: Though our method support discrete or coarsened variables, we still encourage continuous variable in RFQT fitting for more reliable results. One can do the selection for binary variables like gender and ancestry in begining (i.e. prior to `Dat`)  

## RFQT Fitting
The most important function used for RFQT fitting is `RFQTfit()`, where multiple tuning parameters.arguments exist (try `?RFQTfit` for details) but the default settings are ok for most general anlysis. Here assume you split the original data `Dat` into the traning set `odat` and testing set `vdat`. The testing set may not needed to be formed by splitting `Dat`, and one can instead let `odat<-Dat`, and build `vdat` with any user-defined values or simply let `vdat<-odat`[^3]. Now run
```R
ALLRES_real<-RFQTfit(odat,vdat,Nb=200)
```
where `Nb=200` indicates the forest contains 200 Q-trees. `ALLRES_real` contains the fitting results. Note that `RFQTfit` has a built-in multi-core parallel setup and one fit of RFQT with 200 trees and 8 cores may cost half hour or longer. For multiple fitting task like simulation or permutation, it is stongly suggest to use HPC.  

[^3]: If you copy the training data to testing data, note that unlike traditional ML models, the same individual in the tratining set and testing set may have different estimates/predicts in RFQT (the difference should be small and most times exactly zero).

## Results Analysis
You may wish to obtain or present different kinds of results. Some of them can be easily achieved by calling`ALLRES_real` directly; for example

*Plot the individual predicted heterogenous effects
```R
hist( ALLRES_real$Predicts_test[,200] , n=200, main='', xlab='Predicted effect')
```
*Get the variable importance
```R
ALLRES_real$VI2
```
*Draw the traceplot for any (eg the 123th) individual to check the forest stabilization
```R
plot( 1:200, ALLRES_real$Predicts_test[123,],  type='l',xlab='Number of Q trees', ylab='Predicted effect')
```
*Get any single (e.g. the 123th) tree end node information
```R
ALLRES_real$RES[1,123]$end_node_information
```


## Results Analysis - Extensions
There are many other results or futher analysis; for example

*Re-stratification: modelling the individual RFQT predicted effect and the individual covariates to have more insights on the modification pattern. There are many models like the simple linear regression and decision tree
```R
fit  <- lm(  ALLRES_real$Predicts_test[,200] ~ as.matrix(vdat[, 5:ncol(vdat) ])     )#simple linear regression
```
```R
tree <- rpart::rpart(estHTE ~. , data = data.frame(M=vdat[,5:ncol(vdat)],estHTE=ALLRES_real$Predicts_test[,200])) #decision tree
```

*Single variable MR fitting: use only one variable (e.g. the 10-th column variable of `Dat`), which was usually indicated by the variable importance result of RFQT, for further stratification MR analysis. This can be achieved by DRMR package (with quantile/split point detarmined by user) or by tree fitting (splitting point driven by data)
```R
devtools::install_github("HDTian/DRMR");library(DRMR)
single_dat <- Dat[, c(2:4,10)] ; colnames(single_dat)[4]<-'M'
DRMRfit <- getSummaryInf( Stratify(single_dat,onExposure = FALSE)     ) #via DRMR code
```
```R
single_dat <- Dat[, c(1:4,10)]
DTfit(single_dat) #via tree fitting
```

## Permutation Heterogenity Test
You may wish to test whether the heterogenous estimates are more variable than expected due to chance alone. This can be acheived by permutation test: first permute the original data and then calculate the test staitstic for each permuated data (stored in `$VI1` or `$VI2`), and compare them with the real data test staitstic result (i.e. `ALLRES_real$VI1` or `ALLRES_real$VI2`). For more details see the script *illustration_sim_real/Real.R* 

You are encouraged to use HPC as permutation test with RFQT is usually quite time-consuming.

## Special Applications
*Non-linear analysis: if you are doing non-linear MR studies, you can easily fit your model via RFQT (or single tree fitting) as nonlinearity is a particular scneario of heterogenity where the exposure itself act as the effect modifior. That is, 
```R
nonlinear_Dat <- data.frame(I=Dat$I, Z=Dat$Z, X=Dat$X, Y=Dat$Y, M1=Dat$X) #other covariates M2, M3, ... may follow
RFQT_nonlinear <- RFQTfit(odat=nonlinear_Dat,vdat=nonlinear_Dat,Nb=200,rate=1) #RFQT fitting
DT_nonlinear   <- DTfit(odat=nonlinear_Dat,vdat=nonlinear_Dat,rate=1)   #single tree fitting       
```
the further analysis is similar as above. Note that for single tree fitting, you may wish to use the function `GetTree`, which contains more details for the tree nodes information (try `?GetTree` for explanations).

