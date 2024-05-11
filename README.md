# Random Forest of Q Trees
[A data-adaptive method for investigating effect heterogeneity with high-dimensional covariates in Mendelian randomization](https://link.springer.com/article/10.1186/s12874-024-02153-1)

The data-adaptive effect heterogeneity analysis in Mendelian randomization (MR) with high-dimensional covariates and collider-robust techniques. 

Random Forest of Q Tree (RFQT) can estimate individual or subgroup heterogeneous causal effects and identify the strong effect modifications. Additionally, the package aids in conducting effect heterogeneity testing and constructing the decision structure for predicting individual effects.

RFQT is a data-adaptive method, which means there is no need to pre-specify effect modifications, and it accommodates high-dimensional covariates. It also permits covariates to be downstream colliders of the exposure and is able to avoid the collider bias.

This manuscript provides a simple, step-by-step guide to conducting effect heterogeneity analysis in MR using the RFQT package.

## Load the Package
Install the RFQT package from Github with:
```R
devtools::install_github("HDTian/RFQT")
```
```R
library(RFQT)
```


## Data Preparation
Before fitting, it is advisable to perform some basic data checks. Assume you have the data named as `Dat`. Make sure: 

(i) Have the data as a data frame
```R
is.data.frame(Dat) 
```
```R
#TRUE
```

(ii) The column variables, in order, are:: ID, the (one-dimensional[^1]) instrument, the exposure of interest, the outcome of interest, and (possibly high-dimensional) covariates.

(iii) Rename the column variables of the data frame in order: `I`, `Z`, `X`, `Y`, `...`. The name of the covariates can be arbitrary.
```R
colnames(Dat)
```
```R
#"I"     "Z"     "X"     "Y"     "ages"  "sbp"   "dbp"  ... 
```

(iv) Check all variables are numeric[^2]
```R
apply(Dat, 2, is.numeric)
```
```R
#   I        Z        X        Y     ages      sbp      dbp  ...
#TRUE     TRUE     TRUE     TRUE     TRUE     TRUE     TRUE  ...
```



[^1]: When the instruments are high-dimensional, such as multiple genetic variants in MR, first convert them into a single weighted gene score 'Z' and then store them in the data frame 'Dat'.
[^2]: Although our method supports discrete or coarsened variables, we encourage using continuous variables for RFQT fitting to ensure more reliable results. For binary variables like gender and ancestry, you can perform variable selection before processing 'Dat'.

## RFQT Fitting
The most important function used for RFQT fitting is `RFQTfit()`, which offers multiple tuning parameters/arguments (try `?RFQTfit` for details). The default settings are suitable for most general analyses. Assuming you have split the original data `Dat` into the training set `odat` and testing set `vdat`; Note that the testing set doesn't necessarily have to be derived from splitting `Dat`. Alternatively, you can set `odat<-Dat` and construct `vdat` with user-defined values or simply use `vdat<-odat`[^3].

Now run
```R
ALLRES_real<-RFQTfit(odat,vdat,Nb=200)
```
where `Nb=200` indicates the forest contains 200 Q-trees. `ALLRES_real` contains the fitting results. Note that `RFQTfit` has a built-in multi-core parallel setup and one fit of RFQT with 200 trees and 8 cores may cost half hour or longer. For multiple fitting task like simulation or permutation, it is stongly suggested to use HPC.  

[^3]: If you copy the training data to testing data, note that unlike traditional ML models, the same individual in the tratining set and testing set may have different estimates/predictions in RFQT (the difference should be small and often exactly zero).

## Results Analysis
You may have various goals in mind when obtaining or presenting results. Some of these objectives can be readily achieved by directly calling `ALLRES_real`; for example

*Plot the individual predicted heterogenous effects
```R
hist( ALLRES_real$Predicts_test[,200] , n=200, main='', xlab='Predicted effect')
```
*Obtain the variable importance
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

*Re-stratification: modeling the individual RFQT predicted effect alongside individual covariates to gain deeper insights into the modification pattern. Various modeling techniques, such as simple linear regression and decision trees, can be employed for this purpose:
```R
fit  <- lm(  ALLRES_real$Predicts_test[,200] ~ as.matrix(vdat[, 5:ncol(vdat) ])     )#simple linear regression
```
```R
tree <- rpart::rpart(estHTE ~. , data = data.frame(M=vdat[,5:ncol(vdat)],estHTE=ALLRES_real$Predicts_test[,200])) #decision tree
```

*Single variable MR fitting: only one variable, such as the 10th column variable of 'Dat,' typically indicated by the variable importance results of RFQT, for further stratification MR analysis. This can be accomplished using the DRMR package (with quantile/split points determined by the user) or by tree fitting (with splitting points determined by the data)
```R
devtools::install_github("HDTian/DRMR");library(DRMR)
single_dat <- Dat[, c(2:4,10)] ; colnames(single_dat)[4]<-'M'
DRMRfit <- getSummaryInf( Stratify(single_dat,onExposure = FALSE)     ) #via DRMR code
```
```R
single_dat <- Dat[, c(1:4,10)]
DTfit(single_dat) #via tree fitting
```

## Permutation Heterogeneity Test
You may wish to test whether the heterogenous estimates are more variable than expected due to chance alone. This can be acheived by permutation test: first permute the original data and then calculate the test staitstic for each permuated data (stored in `$ts1` or `$ts2` of the RFQT result), and compare them with the real data test staitstic result (i.e. `ALLRES_real$ts1` or `ALLRES_real$ts2`). For more details see the script *illustration_sim_real/Real.R* 

You are encouraged to use HPC as permutation test with RFQT is usually quite time-consuming.

## Special Applications
*Non-linear analysis: if you are doing non-linear MR studies, you can easily fit your model via RFQT (or single tree fitting) as nonlinearity is a particular scneario of heterogenity where the exposure itself act as the effect modifior. That is, 
```R
nonlinear_Dat <- data.frame(I=Dat$I, Z=Dat$Z, X=Dat$X, Y=Dat$Y, M1=Dat$X) #other covariates M2, M3, ... may follow
RFQT_nonlinear <- RFQTfit(odat=nonlinear_Dat,vdat=nonlinear_Dat,Nb=200,rate=1) #RFQT fitting
DT_nonlinear   <- DTfit(odat=nonlinear_Dat,vdat=nonlinear_Dat,rate=1)   #single tree fitting       
```
the further analysis is similar as above. Note that for single tree fitting, you may wish to use the function `GetTree`, which contains more details for the tree nodes and structure information (try `?GetTree` for explanations).

