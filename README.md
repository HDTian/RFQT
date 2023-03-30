# RFQT
This text illustrates the key messages to draw the **Random Forest of Q trees (RFQT)** and its relevant analysis for treatment heterogeneity studies with Mendelian randomization, it can achieve a lots of work, including
- Draw heterogeneity Q test 
- Build single Q tree  
- Build random forest of Q trees 
- Make treatment effect prediction
- Evaluate Variable Importance (VI) 
- Draw permutation test
- Draw marginal heterogenous effect analysis for any covariates  


When drawing your treatment heterogeneity study, define the one-dimensional instrument *Z*, the exposure of interest *X*, the (high-dimensional) covariates *M*, and the outcome of interest *Y*. If the instrument is high-dimensional, first build the single instrument (e.g. by weighted gene score). The standard data structure regonized by RFQT is *>cbind(individual index, Z , X , Y, M)*.

Below is the mind line to understand what RFQT did in each step and how the internal functions work, which is particularly helpful to those working on methodology improvement. If you have a data and wish to implement RFQT to obtain the results with a step-by-step guidance, you may wish to see the final summary section.

## Model and toy example
Below is one heterogenous treatment example used for simulation considering 20 candidiate covariates where the first five covariates are the treatment effect modifiors  
$$X = 0.5Z + 0.5 \sum_{j=1}^{20} U_j  + \epsilon_X $$ 
$$M_j  =  b_j X + U_j \qquad j=1,2,\ldots,20 $$
$$Y = \left(   0.5+ \sum_{j=1}^{5} \gamma_j M_j  \right) X +  0.5 \sum_{j=1}^{20} U_j  + \epsilon_Y  $$
where $b_j$ is the effect of the exposure on the covariate $M_j$; when $b_j$ is not zero, it causes the collider bias conditioning on $M_j$. $\gamma_j$ is the strength of modification for $M_j$. The command `getDat()` helps to build standard toy samples based on the euqations above.  

One need to evaluate the heterogenous treatment effects based on different levels of effect modifiors. RFQT is a collider robust methods combing with data-adaptive methods that can let you estimate reliable individual treatment effects, find out the potential effect modifors from high-dimensional covariates, and many other relevant analysis.
## Build a single Q tree
Constructing a single Q tree caintains two iterative steps: 
1. Determine the splitting covariate (achieved by `GetIndex`)
2. Split based on the candidate covariate
The schematic diagram is
![Q-tree1024_1](https://user-images.githubusercontent.com/127906571/225363783-32754381-27a3-45aa-9591-c2ea56bfd89b.jpg)
Let `dat` be the standard data, simply  running `GetTree(dat)` to build a tree. It will return a dataset containing all individual information necessary for further analysis.
## Build RFQT
One may need build RFQT to reduce variance of any estimator based on a single Q tree. RFQT is a multiple-fitting of Q tree with bootstrap data. The function `BootstrapTreeFitting` or `BootstrapTreeFitting_real` helps to bootstrap samples and return a list containing neccessary results, including testing/OOB individual effect preditons and VI measurements for the present bootstrap. To run RFQT, one need bootstrap many times and this may be time-consuming. One is suggested to use parallel computation like



       cl<-makeCluster(n_cores_used)
       clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat','GetTree', 'GetNindex', 'GetIndex' )  )
       RES<-parSapply(   cl ,  1:100, BootstrapTreeFitting  ) ##RES are the RFQT fitting results

`RES` contains individual results for all bootstrap Q trees, and it will be used as the further analysis.

## Further analysis
One may wish to obtain the MSE (if one is working on simulation data and methodology improvement), this is easily achieved by `getMSE(RES,2)`. To get individual predicted effect of RFQT, use `getPredict(RES,2)`. VI measurements can be derived runing `getVI2(RES,'order')`. There are many other analysis metrics, users can easily define and build their preferred metrics based on `RES`.   

## Permutation test
The permutation test is to assess the effect heterogeneity, accounting for uncertainty from the RFQT algorithm. The null hypothesis is that all the candidate covariates considered do not modify the treatment effect, and so the individual-level predicted causal effects are not more variable than would be expected due to chance alone. The test is a multiple fitting version of RFQT, each time the covariate information is permuted. The relevant RFQT fitting function under the permutation test is `BootstrapTreeFitting_real_bp`.


## Step-by-Step guidance
Below is the step-by-step no-lost guidance for those implementing RFQT with a data.

Depending on the type of data (simulation data or real application data) one wish to analyze, follow the coresponding guidance.

### 0. Load the RFQT package from GitHub
       library(MendelianRandomization)
       library(tidyverse)
       library(data.table) 
       library(parallel) 
       library(devtools)
       install_github("HDTian/RFQT")
       library(RFQT)
       

### 1. Build the standard data
Prepare the standard data `Dat` where the first four colunms are individual index, instrument values, the exposure, the outcome, following the high-dimensional covariate information. If missing data exists, the individual will be removed. If one wish to investigate potential non-linear exposure effect patterns, add the exposure itself as a covariate. 

For simulated data, the data strucuture is the same as the real data, only with the true individual effects added as the final column. One can use the command `getDat()` to obtain standard toy samples.

If one wish to left a testing sub-set for validation and seperate analysis, simply split `Dat` as the training set `odat` and the testing set `vdat`.

### 2. Define any hyperparameters
Assign the hyperparameter values as in function arguments (i.e. in the local environment), for example

       GetTree(dat, S = 5, endsize = 5000, rate = 0.4, Qthreshold = 3,...) 
       
The hyperparameters `S` `endsize` `rate` `Qthreshold` refer to the maximun depth of Q-tree, the minimal size of end node, the proportion of covariates considered in each split and the threshold value for Q statistic, respectively.

### 3. Run RFQT
With all hyperparameters defined and all data structures needed, run the following codes to fit a RFQT (you can use other cluster commands as you prefer)

       library(parallel)
       n_cores_used<-detectCores()-1
       cl<-makeCluster(n_cores_used)
       clusterEvalQ(cl=cl , expr=library(dplyr))
       clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
       clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat','GetTree', 'GetNindex', 'GetIndex' )  )
       Nb<-n_cores_used
       RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  ) 
       
       ## Or if you wish to use user-defined parameters like SoP=20, try
       user_BootstrapTreeFitting<-function(seed){  RES<-BootstrapTreeFitting(seed,SoP=20);  return(RES)  }
       clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat', 'GetTree', 'GetNindex', 'GetIndex' , 'BootstrapTreeFitting')  )
       RES<-parSapply(   cl ,  1:Nb, user_BootstrapTreeFitting  ) 

       stopCluster(cl)
Here we use the doubly-ranked stratification method with the pre-stratum size 10 (defalut values). Other stratification methods including the residual method (`method<-'Residual'`) and the naive method is allowed. 

### 4. Obtain the results you wish
`RES` contains lots of useful information that can be transformed to the metrics one need. Here are some scenarios:

Obtain indidividual predicts histogram

       predict_matrix<-getPredict( RES,2     )
       predict_RFQT<-predict_matrix[, ncol( predict_matrix )]
       hist( predict_RFQT ,100, xlab='Predicted effect')
       
Obtain MSE (only applicable if the true individual effects are known and inputed)

       MSE1<-getMSE( RES , 1 )  #MSE1: OOB error
       MSE2<-getMSE( RES , 2 )  #MSE2: test error
       
Obtain VI meaasurements for all the candidiatw covariates 

       getVI1(RES,'order')  
       getVI2(RES,'order')  #for real data

One can do other analysis, like draw the permutation test

       permutation_test(Dat,100) #100 is the number of permutation
       
or get the marginal covariate plot 

       marginal_covariate_MR( Dat,index,method ) #index is the covariate index one choose; method refers the stratification used

### 5. Integrated commands
One may wish to use an intergrated command to obtain certain results with the data inputed. One could use

       RFQTfit(odat,vdat,...)
       
where `odat` and `vdat` are generally the training set used for RFQT fitting and the testing set used for seperate analysis, divided from `Dat`. This integrated command returns the basic results (`RES`, predictions, VI measurements, etc). Be careful that this command could be time-costing.
