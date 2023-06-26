



###simulation code

###result1: (Figure 3)
#The three X-M model Scenarios with different Strength of Modification (SoM) using RFQT and single-stratification


ATE_MSE<-c()

Gammas<-c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    )#共6中SoM情况
for(i in 1:length(Gammas)){
  set.seed(10*i)
  gamma_used<-Gammas[i]
  Dat<-getDat(scenario='A', SoM=gamma_used)
  odat<-Dat$traning.set  #training set in 全局环境
  vdat<-Dat$testing.set  #testing set in 全局环境

  ##############ATE
  odat#用training set估计
  fitGX<-lm(    odat[,3]~  odat[,2]  )  ; bx<-as.numeric( summary(fitGX)$coef[-1,1]  ); bxse<-as.numeric(  summary(fitGX)$coef[-1,2])
  fitGY<-lm(    odat[,4]~  odat[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
  MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
  ATE_MSE<-c(   ATE_MSE     ,        mean(  (vdat$true_STE - MRres@Estimate )^2 ) )


  ###############Single Stratification
  ###DR
  ALLRES<-RFQTfit(odat,vdat,Nb=200,SingleM=TRUE,)
  DR_S_MSE<-c(  DR_S_MSE ,  SingleMTreeFitting() )

  ###R


  ###N


}





###result2: (Figure 7)
#MSE and individual predicted effect curve changing along the number ofQ trees









###result3: (Figure S1)
#Permutation test statistics
















