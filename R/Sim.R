



###simulation code

###result1: (Figure 3)
#MSE results of the three X-M model Scenarios with different Strength of Modification (SoM) using RFQT and single-stratification

ATE_MSE_<-c()#ATE
DR_S_MSE_<-c();R_S_MSE_<-c();N_S_MSE_<-c()#single stratification
DR_RFQT_MSE_<-c();R_RFQT_MSE_<-c();N_RFQT_MSE_<-c() #RFQT

for(s in c('A','B','C')){
  scenario_used<-s
  Gammas<-c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    )#共6中SoM情况
  ATE_MSE<-c()#ATE
  DR_S_MSE<-c();R_S_MSE<-c();N_S_MSE<-c()#single stratification
  DR_RFQT_MSE<-c();R_RFQT_MSE<-c();N_RFQT_MSE<-c() #RFQT
  for(i in 1:length(Gammas)){
    set.seed(10*i)
    gamma_used<-Gammas[i]
    Dat<-getDat(scenario=scenario_used, SoM=gamma_used)
    odat<-Dat$traning.set  #training set in 全局环境
    vdat<-Dat$testing.set  #testing set in 全局环境
    
    ##############ATE
    #odat#用training set估计
    fitGX<-lm(    odat[,3]~  odat[,2]  )  ; bx<-as.numeric( summary(fitGX)$coef[-1,1]  ); bxse<-as.numeric(  summary(fitGX)$coef[-1,2])
    fitGY<-lm(    odat[,4]~  odat[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
    MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
    ATE_MSE<-c(   ATE_MSE     ,        mean(  (vdat$true_STE - MRres@Estimate )^2 ) )
    
    
    ###############Single Stratification (still random forest)
    ###DR
    ALLRES<-RFQTfit(odat,vdat,Nb=200,SingleM=TRUE,method='DR')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\S_DR_',scenario_used,'_',gamma_used,'.RData'))
    DR_S_MSE<-c(  DR_S_MSE , ALLRES$MSE_test[length(ALLRES$MSE_test)] )
    
    ###R
    ALLRES<-RFQTfit(odat,vdat,Nb=200,SingleM=TRUE,method='Residual')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\S_R_',scenario_used,'_',gamma_used,'.RData'))
    R_S_MSE<-c(  R_S_MSE , ALLRES$MSE_test[length(ALLRES$MSE_test)] )
    
    ###N
    ALLRES<-RFQTfit(odat,vdat,Nb=200,SingleM=TRUE,method='N')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\S_N_',scenario_used,'_',gamma_used,'.RData'))
    N_S_MSE<-c(  N_S_MSE , ALLRES$MSE_test[length(ALLRES$MSE_test)] )
    
    
    ################RFQT
    ###DR
    ALLRES<-RFQTfit(odat,vdat,Nb=200,method='DR')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_DR_',scenario_used,'_',gamma_used,'.RData'))
    DR_RFQT_MSE<-c(  DR_RFQT_MSE , ALLRES$MSE_test[length(ALLRES$MSE_test)] )
    
    ###R
    ALLRES<-RFQTfit(odat,vdat,Nb=200,method='Residual')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_R_',scenario_used,'_',gamma_used,'.RData'))
    R_RFQT_MSE<-c(  R_RFQT_MSE , ALLRES$MSE_test[length(ALLRES$MSE_test)] )
    
    ###N
    ALLRES<-RFQTfit(odat,vdat,Nb=200,method='N')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_N_',scenario_used,'_',gamma_used,'.RData'))
    N_RFQT_MSE<-c(  N_RFQT_MSE , ALLRES$MSE_test[length(ALLRES$MSE_test)] )
    
  } #end of all strength of modification for a specific scenario
  #ATE
  ATE_MSE_<-rbind(ATE_MSE_,ATE_MSE)
  #single stratification
  DR_S_MSE_<-rbind(DR_S_MSE_,DR_S_MSE)
  R_S_MSE_<-rbind(R_S_MSE_,R_S_MSE)
  N_S_MSE_<-rbind(N_S_MSE_,N_S_MSE)
  #RFQT
  DR_RFQT_MSE_<-rbind(DR_RFQT_MSE_,DR_RFQT_MSE)
  R_RFQT_MSE_<-rbind(R_RFQT_MSE_,R_RFQT_MSE)
  N_RFQT_MSE_<-rbind(N_RFQT_MSE_,N_RFQT_MSE)
}






###result2: (Figure 7)
#MSE and individual predicted effect curve changing along the number ofQ trees

set.seed(1123)
Dat<-getDat(scenario='A', SoM=0.5) #simulated data  #the default setting: scenario='A' and SoM=0.5
odat<-Dat$traning.set  #training set in 全局环境
vdat<-Dat$testing.set

ALLRES<-RFQTfit(odat,vdat,Nb=500,method='DR')#Nb=500

saveRDS(ALLRES,file='D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\result2and3.RData')
#ALLRES<-readRDS('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\result2and3.RData')

RES<-ALLRES$RES
MSE1<-getMSE( RES , 1 ) ; MSE2<-getMSE( RES , 2 )  #MSE1: OOB error; MSE2: test error
#getMSE()中odat vdat使用全局环境
plot( 1:length(MSE1) , MSE1, type='l',xlab='Number of Q trees', ylab='MSE' )
lines(     1:length(MSE2) , MSE2, col='blue'   )
legend('topright', legend=c('OOB error', 'Test error'),
       col=c("black", "blue"), lty=1, cex=1.0)
#Figure7_left  500*450

predict_matrix2<-getPredict( RES,2     )#2: test set predicts
plot(   1:ncol( predict_matrix2   )   ,  predict_matrix2[1,],  type='l', ylim=c(0.1,1.2),
        xlab='Number of Q trees', ylab='Predicted effect')
lines(     1:ncol( predict_matrix2   )   ,  predict_matrix2[2,], col='blue'   )
lines(     1:ncol( predict_matrix2   )   ,  predict_matrix2[3,], col='red'   )
lines(     1:ncol( predict_matrix2   )   ,  predict_matrix2[4,], col='green'   )

legend('topright', legend=c('Individual One', 'Individual Two', 'Individual Three', 'Individual Four'),
       col=c("black", "blue",'red','green'), lty=1, cex=0.7)
#Figure7_right 500*450



##Result3: (new Figure)
#according to vdat predicted HTE (i.e. v_predict), refit HTE ~ M to visualize the results or for future guidance
#ALLRES就用results2 里面的即可

ALLRES<-RFQTfit(odat,vdat,Nb=50,method='DR',S=7,endsize=1000)
RES<-ALLRES$RES
predict_matrix2<-getPredict( RES,2     )
dim(predict_matrix2)
dim(vdat)

#simple linear regression
fit<-lm(  predict_matrix2[,ncol(predict_matrix2)]~ as.matrix(vdat[,5:24])     )
summary(fit)













