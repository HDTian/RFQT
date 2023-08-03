#sim_supplementary


#agreed to: increase the iteration and use sampling way to reduce varibaility




scenario_used<-'C'
Gammas<-c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    )#共6中
i=3 #gamma_used = 0.2

MSE<-c()
for(iii in 1:10){
  set.seed(10*i*100+iii)  #set.seed(10*i)
  gamma_used<-Gammas[i]
  
  #get data
  Dat<-getDat(scenario=scenario_used, SoM=gamma_used)
  odat<-Dat$traning.set  #training set in 全局环境
  vdat<-Dat$testing.set  #testing set in 全局环境
  
  
  ###############Decision Tree
  DTRES<-DTfit(odat,vdat,method='Residual')
  MSE<-c( MSE, DTRES$MSE )
}


###DT paralle working
DTRES_par<-function(seed){
  set.seed(seed)  #set.seed(10*i)
  
  #get data
  Dat<-getDat(scenario=scenario_used, SoM=gamma_used)
  odat<-Dat$traning.set  #training set in 全局环境
  vdat<-Dat$testing.set  #testing set in 全局环境
  
  
  ###############Decision Tree
  DTRES<-DTfit(odat,vdat,method='Residual')
  
  #
  return( DTRES$MSE )
}
cl<-makeCluster(detectCores()-1)
clusterEvalQ(cl=cl , expr=library(dplyr))
clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
clusterExport(  cl=cl ,  varlist=c( 'scenario_used', 'gamma_used' ,
                                    'getDat','GetIndex', 'GetTree' , 'GetNindex' , 'DTfit'))
DTMSERES<-parSapply(   cl , 10*i*100+(1:100) , DTRES_par )
stopCluster(cl)

  












#for scenarioB&C, and for Residual and Naive, we see the traceplot is not stable as DR 
#try Halve=TRUE to see whether has stable traceplot

#HT: Halve=True

for(s in c('B','C')){
  scenario_used<-s
  Gammas<-c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    )#共6中SoM情况
  for(i in 1:length(Gammas)){
    set.seed(10*i)
    gamma_used<-Gammas[i]
    
    #get data
    Dat<-getDat(scenario=scenario_used, SoM=gamma_used)
    odat<-Dat$traning.set  #training set in 全局环境
    vdat<-Dat$testing.set  #testing set in 全局环境
    
    
    ###############Decision Tree
    ###DR
    #DTRES<-DTfit(odat,vdat,method='DR',Halve=TRUE)
    #saveRDS(DTRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_DR_',scenario_used,'_',gamma_used,'HT.RData'))
    #DTRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_DR_',scenario_used,'_',gamma_used,'HT.RData'))
    
    ###R
    DTRES<-DTfit(odat,vdat,method='Residual',Halve=TRUE)
    saveRDS(DTRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_R_',scenario_used,'_',gamma_used,'HT.RData'))
    
    ###N
    DTRES<-DTfit(odat,vdat,method='N',Halve=TRUE)
    saveRDS(DTRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_N_',scenario_used,'_',gamma_used,'HT.RData'))
    
    
    ################RFQT
    ###DR
    #ALLRES<-RFQTfit(odat,vdat,Nb=200,method='DR',Halve=TRUE)
    #saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_DR_',scenario_used,'_',gamma_used,'HT.RData'))
    #ALLRES<-readRDS(paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_DR_',scenario_used,'_',gamma_used,'HT.RData'))
    #ALLRES$MSE_test[length(ALLRES$MSE_test)] 
    
    ###R
    ALLRES<-RFQTfit(odat,vdat,Nb=200,method='Residual',Halve=TRUE)
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_R_',scenario_used,'_',gamma_used,'HT.RData'))
    
    ###N
    ALLRES<-RFQTfit(odat,vdat,Nb=200,method='N',Halve=TRUE)
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_N_',scenario_used,'_',gamma_used,'HT.RData'))
    
  } #end of all strength of modification for a specific scenario
  
}



ggdata<-c()
ggATE<-c()
scenario_used<-'C'
Gammas<-c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    )#共6中SoM情况
for(i in 1:length(Gammas)){
  set.seed(10*i)
  gamma_used<-Gammas[i]
  
  #get data
  Dat<-getDat(scenario=scenario_used, SoM=gamma_used)
  odat<-Dat$traning.set  #training set in 全局环境
  vdat<-Dat$testing.set  #testing set in 全局环境
  
  ##############ATE
  fitGX<-lm(    odat[,3]~  odat[,2]  )  ; bx<-as.numeric( summary(fitGX)$coef[-1,1]  ); bxse<-as.numeric(  summary(fitGX)$coef[-1,2])
  fitGY<-lm(    odat[,4]~  odat[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
  MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
  
  ggATE<-rbind(ggATE ,  
               c( mean(  (vdat$true_STE - MRres@Estimate )^2 )  , gamma_used , scenario_used    )
  )
  
  ##############Decision Tree
  DTRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_DR_',scenario_used,'_',gamma_used,'.RData'))
  ggdata<-rbind(ggdata,
                c(  mean( ( DTRES$v_predict - vdat$true_STE   )^2   ),gamma_used ,scenario_used, 'DR', 'DT'  )
  )
  DTRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_R_',scenario_used,'_',gamma_used,'HT.RData'))
  ggdata<-rbind(ggdata,
                c(  mean( ( DTRES$v_predict - vdat$true_STE   )^2   ),gamma_used ,scenario_used, 'R', 'DT'  )
  )
  DTRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_N_',scenario_used,'_',gamma_used,'HT.RData'))
  ggdata<-rbind(ggdata,
                c(  mean( ( DTRES$v_predict - vdat$true_STE   )^2   ),gamma_used ,scenario_used, 'N', 'DT'  )
  )
  
  
  ##############Random Forest
  ALLRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_DR_',scenario_used,'_',gamma_used,'.RData'))
  ggdata<-rbind(ggdata,
                c(  ALLRES$MSE_test[length(ALLRES$MSE_test)] ,gamma_used ,scenario_used, 'DR', 'RFQT'  )
  )
  ALLRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_R_',scenario_used,'_',gamma_used,'HT.RData'))
  ggdata<-rbind(ggdata,
                c(  ALLRES$MSE_test[length(ALLRES$MSE_test)] ,gamma_used ,scenario_used, 'R', 'RFQT'  )
  )
  ALLRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_N_',scenario_used,'_',gamma_used,'HT.RData'))
  ggdata<-rbind(ggdata,
                c(  ALLRES$MSE_test[length(ALLRES$MSE_test)] ,gamma_used ,scenario_used, 'N', 'RFQT'  )
  )
  
}




ggATE<-as.data.frame( ggATE )
ggdata<-as.data.frame( ggdata  )
names(ggATE)<- c('MSE', 'Gammas','Scenario')
names(ggdata)<-c('MSE', 'Gammas','Scenario', 'de_collider','method_type')
ggATE$MSE<-as.numeric( ggATE$MSE   ); ggATE$Gammas<-as.numeric( ggATE$Gammas   )
ggdata$MSE<-as.numeric( ggdata$MSE   ); ggdata$Gammas<-as.numeric( ggdata$Gammas   )
ggdata$method_type<-factor(   ggdata$method_type , levels =c( 'RFQT','DT'   )      )
ggdata$de_collider<-factor(   ggdata$de_collider , levels =c( 'DR','R', 'N'   )      )


library(facetscales)
pd <- position_dodge(0.001)
scales_y <- list(
  `A` = scale_y_continuous(limits = c(0, 1.5)),
  `B` = scale_y_continuous(limits = c(0, 8)),
  `C` = scale_y_continuous(limits = c(0, 10))
)

p<- ggplot(data = ggdata, mapping = aes(x = Gammas, y = MSE))+
  geom_line(mapping = aes(color= de_collider , linetype = method_type),size=1.0,alpha=1.0,position = pd)+
  geom_point(mapping = aes(color= de_collider ,  shape=de_collider ,  group = de_collider ) ,size = 2.5,alpha=1.0,position = pd)+
  geom_line(  data=ggATE ,  mapping = aes(x = Gammas, y = MSE), linetype=3 ,linewidth = 1.0 )+
  geom_point(  data=ggATE ,  mapping = aes(x = Gammas, y = MSE),size=2 )+
  labs(x = "Strength of modification", y = "MSE")+
  facet_grid_sc(rows=vars(Scenario), scales = list(y = scales_y))+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  guides(  color=guide_legend(title='Stratification Method'),
           shape=guide_legend(title='Stratification Method'),
           group=guide_legend(title='Stratification Method')  , linetype=guide_legend(title='Estimation Method')  )  
p



#######################------------------------------------------------------------------------------------
#######################------------------------------------------------------------------------------------
#######################------------------------------------------------------------------------------------
#######################------------------------------------------------------------------------------------
#######################------------------------------------------------------------------------------------
#######################------------------------------------------------------------------------------------
#######################------------------------------------------------------------------------------------
#endsize = 5000

for(s in c('B','C')){
  scenario_used<-s
  Gammas<-c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    )#共6中SoM情况
  for(i in 1:length(Gammas)){
    set.seed(10*i)
    gamma_used<-Gammas[i]
    
    #get data
    Dat<-getDat(scenario=scenario_used, SoM=gamma_used)
    odat<-Dat$traning.set  #training set in 全局环境
    vdat<-Dat$testing.set  #testing set in 全局环境
    
    
    ###############Decision Tree
    ###DR
    #DTRES<-DTfit(odat,vdat,method='DR',endsize=5000)
    #saveRDS(DTRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_DR_',scenario_used,'_',gamma_used,'_5000.RData'))
    #DTRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_DR_',scenario_used,'_',gamma_used,'_5000.RData'))
    
    ###R
    #DTRES<-DTfit(odat,vdat,method='Residual',endsize=5000)
    #saveRDS(DTRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_R_',scenario_used,'_',gamma_used,'_5000.RData'))
    
    ###N
    DTRES<-DTfit(odat,vdat,method='N',endsize=5000)
    saveRDS(DTRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_N_',scenario_used,'_',gamma_used,'_5000.RData'))
    
    
    ################RFQT
    ###DR
    #ALLRES<-RFQTfit(odat,vdat,Nb=200,method='DR',endsize=5000)
    #saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_DR_',scenario_used,'_',gamma_used,'_5000.RData'))
    #ALLRES<-readRDS(paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_DR_',scenario_used,'_',gamma_used,'_5000.RData'))
    #ALLRES$MSE_test[length(ALLRES$MSE_test)] 
    
    ###R
    #ALLRES<-RFQTfit(odat,vdat,Nb=200,method='Residual',endsize=5000)
    #saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_R_',scenario_used,'_',gamma_used,'_5000.RData'))
    
    ###N
    ALLRES<-RFQTfit(odat,vdat,Nb=200,method='N',endsize=5000)
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_N_',scenario_used,'_',gamma_used,'_5000.RData'))
    
  } #end of all strength of modification for a specific scenario
  
}



ggdata<-c()
ggATE<-c()
scenario_used<-'C'
Gammas<-c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    )#共6中SoM情况
for(i in 1:length(Gammas)){
  set.seed(10*i)
  gamma_used<-Gammas[i]
  
  #get data
  Dat<-getDat(scenario=scenario_used, SoM=gamma_used)
  odat<-Dat$traning.set  #training set in 全局环境
  vdat<-Dat$testing.set  #testing set in 全局环境
  
  ##############ATE
  fitGX<-lm(    odat[,3]~  odat[,2]  )  ; bx<-as.numeric( summary(fitGX)$coef[-1,1]  ); bxse<-as.numeric(  summary(fitGX)$coef[-1,2])
  fitGY<-lm(    odat[,4]~  odat[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
  MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
  
  ggATE<-rbind(ggATE ,  
               c( mean(  (vdat$true_STE - MRres@Estimate )^2 )  , gamma_used , scenario_used    )
  )
  
  ##############Decision Tree
  DTRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_DR_',scenario_used,'_',gamma_used,'_5000.RData'))
  ggdata<-rbind(ggdata,
                c(  mean( ( DTRES$v_predict - vdat$true_STE   )^2   ),gamma_used ,scenario_used, 'DR', 'DT'  )
  )
  DTRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_R_',scenario_used,'_',gamma_used,'_5000.RData'))
  ggdata<-rbind(ggdata,
                c(  mean( ( DTRES$v_predict - vdat$true_STE   )^2   ),gamma_used ,scenario_used, 'R', 'DT'  )
  )
  DTRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_N_',scenario_used,'_',gamma_used,'_5000.RData'))
  ggdata<-rbind(ggdata,
                c(  mean( ( DTRES$v_predict - vdat$true_STE   )^2   ),gamma_used ,scenario_used, 'N', 'DT'  )
  )
  
  
  ##############Random Forest
  ALLRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_DR_',scenario_used,'_',gamma_used,'_5000.RData'))
  ggdata<-rbind(ggdata,
                c(  ALLRES$MSE_test[length(ALLRES$MSE_test)] ,gamma_used ,scenario_used, 'DR', 'RFQT'  )
  )
  ALLRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_R_',scenario_used,'_',gamma_used,'_5000.RData'))
  ggdata<-rbind(ggdata,
                c(  ALLRES$MSE_test[length(ALLRES$MSE_test)] ,gamma_used ,scenario_used, 'R', 'RFQT'  )
  )
  ALLRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_N_',scenario_used,'_',gamma_used,'_5000.RData'))
  ggdata<-rbind(ggdata,
                c(  ALLRES$MSE_test[length(ALLRES$MSE_test)] ,gamma_used ,scenario_used, 'N', 'RFQT'  )
  )
  
}




ggATE<-as.data.frame( ggATE )
ggdata<-as.data.frame( ggdata  )
names(ggATE)<- c('MSE', 'Gammas','Scenario')
names(ggdata)<-c('MSE', 'Gammas','Scenario', 'de_collider','method_type')
ggATE$MSE<-as.numeric( ggATE$MSE   ); ggATE$Gammas<-as.numeric( ggATE$Gammas   )
ggdata$MSE<-as.numeric( ggdata$MSE   ); ggdata$Gammas<-as.numeric( ggdata$Gammas   )
ggdata$method_type<-factor(   ggdata$method_type , levels =c( 'RFQT','DT'   )      )
ggdata$de_collider<-factor(   ggdata$de_collider , levels =c( 'DR','R', 'N'   )      )


library(facetscales)
pd <- position_dodge(0.001)
scales_y <- list(
  `A` = scale_y_continuous(limits = c(0, 1.5)),
  `B` = scale_y_continuous(limits = c(0, 8)),
  `C` = scale_y_continuous(limits = c(0, 10))
)

p<- ggplot(data = ggdata, mapping = aes(x = Gammas, y = MSE))+
  geom_line(mapping = aes(color= de_collider , linetype = method_type),size=1.0,alpha=1.0,position = pd)+
  geom_point(mapping = aes(color= de_collider ,  shape=de_collider ,  group = de_collider ) ,size = 2.5,alpha=1.0,position = pd)+
  geom_line(  data=ggATE ,  mapping = aes(x = Gammas, y = MSE), linetype=3 ,linewidth = 1.0 )+
  geom_point(  data=ggATE ,  mapping = aes(x = Gammas, y = MSE),size=2 )+
  labs(x = "Strength of modification", y = "MSE")+
  facet_grid_sc(rows=vars(Scenario), scales = list(y = scales_y))+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  guides(  color=guide_legend(title='Stratification Method'),
           shape=guide_legend(title='Stratification Method'),
           group=guide_legend(title='Stratification Method')  , linetype=guide_legend(title='Estimation Method')  )  
p



###----------------------------------------------------------------------------------
###----------------------------------------------------------------------------------
###focus on scenarioC and i=6 (SoM=0.5) try multiple iterations and see the RFQTfit results
#using Residual stratification and Nb=200 endsize=1000 
scenario_used<-'C'
Gammas<-c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    )#共6中SoM情况
i<-6
for(iii in 1:100){
  set.seed(10*i*100+iii)  #set.seed(10*i)
  gamma_used<-Gammas[i]
  #get data
  Dat<-getDat(scenario=scenario_used, SoM=gamma_used)
  odat<-Dat$traning.set  #training set in 全局环境
  vdat<-Dat$testing.set  #testing set in 全局环境
  
  
  
  #RFQT  -  Residual
  ALLRES<-RFQTfit(odat,vdat,Nb=200,method='Residual',Cores=8) #full cores
  saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_R_',scenario_used,'_',gamma_used,'_iter',iii,'.RData'))
}

MSEs<-c()
for(iii in 1:66){
  ALLRES<-readRDS(paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_R_',scenario_used,'_',gamma_used,'_iter',iii,'.RData'))
  MSEs<-c(  MSEs ,ALLRES$MSE_test[200]  )
}
var(MSEs)
var(MSEs)/66
hist(MSEs,n=66)



