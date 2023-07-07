


#Real application code

Dat<-read.csv(  paste0('D:\\files\\R new\\Precison_Medicine\\real_data\\Dat.csv'),header=T,na.strings ="?"  )

#Dat: the real data with the same form as the simulated data (i.e ID+Z+X+Y+M) no HTE labels
dim(Dat) #142186     32

apply(Dat, 2, is.numeric)
apply(Dat, 2, is.character)

NNNN<-dim(Dat)[1]
NNN<-round(NNNN*2/3)
NNNN;NNN#For male: 142186(NNNN) 94791(NNN)

Dat<-as.data.frame(Dat)
names(Dat)[1]<-'I'
names(Dat)[3]<-'X'
names(Dat)[4]<-'Y'
colnames(Dat)

#同一个Dat下，odat和vdat永远不会改变 (odat vdat在全局环境)
vdat<-Dat[(NNN+1):NNNN, ]  #validation dat
odat<-Dat[1:NNN,]   #original (training) dat


#RFQT fitting
ALLRES_real<-RFQTfit(odat,vdat,Nb=200,method='DR')

#saveRDS(ALLRES_real,file='D:\\files\\R new\\Precison_Medicine\\real_data\\ALLRES_real.RData')
#ALLRES_real<-readRDS('D:\\files\\R new\\Precison_Medicine\\real_data\\ALLRES_real.RData')


###result1: (Figure 4)-------------------------------------------------------------------------------------------------------------------------
#histogram plot of the testing set predicts using (DR Nb=200)----------------------------------------------------------------------------------
summary(  ALLRES_real$Predicts_test[,200]  )
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.   (honest)
#-0.717172 -0.028483 -0.019300 -0.020037 -0.009268  0.531695 
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.   (nonhonest)
#-0.065441 -0.027532 -0.017536 -0.014073 -0.003499  0.061977

hist( ALLRES_real$Predicts_test[,100]  , n=200, main='', xlab='Predicted effect')
#hist( ALLRES_real$Predicts_test[,200]  , n=500, xlim=c(-0.1,0.05) , main='', xlab='Predicted effect')

ggdata<-data.frame(Predicted_effect= ALLRES_real$Predicts_test[,100] )
p<-ggplot(ggdata) +
  geom_bar(aes(Predicted_effect)) +
  scale_x_binned(n.breaks=30)+
  geom_abline(intercept = 0, slope = 1000000000000000000, colour='grey', linewidth=1.0 ,linetype = 2)+
  xlab('Predicted effect')+ylab('Count')+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p
ggsave(paste0('Fig4.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots'), 
       height = 5, width = 15, units = "in",limitsize=TRUE)


###results2: (Figure S2)------------------------------------------------------------------------------------------------------------------
#Variable importance plot-----------------------------------------------------------------------------------------------------------------
VIindex<-order(ALLRES_real$VI2)#得出一组position index用以表示最小到最大unit的各自的位置 

ggdata_VI<-data.frame(candidate=factor(names(Dat)[VIindex+4],levels =(names(Dat)[VIindex+4])[length(VIindex):1] ),  #逆向排序
                      VI=( ALLRES_real$VI2 )[VIindex]   )  
p<-ggplot(ggdata_VI,aes(x=candidate,y=VI))+geom_col( )+
  xlab('Candidate covariate')+ylab('Variable Importance')+#labs(title=c('Doubly-ranked method', 'Residual method' , 'Naive method'   )[RESindex] )+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  coord_flip()
p
ggsave(paste0('Fig_S2.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots'), 
       height = 6, width = 4, units = "in",limitsize=TRUE)


###result3: (Figure 5)------------------------------------------------------------------------------------------------------------------------
#important variable -> marginal stratum analysis using pre-determined (tenth) quantile points-------------------------------------------------
#one can use DRMR package
#Q statistic use average; point and CIS use Rubin's Rule (RR)?

VInames<-(names(Dat)[VIindex+4])[length(VIindex):1]
VInames[1:6]  #"dbp"    "hip"    "monos"  "wt"     "ages"   "vitcap"
VInames_chosen<- c(VInames[1:6], 'bmi', 'ht') #不再使用bmi.1这个string name

dim(Dat)
names(Dat)
names(Dat)[10]
names(Dat)[10]<-'bmi'

Label_names<-c('Diastolic blood pressure (mmHg)',
               'Hip circumference (cm)',
               'Monocyte count' ,
               'Weight (kg)',
               'Age at survey (years)',
               'Vital capacity (litres)',
               'Body mass index' ,
               'Height (cm)'
               )

#devtools::install_github("HDTian/DRMR")
library(DRMR)

for(i in c(1,2,4,5,6,8)){
  Mname<-VInames_chosen[i]  #M: marginal specific covariate name
  print(Mname)
  dat<-Dat[,     c(   2 ,3,4, match(Mname,names(Dat)   ) ) ]
  names(dat)<-c(  'Z', 'X', 'Y', 'M' ) #for DRMR
  apply(dat, 2, is.numeric) #TRUE TRUE TRUE TRUE
  
  
  #Due to the rank variability: first randomly remove 10 individuals
  EST_ZX<-c();EST_ZY<-c() #dim(EST) 10*1000
  SE_ZX<-c();SE_ZY<-c() #dim(SE) 10*1000
  MEANX<-c() #dim(MEANX) 10*1000
  mtimes<-500
  for(m in 1:mtimes){
    dat_used<-dat[-sample(1:nrow(dat), 10)  ,  ]
    
    #stratification via DRMR and get the stratum-specific (est and s.e.) and Xmean
    rdat<-Stratify(dat_used,onExposure = FALSE)#augmented dataset
    SRES<-getSummaryInf(rdat) #Stratification RESults
    
    EST_ZX<-cbind( EST_ZX, SRES$DRres$bx );EST_ZY<-cbind( EST_ZY, SRES$DRres$by )
    SE_ZX<-cbind(SE_ZX, SRES$DRres$bxse  );SE_ZY<-cbind(SE_ZY, SRES$DRres$byse  )
    MEANX<-cbind( MEANX , SRES$DRres$mean)
  }
  ###Rubin's Rule
  #for ZX association
  Qbar_ZX<-apply( EST_ZX , 1, mean )
  TT1_ZX<-apply( SE_ZX^2 , 1, mean  ) #i.e. Ubar  #a vector
  TT2_ZX<-(1+ 1/mtimes)*(  apply(  ( EST_ZX  -   Qbar_ZX  )^2 ,1 , sum    )/(mtimes-1)  ) #a vector 
  TT_ZX<-TT1_ZX+TT2_ZX#a vector
  Fdf_ZX<-(mtimes -1 )*(1+ TT1_ZX/TT2_ZX)^2   #the second df for F statistic (越大，使得sqrt(F_1,df))越接近Z statistic
  #qf(0.95, 1,Fdf_ZX )#still a vector
  
  #for ZY association
  Qbar_ZY<-apply( EST_ZY , 1, mean )
  TT1_ZY<-apply( SE_ZY^2 , 1, mean  ) #i.e. Ubar  #a vector
  TT2_ZY<-(1+ 1/mtimes)*(  apply(  ( EST_ZY  -   Qbar_ZY  )^2 ,1 , sum    )/(mtimes-1)  ) #a vector 
  TT_ZY<-TT1_ZY+TT2_ZY#a vector
  Fdf_ZY<-(mtimes -1 )*(1+ TT1_ZY/TT2_ZY)^2   #the second df for F statistic (越大，使得sqrt(F_1,df))越接近Z statistic
  #qf(0.95, 1,Fdf_ZY )#still a vector
  
  #for MR est: simply ignore the uncertainty of the estimated ZX association
  MRest<-Qbar_ZY/Qbar_ZX
  
  ##MR estimand CI_low and CI_up based on F test statistic with ignoreing the estimated ZX association
  ggdata<-data.frame( X = apply( MEANX ,1 ,mean)  , 
                      Est = MRest    , 
                      CI_low = (  Qbar_ZY- sqrt( qf(0.95, 1,Fdf_ZY ) )*sqrt( TT_ZY )   )/Qbar_ZX, 
                      CI_up =  (   Qbar_ZY+ sqrt( qf(0.95, 1,Fdf_ZY ) )*sqrt( TT_ZY )  )/Qbar_ZX 
                      )
  ##基于Multiple imputation后的stratum specific结果算Q statistic？
  #here each stratum-specific ZX ZY associations are assumed to be Gaussian distributed
  #that is, Qbar_ZX ~ N(  Q_ZX , sqrt(TT_ZX) ^2 ) where Q_ZX represents the LACE estimand
  
  Bx<-Qbar_ZX;Bxse<-sqrt(TT_ZX)
  By<-Qbar_ZY;Byse<-sqrt(TT_ZY)
  #先算effect via naive IVW under null (no effect difference) - i.e. ignore the ZX estimate uncertainty
  Sr<-mr_ivw(mr_input(Bx, (1:10) , By, Byse))
  delta<-Sr@Estimate #一维的
  #Cochran's Q statistic
  ssigma_square <- Byse^2 + delta^2*Bxse^2
  QQ<-sum( ( By - delta*Bx     )^2 /ssigma_square      )
  pvalue<-1-pchisq(QQ,9)
  
  
  #visualization
  
  if(!i%in%c(3,7)){
    p <- ggplot(ggdata, aes(X, Est))+
      geom_point(ggdata, mapping =aes(X, Est), alpha=1,size=2  )+
      geom_errorbar(data=ggdata,mapping=aes(x=X,ymin=CI_low,ymax=CI_up),size = 0.25) +
      geom_hline(aes(yintercept = 0),linetype='dashed',alpha=1,linewidth=1,color='grey')+
      labs(x=paste0(Label_names[i], " (Q statistic: ", sprintf("%.1f", round(QQ,1)) , "; pvalue: ",  sprintf("%.1f",round(pvalue,3) ) ,")"  ),
           y='Stratum-specific estimates')+
      coord_cartesian(ylim = c(-0.07,0.07) )+ 
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  }
  
  if(i==3){
    p <- ggplot(ggdata, aes(X, Est))+
      geom_point(ggdata, mapping =aes(X, Est), alpha=1,size=2  )+
      geom_errorbar(data=ggdata,mapping=aes(x=X,ymin=CI_low,ymax=CI_up),size = 0.5) +
      geom_hline(aes(yintercept = 0),linetype='dashed',alpha=1,linewidth=1,color='grey')+
      labs(x=expression( paste('Monocyte count (',10^9,'cells/Litre)', " (Q statistic: 4.0; pvalue: 0.910)"  ) ),
           y='Stratum-specific estimates')+
      coord_cartesian(ylim = c(-0.07,0.07) )+   #expression会使得所有目标都直意表达
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  }
  
  
  if(i==7){
    p <- ggplot(ggdata, aes(X, Est))+
      geom_point(ggdata, mapping =aes(X, Est), alpha=1,size=2  )+
      geom_errorbar(data=ggdata,mapping=aes(x=X,ymin=CI_low,ymax=CI_up),size = 0.5) +
      geom_hline(aes(yintercept = 0),linetype='dashed',alpha=1,linewidth=1,color='grey')+
      labs(x=expression( paste('Body mass index (' , kg/m^2  , ')', " (Q statistic: 9.04; pvalue: 0.434)"  ) ),y='Stratum-specific estimates')+
      coord_cartesian(ylim = c(-0.07,0.07) )+   #expression会使得所有目标都直意表达
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  }
  print(p)
  ggsave(paste0('Fig5_',  i  , '.eps' ), plot = p , 
          path='C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots', height = 4, width = 5, units = "in",limitsize=TRUE)
}







###result4: (new Figure)---------------------------------------------------------------------------------------------------------
#according to vdat predicted HTE (i.e. v_predict), refit HTE ~ M to visualize the results or for future guidance-----------------




###result5: (Figure S1)----------------------------------------------------------------------------------------------------------
#permutation test results--------------------------------------------------------------------------------------------------------




















