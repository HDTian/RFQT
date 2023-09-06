


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

hist( ALLRES_real$Predicts_test[,200]  , n=200, main='', xlab='Predicted effect')
#hist( ALLRES_real$Predicts_test[,200]  , n=500, xlim=c(-0.1,0.05) , main='', xlab='Predicted effect')

ggdata<-data.frame(Predicted_effect= ALLRES_real$Predicts_test[,200] )
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

for(i in c(1:8)){
  Mname<-VInames_chosen[i]  #M: marginal specific covariate name
  print(Mname)
  dat<-Dat[,     c(   2 ,3,4, match(Mname,names(Dat)   ) ) ]
  names(dat)<-c(  'Z', 'X', 'Y', 'M' ) #for DRMR
  apply(dat, 2, is.numeric) #TRUE TRUE TRUE TRUE
  
  
  #Due to the rank variability: first randomly remove 10 individuals
  set.seed(i)
  
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
  
  ##also store some other metrics (possible for some trend test)
  ggdata$ZXEst<-Qbar_ZX
  ggdata$ZYEst<-Qbar_ZY
  ggdata$TT_ZX<-TT_ZX
  ggdata$TT_ZY<-TT_ZY
  ggdata$Fdf_ZX<-Fdf_ZX
  ggdata$Fdf_ZY<-Fdf_ZY
  
  
  ##基于Multiple imputation后的stratum specific结果算Q statistic
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
  print(c(QQ,pvalue)  )
  
  
  write.csv(ggdata, paste0('D:\\files\\R new\\Precison_Medicine\\real_data\\',"ggdata_",i, '_',QQ,'_',pvalue  , '.csv'), row.names=F)
  #visualization
  
  if(!i%in%c(3,7)){
    p <- ggplot(ggdata, aes(X, Est))+
      geom_point(ggdata, mapping =aes(X, Est), alpha=1,size=2  )+
      geom_errorbar(data=ggdata,mapping=aes(x=X,ymin=CI_low,ymax=CI_up),size = 0.5) +
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
  #ggsave(paste0('newFig5_',  i  , '.eps' ), plot = p , 
  #        path='C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots', height = 4, width = 5, units = "in",limitsize=TRUE)
}

#or direct draw the plots:
Label_names<-c('Diastolic blood pressure (mmHg)',
               'Hip circumference (cm)',
               'Monocyte count' ,
               'Weight (kg)',
               'Age at survey (years)',
               'Vital capacity (litres)',
               'Body mass index' ,
               'Height (cm)'
)

refmatrix<-rbind(
  c(9.63049607371929,0.381217953296086),
  c(23.1870338714598,0.00578965359286898),
  c(4.06792187403348,0.906885415485453),
  c(12.8039869844824,0.171678073383885),
  c(7.430776921836,0.592362972147522),
  c(12.6542435402473,0.178885527324497),
  c(8.86472403604516,0.449854028150605),
  c(4.44112943737513,0.880058735068154)
)


for(i in c(1:8)){
  QQ<-refmatrix[i,1]; pvalue<-refmatrix[i,2]
  #write.csv(ggdata, paste0('D:\\files\\R new\\Precison_Medicine\\real_data\\',"ggdata_",i, '_',QQ,'_',pvalue  , '.csv'), row.names=F)
  #visualization
  
  ggdata<-read.csv( paste0('D:\\files\\R new\\Precison_Medicine\\real_data\\',"ggdata_",i, '_',QQ,'_',pvalue  , '.csv'), header=T,na.strings ="?" )
  
  if(!i%in%c(3,7)){
    p <- ggplot(ggdata, aes(X, Est))+
      geom_point(ggdata, mapping =aes(X, Est), alpha=1,size=2  )+
      geom_errorbar(data=ggdata,mapping=aes(x=X,ymin=CI_low,ymax=CI_up),size = 0.5) +
      geom_hline(aes(yintercept = 0),linetype='dashed',alpha=1,linewidth=1,color='grey')+
      labs(x=paste0(Label_names[i], " (Q statistic: ", sprintf("%.1f", round(QQ,1)) , "; pvalue: ",  sprintf("%.3f",round(pvalue,3) ) ,")"  ),
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
      labs(x=expression( paste('Body mass index (' , kg/m^2  , ')', " (Q statistic: 8.9; pvalue: 0.450)"  ) ),y='Stratum-specific estimates')+
      coord_cartesian(ylim = c(-0.07,0.07) )+   #expression会使得所有目标都直意表达
      theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  }
  print(p)
  ggsave(paste0('newFig5_',  i  , '.eps' ), plot = p , 
         path='C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots', height = 4, width = 5, units = "in",limitsize=TRUE)
}




#calculate according to stored ggdata
library(metafor)
for(i in c(1:8)){
  ggdata<-read.csv( paste0('D:\\files\\R new\\Precison_Medicine\\real_data\\',"ggdata_",i  , '.csv'), header=T,na.strings ="?" )
  #meta test
  ggdata$se<-sqrt(  ggdata$TT_ZY )/ggdata$ZXEst
  metafit <- rma(   yi=Est, sei=se, mods =X, data = ggdata)#Meta-Analysis via Linear (Mixed-Effects) Models
  metafit$pval[2]
  meta_res<-c( metafit$tau2  , metafit$b[2],metafit$pval[2])
  names(meta_res)<-c( 'tau^2' , 'est_parameter ' , ' p-value'  )
  print(   meta_res  )
}

#     tau^2    est_parameter     p-value  
#0.000000000    0.000835455    0.031582901  "dbp" 
#0.000000e+00  -2.191135e-03   1.154498e-05   "hip"  
#0.0000000000  -0.0008223616   0.6296620086 "monos"  
#0.0000000000  -0.0008337395   0.0017211867  "wt"
#1.980148e-06  -1.632990e-04   7.190750e-01 "ages"
#0.000000000    0.009122741    0.002549079  "vitcap" 
#0.000000000   -0.002419915    0.006404262 "bmi"
#0.000000e+00   6.200019e-05   9.065071e-01 "ht"


###result4: (newFig_12)---------------------------------------------------------------------------------------------------------
#according to vdat predicted HTE (i.e. v_predict), refit HTE ~ M to visualize the results or for future guidance-----------------

predict_matrix2<-getPredict( ALLRES_real$RES,2     )# ALLRES_real$Predicts_test  一样的
length( predict_matrix2[,200] ) #47395
dim(vdat)#47395    32
names(vdat)

apply(vdat, 2, is.numeric)
fit<-lm(  predict_matrix2[,ncol(predict_matrix2)]~ as.matrix(vdat[,5:32])     ) #共28个covariates
summary(fit)


dt_data<-cbind(vdat[,5:32],predict_matrix2[,200] ) #decision tree data
apply(dt_data, 2, is.numeric)
names(dt_data)[29]<-'estHTE'

names(dt_data)[8]<-'hip circumference'
names(dt_data)[3]<-'diastolic blood pressure'
names(dt_data)[13]<-'leucocyte count'
names(dt_data)[7]<-'waist circumference'
names(dt_data)[28]<-'vital capacity'
names(dt_data)[25]<-'platelet count'
names(dt_data)[22]<-'monocyte count'
names(dt_data)[5]<-'weight'

library(rpart)
library(rpart.plot)

tree <- rpart(estHTE ~. , data = dt_data,control = rpart.control(minsplit  = 2,cp=0.005))#cp可以控制tree的精度
rpart.plot(tree)#共14 leaf/end-node 
#newFig12 650 500
#newnewFig12 680 450

###result5: (Fig10)---------------------------------------------------------------------------------------------------
#individual predicted effect curve changing along the number of Q trees-----------------------------------------------
RES_real<-ALLRES_real$RES
predict_matrix2<-getPredict( RES_real,2     )#2: test set predicts
plot(   1:ncol( predict_matrix2   )   ,  predict_matrix2[10,],  type='l', #ylim=c(0.1,1.2),
        xlab='Number of Q trees', ylab='Predicted effect')
lines(     1:ncol( predict_matrix2   )   ,  predict_matrix2[20,], col='blue'   )
lines(     1:ncol( predict_matrix2   )   ,  predict_matrix2[30,], col='red'   )
lines(     1:ncol( predict_matrix2   )   ,  predict_matrix2[40,], col='green'   )

legend('topright', legend=c('Individual One', 'Individual Two', 'Individual Three', 'Individual Four'),
       col=c("black", "blue",'red','green'), lty=1, cex=0.7)

#ggplot version
ggdata<-data.frame( trees= 1:ncol( predict_matrix2)  , pe=c(predict_matrix2[10,],predict_matrix2[20,],predict_matrix2[30,],predict_matrix2[40,]),   
                    type=  factor(c( rep( 'Individual One' , ncol(predict_matrix2)  ) ,
                                     rep( 'Individual Two' , ncol(predict_matrix2)  ) ,
                                     rep( 'Individual Three' , ncol(predict_matrix2)  ) ,
                                     rep( 'Individual Four' , ncol(predict_matrix2)  ) )   , 
                                  levels=c('Individual One','Individual Two','Individual Three','Individual Four'  )  ) )
p<-ggplot(ggdata,aes(trees, pe,colour=type) ) +
  geom_line( linewidth=1.0 ) +
  xlab('Number of Q trees')+ylab('Predicted effect')+
  #geom_vline(xintercept = 200, colour='grey', linewidth=1.0 ,linetype = 2)+
  guides(colour = guide_legend(title = " "))+
  theme(panel.background = element_rect(linewidth=1,colour = "black"),
        legend.position = c(0.82, 0.78),legend.title = element_blank(),
        legend.background = element_rect( linewidth=0.5, linetype="solid",colour = "black"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p
ggsave(paste0('Fig10.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots'), 
       height = 5, width = 7, units = "in",limitsize=TRUE)


###result6: (Figure S1)----------------------------------------------------------------------------------------------------------
#permutation test results--------------------------------------------------------------------------------------------------------

#single >ALLRES_real<-RFQTfit(odat,vdat,Nb=200,method='DR') with cores=7 是一个半小时


ALLRES_real$TS
#0.003159107 165.9386


#请确保是real data的odat! 不要是simulation的odat
dim(odat) #94791    32  = I Z X Y + 28
dim(vdat) #47395    32  = I Z X Y + 28

TTS<-read.csv(  paste0('D:\\files\\R new\\Precison_Medicine\\real_data\\TTS.csv'),header=T, na.strings ="?")
nrow(  TTS )
#TTS<-c() #如果是半途再跑，小心不要设空！


for(np in ((nrow(TTS)+1):500)){ #pb: permutation+bootstrap  #np in 1:Np
  
  #针对odat permutate 一次
  Mpart<-odat[,5:ncol(odat)]
  set.seed(np)
  Mpart_p<-Mpart[ sample( 1:nrow(odat),nrow(odat)  ) ,  ]  #permutation 一次
  odat_p<-cbind(  odat[, 1:4],  Mpart_p )    
  #odat 必须要是一个data.frame  #is.data.frame(odat)  #TRUE 没毛病！
  
  ###same RFQT fitting for odat 
  #注意不需要vdat，因为ts1 ts2和vdat乃至OOBdat都没关系
  ALLRES_pb<-RFQTfit(odat_p,Nb=200,method='DR')

  TTS<-rbind(TTS ,ALLRES_pb$TS   )

  write.csv(TTS, paste0('D:\\files\\R new\\Precison_Medicine\\real_data\\TTS.csv'), row.names=F)
  print(  nrow(  TTS ) )
  print(  as.numeric( ALLRES_pb$TS  )  )
}
nrow(  TTS )
###nonparametric testing analysis - Kernal density estimation
bw.nrd0( TTS[,1] )  
plot(density(  TTS[,1]  ),main=expression(paste("Permutation test statistic  ", S[1], " values")) , xlim=c(0.002, 0.005))
abline(v= 0.003159107,col='red' ) #FigS1_1 500*400

bw.nrd0(TTS[,2] )  # 4.547092  3.867533
plot(density( TTS[,2]  ),main=expression(paste("Permutation test statistic  ", S[2], " values")),xlim=c(  115, 190 ))
abline(v= 165.9386,col='red' ) #FigS1_2 500*400

#其实没必要依据Kernel smoohting结果来算p-value，因为这会包含smoothing的error
#直接用nonparameteric! review: Kernel本身就是一种parametric！


Bernoulli_ts1<-(TTS[,1]>=0.003159107)
Bernoulli_ts2<-(TTS[,2]>=165.9386)

#由于样本数不够大！不够接近continuous的Normal distribution！还是用GLM吧
myf<-function(x){   exp(x)/(1+ exp(x)    ) }
fit1<-summary(glm(Bernoulli_ts1 ~ 1, family =binomial(link = "logit") ))

myf(  coef(fit1)[1] ) ;   myf(  coef(fit1)[1]+c(-1,1)*1.96*coef(fit1)[2] )   
#0.06923077     
#0.0364122 0.1277086   #130 





###*if HPC available; try:

data0<-read.csv(  'C:\\Users\\Haodong\\Desktop\\newPER.txt' , header = FALSE, sep = "", dec = ".")

okID<-sort(data0[,1])
(1:1000)[!(1:1000)%in%okID]#the killed array task ID
#741 742 901 902

#HPC补跑并把新的结果加进newPER,txt后

PERdata<-read.csv(  'C:\\Users\\Haodong\\Documents\\RFQT\\RFQT\\Data\\newPER.txt'   , header = FALSE, sep = "", dec = ".")
dim(PERdata)

okID<-sort(PERdata[,1])
(1:1000)[!(1:1000)%in%okID] #integer(0)

PERdata<-PERdata[order(PERdata[,1]),]

#same code as above
TTS<-PERdata[,-1]
  
nrow(  TTS )
###nonparametric testing analysis - Kernel density estimation
bw.nrd0( TTS[,1] )  
plot(density(  TTS[,1],from=0.001 ,to=0.005  ),  #truncate the extreme values
     main=expression(paste("Permutation test statistic ", S[1], " values")) , xlim=c(0.002, 0.005))
abline(v= 0.003159107,col='red' ) #FigS1_1 500*400

bw.nrd0(TTS[,2] )  # 4.547092  3.867533
plot(density( TTS[,2]  ),main=expression(paste("Permutation test statistic ", S[2], " values")),xlim=c(  115, 190 ))
abline(v= 165.9386,col='red' ) #FigS1_2 500*400

#其实没必要依据Kernel smoohting结果来算p-value，因为这会包含smoothing的error
#直接用nonparameteric! review: Kernel本身就是一种parametric！


Bernoulli_ts1<-(TTS[,1]>=0.003159107)
Bernoulli_ts2<-(TTS[,2]>=165.9386)

#由于样本数不够大！不够接近continuous的Normal distribution！还是用GLM吧

#for S1:
myf<-function(x){   exp(x)/(1+ exp(x)    )  }
fit1<-summary(glm(Bernoulli_ts1 ~ 1, family =binomial(link = "logit") ))

myf(  coef(fit1)[1] ) ;   myf(  coef(fit1)[1]+c(-1,1)*1.96*coef(fit1)[2] )   
#0.058     
#0.04510002 0.07430265   #1000 


#for S2:
c(0,3/1000)









