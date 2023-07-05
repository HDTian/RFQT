


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


###result1: (Figure 4)
#histogram plot of the testing set predicts using (DR Nb=200)
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


###results2: (Figure S2)
#Variable importance plot
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


###result3: (Figure 5)
#important variable -> marginal stratum analysis using pre-determined (tenth) quantile points
#one can use DRMR package
#Q statistic use average; point and CIS use Rubin's Rule (RR)?

VInames<-(names(Dat)[VIindex+4])[length(VIindex):1]
VInames[1:6]  #"dbp"    "hip"    "monos"  "wt"     "ages"   "vitcap"

dim(Dat)
names(Dat)
#devtools::install_github("HDTian/DRMR")
library(DRMR)

for(i in 1:6){
  Mname<-VInames[i]  #M: marginal specific covariate name
  dat<-Dat[,     c(   2 ,3,4, match(Mname,names(Dat)   ) ) ]
  names(dat)<-c(  'Z', 'X', 'Y', 'M' ) #for DRMR
  apply(dat, 2, is.numeric) #TRUE TRUE TRUE TRUE
  
  #stratification via DRMR
  rdat<-Stratify(dat,onExposure = FALSE)#augmented dataset
  SRES<-getSummaryInf(rdat) #Stratification RESults
  
  #visualization
}
dat<-DAT




###result4: (new Figure)
#according to vdat predicted HTE (i.e. v_predict), refit HTE ~ M to visualize the results or for future guidance




###result5: (Figure S1)
#permutation test results




















