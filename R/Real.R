


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
ALLRES_real<-RFQTfit(odat,vdat,Nb=200,method='DR')\

#saveRDS(ALLRES_real,file='D:\\files\\R new\\Precison_Medicine\\real_data\\ALLRES_real.RData')
#ALLRES_real<-readRDS('D:\\files\\R new\\Precison_Medicine\\real_data\\ALLRES_real.RData')


###result1: (Figure 4)
#histogram plot of the testing set predicts using (DR Nb=200)
summary(  ALLRES_real$Predicts_test[,200]  )
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.717172 -0.028483 -0.019300 -0.020037 -0.009268  0.531695 


hist( ALLRES_real$Predicts_test[,200]  , n=500, xlim=c(-0.1,0.05) , main='', xlab='Predicted effect')


###results2: (Figure S2)
#Variable importance plot



###result3: (Figure 5)
#important variable -> marginal stratum analysis using pre-detetmined quantile points
#one can use DRMR package
#Q statistic use average; point and CIS use Rubin's Rule



###result4: (new Figure)
#according to vdat predicted HTE (i.e. v_predict), refit HTE ~ M to visualize the results or for future guidance



###result5: (Figure S1)
#permutation test results




















