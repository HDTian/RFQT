

#一系列的functions for HDTian/RFQT 最后打包成一个package

#这些.R scripts的存储位置不重要；因为最后package会有特定的目???

#关于中间需要使用的packages，直接在这里???.R scripts里写library(...)就好


#encoding: UTF-8 否则中文全部乱码

library(MendelianRandomization)#naive IVW
library(tidyverse) #数据处理??? DR 需???
library(data.table) #数据处理??? DR 需???

library(parallel) #并行运算

library(corrplot) #data cleaning需要，RFQT可以不需???



rm(list=ls())



###my function order:
#getDat  (GetIndex->GetTree->GetNindex->BootstrapTreeFitting) (getMSE getPredict getVI)





######下面的内容是模拟使用这个包的人：(用在网页readme介绍)

###simulation部分
set.seed(60)
res<-getDat() #simulated data
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set      #确实和源代码是一样的

#RFQT

#method<-'DR'; SoP<-20; rate<-2/5
getPars() #check the parameters

#fit a single Q-tree
tree_result<-GetTree(odat) #return a augmented data of the training data containing individual stratification information within the Q tree

#or use RFQT with one-time bootstrap
RES<-BootstrapTreeFitting(seed=1,odat)
#$OOB_predict   $v_predict  $vi1  $vi2


#When running RFQT with mutiple Q trees/bootstrap
n_cores_used<-detectCores()-1
Nb<-7
cl<-makeCluster(n_cores_used)#规定用多少并行clusters/线程
clusterEvalQ(cl=cl , expr=library(dplyr))  #给各个cluster中的运行一些表达式expression (例如导入一些包；做一些基础操作)
clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat',
                                    'GetTree', 'GetNindex', 'GetIndex' )  )#or any other arguments
RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  ) #本脚本里parSapply的结果常用RES表示
stopCluster(cl)

#If you wish to use your own parameters rather than the default parameters, try
user_BootstrapTreeFitting<-function(seed){
  RES<-BootstrapTreeFitting(seed,
                            S=my.S,
                            JJ=my.JJ,   #or any partial of them
                            rate=my.rate,
                            Qthreshold=my.Qthreshold,
                            method=my.method,
                            SoP=my.SoP,
                            howGX=my.howGX,
                            endsize=my.endsize)
  return(RES)
}

user_BootstrapTreeFitting<-function(seed){
  RES<-BootstrapTreeFitting(seed,SoP=20)
  return(RES)
}
clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat',
                                    'GetTree', 'GetNindex', 'GetIndex' , 'BootstrapTreeFitting')  )
RES<-parSapply(   cl ,  1:Nb, user_BootstrapTreeFitting  ) #本脚本里parSapply的结果常用RES表示
stopCluster(cl)
#已经验证！和原代码是一摸一样的!!! - updated/check data: 2023/3/24
#updated/check data: 2023/4/3: 修改了stratification logic后，
#结果没有完全复现是因为原始代码中的rank( dat_order$residual,ties.method ='random' )这个random 引入了新的randomness ；
#导致rate=/=1时，每次的candidate covariate set不一样;这其实没太大问题
#将rate==1时，结果完美复现了
#注意：DR由于rank的性质，自然无法完全复现; (exclusion variability)


#Result analysis
getMSE(RES,2)  #1 for OOB; 2 for testing
getPredict(RES,2)  #已check??? MSE和Predict是匹配的
getVI(RES)

#One function
allres<-RFQTfit(odat,vdat,Nb=5,SoP=20)

###real data part--------------------------------------------------------------

#real data多出来了一个quality control: (1)检测是否合规且课识??? (2)removing missing individuals



path<-'D:\\files\\Paper_Rcode'  #please use your path
master_cut<-read.csv('D:\\files\\real data\\master_cut.csv')
bmi_vars<-read.csv('D:\\files\\real data\\bmi_vars.csv')
#data.table起手操作
dtp<- as.data.table(master_cut)#dt phenotype: 各种variables; UKB_sample_ID为UKB ID
dtg<- as.data.table(bmi_vars)#dt genotype:和BMI相关的SNPs；SNPs(或者Z/gene score); X为UKB ID
nSNPs<-dim(dtg)[2]-1
#initial selection
dtp<-dtp[dtp$euro==1,]
dtp<-dtp[dtp$race=='White',]
dtp<-dtp[dtp$sex=='Male',]
dim(dtp)  # 197128  for Female ; 167121  For Male
dtpp<- dtp[ , c(3,   48,     56 ,     7, 44:55,57:61, 73:82) ] #BMI on fev1
names(dtpp)#确认一下有无问???
dim(dtpp) #488366     31
dt<-dtpp[dtg, on = .(UKB_sample_ID = X)] #按照UKB_sample_ID和X这两个variable相等条件合并phenotype和genotype两个???
dim(dt)#367643    127   #127=31+ 96 (???96==nSNPs个BMI-related SNPs)
DT<-na.omit(dt)
index<-(   (  dim(DT)[2]-nSNPs+1    )  :   (dim(DT)[2])      )
DTT<-round(    DT[,   ..index  ]      )  #dim(DTT)  #168274          96  #names(DTT) 没问???
###remove the correlated SNPs and then build GXM
GXM<-as.matrix(      cbind(   DTT , DT[ , 2]     )             )#这次加上X
GXM<-GXM[,-c(2,77)]
dim(GXM)#  168274       95= 94 SNPs + BMI
#single instrument (Z/gene score) building  ###也提供一下有外部dataset的genetic score；因为正文中是说用的第三方weights for gene score
posi<-ncol(GXM)
fitGX<-lm( GXM[,posi ]~ GXM[,-posi]   ) #SNPs参与的数量没问题
Z<- as.numeric(   fitted( fitGX   ) )
Z<- Z -mean(Z)  #centrized
###now formally build the data with the same style as the sim
index2<-(   2:(  dim(DT)[2]-nSNPs    )    ) #DT是没有Missing data的data
Dat<-cbind(   (1:dim(DT)[1]) ,Z,  DT[ ,..index2]  )  #Dat<-cbind( 1:NNNN , Z , X ,Y ,MM  )
dim(Dat)  # 313329     33
names(Dat)  #index + Z + BMI(exposure) + fev1(outcome) + Ms
fitX<-lm( lm(  Dat$bmi ~ Dat$Z )     )
RR<-as.numeric(summary(fitX)$r.squared[1])#RR 0.01911255  意料之中 (past paper with other sample size is 0.017)
RR  #0.01894159
#先打乱一下，以防出现UKB_ID和其他variable存在关联???
set.seed(1123)
index<-sample(     1:  dim(Dat)[1] )
Dat<-Dat[ index,   ]
Dat$V1<-1:  dim(Dat)[1]  #重新安排individual index
#规定好training test data size; 并且training - testing data size ratio 大致???2:1
NNNN<-dim(Dat)[1]; NNN<-round(NNNN*2/3)    ; endsize<-5000 #5000应该可以???
NNNN;NNN#for female: 168274(NNNN)  112183(NNN)   ;For male: 142186(NNNN) 94791(NNN)
#同一个Dat下，odat和vdat永远不会改变
vdat<-Dat[(NNN+1):NNNN, ]  #validation dat
odat<-Dat[1:NNN,]   #original (training) dat
odat<-as.data.frame(odat);vdat<-as.data.frame(vdat);Dat<-as.data.frame(Dat)
names(odat)[1]<-'I'
names(odat)[3]<-'X'
names(odat)[4]<-'Y'

dim(vdat); dim(odat)
names(vdat)
names(odat)
#odat, vdat一直是外部环境

#Result analysis
getMSE(RES,2)  #1 for OOB; 2 for testing
getPredict(RES,2)  #已check??? MSE和Predict是匹配的
getVI(RES)

allres<-RFQTfit(odat,vdat,Nb=5,SoP=20)

#The variable ordering (the end variables are more important)
(names(odat)[-(1:4)] )[order(  allres$VI )]
