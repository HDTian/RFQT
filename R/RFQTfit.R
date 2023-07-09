###one single function to run RFQT

RFQTfit<-function(odat, #training set
                  vdat=NA, #validation set (can be empty)
                  Nb=5,  # the number of Q-trees
                  S=5, #the largest depth
                  honest=FALSE, #use honest estimation or not?
                  rate=0.4,# the proportion of candidate variables Ms considered
                  SingleM=FALSE,#whether to ue single stratification style?(i,e, use a fixed M determined in the begining)
                  Qthreshold=3.84, ##the threshold for Q heterogneity assessment
                  method='DR',#stratification method used: 'DR' 'Residual' others
                  SoP=10,##size of pre-stratum #only make sense to DR stratification
                  howGX='SpecificGX',##'const' means use extra constant; otherwise estimated by stratum data (stratum-specific GXeffect)
                  endsize=1000,#the minimal size of the node of Q-tree allowed to exist
                  const=NA, #the pre-given fixed GX effect #only make sense when howGX='const'
                  Cores=NA, #how many cores to be used
                  trackfile=NA #if to track the seed results? if TRUE, need to specify the redirector path (default is no path file to track; i.e. NA)
){
  ###data check
  if( is.null(odat$true_STE[1]) ){JJ<-ncol( odat )-4}else{JJ<-ncol( odat )-5}
  if( JJ<1 ){ stop('No candidate covariate, or the data is not regonized')  }
  
  ###parameter define
  my.honest<-honest
  my.S<-S
  my.rate<-rate
  my.SingleM<-SingleM
  my.Qthreshold<-Qthreshold
  my.method<-method
  my.SoP<-SoP
  my.howGX<-howGX
  my.endsize<-endsize
  my.const<-const
  if(as.matrix(is.na(vdat))[1,1]){  vdat_information<-'NA'  }else{ vdat_information<-nrow(vdat)}
  results<-c( JJ, nrow(odat), vdat_information , honest , method , SoP , rate, SingleM, S , howGX, const, endsize, Qthreshold )
  print('Below is the summary of the parameters used for RFQT fitting')
  names(results)<-c( 'number.of.candidates.covariate',
                     'training.data.size',
                     'testing.data.size',
                     'honest.or.not',
                     'stratification.method',
                     'size.of.pre-stratum',
                     'random.proportion',
                     'use.single.stratification',
                     'max.tree.deep' ,
                     'instrument-exposure.style',
                     'GXeffect.value',
                     'min.end.node.size',
                     'Q.value.threshold')
  print(results)
  
  #RFQT fitting
  if(is.na(Cores)){   
    n_cores_used<-detectCores()-1
  }else{
    n_cores_used<-Cores
    }
  print(paste0('The number of cores used: ',n_cores_used))
  
  if( is.na(trackfile)  ){
    cl<-makeCluster(n_cores_used) 
  }else{    cl<-makeCluster(n_cores_used,outfile=trackfile)      }

  
  clusterEvalQ(cl=cl , expr=library(dplyr))
  clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
  
  user_BootstrapTreeFitting<-function(seed){
    RES<-BootstrapTreeFitting(seed,
                              Odat=odat,#其实可以不用这一行，因为BootstrapTreeFitting的Odat是有默认值的
                              Vdat=vdat, #vdat需要去环境中找
                              honest=my.honest,#一旦my.honest找不到，就应该报错；和honest的默认值没有关系
                              S=my.S,
                              rate=my.rate,
                              SingleM=my.SingleM,
                              Qthreshold=my.Qthreshold,
                              method=my.method,
                              SoP=my.SoP,
                              howGX=my.howGX,
                              endsize=my.endsize,
                              const=my.const)
    return(RES)
  }
  #clusterExport里 varlist会去全局环境中匹配，而不是当前环境
  #clusterExport assigns the values on the master R process of the variables named in varlist to variables of the same names in the global environment (aka ‘workspace’) of each node.
  #The environment on the master from which variables are exported defaults to the global environment.
  # clusterExport(  cl=cl ,  varlist=c( 'odat',#in 全局环境
  #                                     'vdat',#会是默认的vdat=NA，不会是全局变量的vdat
  #                                     'my.honest','my.S','my.rate','my.SingleM','my.Qthreshold','my.method','my.SoP','my.howGX','my.endsize','my.const',
  #                                     'GetTree', 'GetNindex', 'GetIndex' , 'BootstrapTreeFitting')  )
  clusterExport(  cl=cl ,  varlist=c( 'GetTree', 'GetNindex', 'GetIndex' , 'BootstrapTreeFitting'))
  RES<-parSapply(   cl ,  1:Nb, user_BootstrapTreeFitting )
  stopCluster(cl)
  
  ###results--------------------------------------------------------------------
  ALLRES<-list()
  #RES: $end_node_information $OOB_predict $v_predict $vi1 $vi2 $ts1 $ts2
  ALLRES$RES<-RES
  
  #OOB and testing MSE
  if( as.matrix(is.na(vdat))[1,1] ){  ALLRES$MSE_OOB<-getMSE(RES,1)  #1 for OOB; 2 for testing
  }else{
    ALLRES$MSE_OOB<-getMSE(RES,1)  #注意，function内部的odat vdat是在全局环境中寻找的
    ALLRES$MSE_test<-getMSE(RES,2)
  }
  
  #OOB and testing set individual prediction
  if( as.matrix(is.na(vdat))[1,1] ){  ALLRES$Predicts_OOB<-getPredict(RES,1)  #1 for OOB; 2 for testing
  }else{
    ALLRES$Predicts_OOB<-getPredict(RES,1)  # #注意，function内部的odat vdat是在全局环境中寻找的
    ALLRES$Predicts_test<-getPredict(RES,2)
  }
  
  #variable importance(VI只和OOB有关系，和vdat无关)
  ALLRES$VI1<-getVI(RES,VItype=1)#vi1: label known
  ALLRES$VI2<-getVI(RES,VItype=2)#vi2: label unknown
  
  
  #permutation test statistics ts1 and ts2
  BN<-dim(RES)[2]
  ts_res<-c()
  for(i in 1:BN){
    ts_res<-rbind(ts_res ,   c(RES[6,i]$ts1, RES[7,i]$ts2 )  )
  }
  TS_res<-apply( ts_res , 2 , mean  )
  names(TS_res)<-c('ts1','ts2')
  ALLRES$TS<-TS_res
  
  
  return(ALLRES)
}

###example:
set.seed(60)
res<-getDat() #simulated data  #the deflaut setting: scenario='A' and SoM=0.5
odat<-res$traning.set  #training set in 全局环境
vdat<-res$testing.set  #testing set in 全局环境
#or: vdat <- NA

ALLRES<-RFQTfit(odat)

ALLRES<-RFQTfit(odat,vdat)

ALLRES<-RFQTfit(odat,SingleM=TRUE)#single stratification style



saveRDS(ALLRES,file='D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\test.RData')
ALLRES_<-readRDS('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\test.RData')

#Error in get(name, envir = envir) : object 'my.honest' not found 因为全局环境中没有my.honest这个量


