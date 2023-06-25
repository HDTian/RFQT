###one single function to run RFQT

RFQTfit<-function(odat, #training set
                  vdat=NA, #validation set (can be empty)
                   Nb=5,  # the number of Q-trees
                   S=5, #the largest depth
                   honest=TRUE, #use honest estimation or not?
                   rate=0.4,# the proportion of candidate variables Ms considered
                   Qthreshold=3.84, ##the threshold for Q heterogneity assessment
                   method='DR',#stratification method used
                   SoP=10,##size of pre-stratum #only make sense to DR stratification
                   howGX='SpecificGX',##'const' means use extra constant; otherwise estimated by stratum data (stratum-specific GXeffect)
                   endsize=5000,#the minimal size of the node of Q-tree allowed to exist
                   const=NA #the pre-given fixed GX effect #only make sense when howGX='const'
                  ){
  ###data check
  if( is.null(odat$true_STE[1]) ){JJ<-ncol( odat )-4}else{JJ<-ncol( odat )-5}
  if( JJ<1 ){ stop('No candidate covariate, or the data is not regonized')  }

  ###parameter define
  my.honest<-honest
  my.S<-S
  my.rate<-rate
  my.Qthreshold<-Qthreshold
  my.method<-method
  my.SoP<-SoP
  my.howGX<-howGX
  my.endsize<-endsize
  my.const<-const
  if(as.matrix(is.na(vdat))[1,1]){  vdat_information<-'NA'  }else{ vdat_information<-nrow(vdat)}
  results<-c( JJ, nrow(odat), vdat_information , honest , method , SoP , rate, S , howGX, const, endsize, Qthreshold )
  print('Below is the summary of the parameters used for RFQT fitting')
  names(results)<-c( 'number.of.candidates.covariate',
                     'training.data.size',
                     'testing.data.size',
                     'honest.or.not',
                     'stratification.method',
                     'size.of.pre-stratum',
                     'random.proportion',
                     'max.tree.deep' ,
                     'instrument-exposure.style',
                     'GXeffect.value',
                     'min.end.node.size',
                     'Q.value.threshold')
  print(results)

  #RFQT fitting
  n_cores_used<-detectCores()-1
  print(paste0('The number of cores used: ',n_cores_used))
  cl<-makeCluster(n_cores_used)
  clusterEvalQ(cl=cl , expr=library(dplyr))
  clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
  user_BootstrapTreeFitting<-function(seed){
    RES<-BootstrapTreeFitting(seed,
                              #Vdat=vdat, #一旦BootstrapTreeFitting内部找不到vdat，就会去全局环境中找
                              honest=my.honest,
                              S=my.S,
                              rate=my.rate,
                              Qthreshold=my.Qthreshold,
                              method=my.method,
                              SoP=my.SoP,
                              howGX=my.howGX,
                              endsize=my.endsize,
                              const=my.const)
    return(RES)
  }
  clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat',
                                      'my.honest','my.S','my.rate','my.Qthreshold','my.method','my.SoP','my.howGX','my.endsize','my.const',
                                      'GetTree', 'GetNindex', 'GetIndex' , 'BootstrapTreeFitting')  )
  RES<-parSapply(   cl ,  1:Nb, user_BootstrapTreeFitting  )
  stopCluster(cl)

  ###results
  ALLRES<-list()
  #RES: $end_node_information $OOB_predict $v_predict $vi1 $vi2 $ts1 $ts2
  ALLRES$RES<-RES

  #OOB and testing MSE
  if( is.na(vdat) ){  ALLRES$MSE_OOB<-getMSE(RES,1)  #1 for OOB; 2 for testing
  }else{
    ALLRES$MSE_OOB<-getMSE(RES,1)
    ALLRES$MSE_test<-getMSE(RES,2)
  }

  #OOB and testing set individual prediction
  if( is.na(vdat) ){  ALLRES$Predicts_OOB<-getPredict(RES,1)  #1 for OOB; 2 for testing
  }else{
    ALLRES$Predicts_OOB<-getPredict(RES,1)
    ALLRES$Predicts_test<-getPredict(RES,2)
  }

  #variable importance(VI只和OOB有关系，和vdat无关)
  ALLRES$VI1<-getVI(RES,VItype=1)#vi1: label known
  ALLRES$VI2<-getVI(RES,VItype=2)#vi2: label unknown

  return(ALLRES)
}

###example:
set.seed(60)
res<-getDat() #simulated data  #the deflaut setting: scenario='A' and SoM=0.5
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set

ALLRES<-RFQTfit(odat)

