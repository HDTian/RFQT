#���ϵ�RFQT������ �����˲�������,ֻ������real data

RFQTfit<-function(odat,vdat,
                   Nb=1,  # the number of Q-trees
                   S=5,
                   rate=0.4,
                   Qthreshold=3.0,
                   method='DR',
                   SoP=10,
                   howGX='SpecificGX',
                   endsize=5000,
                   const=NA){
  #data check
  if( is.null(odat$true_STE[1]) ){JJ<-ncol( odat )-4}else{JJ<-ncol( odat )-5}
  if( JJ<1 ){ stop('No candidate covariate, or the dat is not regonized')  }
  #parameter define
  my.S<-S
  my.rate<-rate
  my.Qthreshold<-Qthreshold
  my.method<-method
  my.SoP<-SoP
  my.howGX<-howGX
  my.endsize<-endsize
  my.const<-const
  results<-c( JJ, nrow(odat), nrow(vdat) , method , SoP , rate, S , howGX, const, endsize,Qthreshold )
  print('Below is the summary of the parameters used for RFQT fitting')
  names(results)<-c( 'number.of.candidates.covariate',
                     'training.data.size',
                     'testing.data.size',
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
  cl<-makeCluster(n_cores_used)#�涨�ö��ٲ���clusters/�߳�
  clusterEvalQ(cl=cl , expr=library(dplyr))  #������cluster�е�����һЩ����ʽexpression (���絼��һЩ������һЩ��������)
  clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
  user_BootstrapTreeFitting<-function(seed){
    RES<-BootstrapTreeFitting(seed,  #��������û���⣬function�ڲ���my.rateȷʵ����RFQTreal����Լ������argument value(����Ĭ��ֵ)
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
                                      'GetTree', 'GetNindex', 'GetIndex' , 'BootstrapTreeFitting')  )
  RES<-parSapply(   cl ,  1:Nb, user_BootstrapTreeFitting  )
  
  ###results
  ALLRES<-list()
  
  ALLRES$RES<-RES
  
  ALLRES$MSE<-getMSE(RES,2)  #1 for OOB; 2 for testing
  ALLRES$Predicts<-getPredict(RES,2)  #��check�� MSE��Predict��ƥ���
  ALLRES$VI<-getVI(RES)
  
  return(ALLRES)
}