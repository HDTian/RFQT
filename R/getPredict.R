
#���������ʵ������һЩ

getPredict<-function(RES,indicator=1){  #RESΪparSapply(   cl ,  1:100, BootstrapTreeFitting  )�Ľ��  ;
  #1: for OOB predicts    2: for test predicts
  if(is.null(  dim(RES)  )){ #���� RESֻ��һ��������������Sapply�Ľ��
    if(  indicator==2    ){
      if(  is.na(  RES$v_predict )){stop('No testing set used, consider OOB error only')}
      Predict<-RES$v_predict
    }else{Predict<- RES$OOB_predict  ;Predict[RES$OOB_predict==0  ]<-NA    }
  }else{
    BN<-dim(RES)[2]  #dim(RES): 5 100 (6Ϊ�����6��$���: $OOB_predict  $v_predict $v_Epredict   $vi2 $ts1  $ts2 ;
    #100Ϊboostrapped����/��RFQT��size/��BN)
    if(  indicator==2    ){
      if(  is.na(  RES[2,1]$v_predict[1] )){stop('No testing set used, consider OOB error only')}
      v_predict<-c()
      for(i in 1:BN){
        v_predict<-cbind( v_predict ,  RES[2,i]$v_predict )  #v_predict
      }
      total_v_predict<-apply(v_predict, 1  , cumsum )  #ע�⣺����һ��������
      v_Predict<-c()
      for(i in 1:BN){
        v_Predict<-cbind(     v_Predict,   total_v_predict[i,]/i  ) #ע��apply(v_predict, 1  , cumsum )�Ľ�����Զ����ת�þ���
      }
      Predict<-v_Predict #��ͬ��֮ǰ��MSE ���Ϊһ��vector����ʱPredict�����ȻΪһ������
    }else{
      oob_predict<-c()
      oob_times<-c()
      for(i in 1:BN){
        oob_predict<-cbind( oob_predict ,  RES[1,i]$OOB_predict )  #OOB_predict
        oob_times<-cbind( oob_times ,    as.numeric(  abs( RES[1,i]$OOB_predict-0) > .Machine$double.eps^0.5 )  )
      }
      total_oob_predict<-apply(oob_predict, 1  , cumsum )
      total_oob_times<-apply(oob_times, 1  , cumsum )
      oob_Predict<-c()
      for(i in 1:BN){
        oob_Predict<-cbind(     oob_Predict,              total_oob_predict[i,]/total_oob_times[i,]           )
        oob_Predict[total_oob_times[i,]==0]<-NA
      }
      Predict<-oob_Predict #��ͬ��֮ǰ��MSE ���Ϊһ��vector����ʱPredict�����ȻΪһ������
    }
  }





  return(Predict)  #Predict��һ������nrow=vdat sample size; ncol=Bootstrap�Ĵ���/��size of RFQT/��BN/��Nb
  #�������Կ�ĳ��individual��predict value����bootstrap�������Ӷ��仯���ȶ���������Nb�Ƿ�ѡȡ���㹻���ʴ��ˣ�
  #����RFQT predict�������ȡ��������һ�н��(������cumsum�����ˣ�)
}