
#这个函数其实更泛用一些

getPredict<-function(RES,#RES为parSapply(   cl ,  1:100, BootstrapTreeFitting  )的结果  ;
                     indicator=1 #1: for OOB predicts    2: for test predicts
                     ){
  if(is.null(  dim(RES)  )){ #即， RES只是一次输出结果，不是Sapply的结果
    if(  indicator==2    ){
      if(  is.na(  RES$v_predict )){stop('No testing set used, consider OOB error only')}
      Predict<-RES$v_predict
    }else{Predict<- RES$OOB_predict  ;Predict[RES$OOB_predict==0  ]<-NA    }
  }else{
    BN<-dim(RES)[2]  #dim(RES): 7 100 (7为输出的7种$结果: $end_node_information $OOB_predict $v_predict $vi1 $vi2 $ts1 $ts2;
    #100为boostrapped次数/即RFQT的size/即BN)
    if(  indicator==2    ){
      if(  is.na(  RES[3,1]$v_predict[1] )){stop('No testing set used, consider OOB error only')}
      v_predict<-c()
      for(i in 1:BN){
        v_predict<-cbind( v_predict ,  RES[3,i]$v_predict )  #v_predict
      }
      total_v_predict<-apply(v_predict, 1  , cumsum )  #注意：还是一个矩阵噢
      v_Predict<-c()
      for(i in 1:BN){
        v_Predict<-cbind(     v_Predict,   total_v_predict[i,]/i  ) #注意apply(v_predict, 1  , cumsum )的结果会自动变成转置矩阵！
      }
      Predict<-v_Predict #不同于之前的MSE 结果为一个vector；此时Predict结果自然为一个矩阵！
    }else{
      oob_predict<-c()
      oob_times<-c()
      for(i in 1:BN){
        oob_predict<-cbind( oob_predict ,  RES[2,i]$OOB_predict )  #OOB_predict
        oob_times<-cbind( oob_times ,    as.numeric(  abs( RES[2,i]$OOB_predict-0) > .Machine$double.eps^0.5 )  )
      }
      total_oob_predict<-apply(oob_predict, 1  , cumsum )
      total_oob_times<-apply(oob_times, 1  , cumsum )
      oob_Predict<-c()
      for(i in 1:BN){
        oob_Predict<-cbind(     oob_Predict,              total_oob_predict[i,]/total_oob_times[i,]           )
        oob_Predict[total_oob_times[i,]==0]<-NA
      }
      Predict<-oob_Predict #不同于之前的MSE 结果为一个vector；此时Predict结果自然为一个矩阵！
    }
  }






  return(Predict)  #Predict是一个矩阵！nrow=vdat sample size; ncol=Bootstrap的次数/即size of RFQT/即BN/即Nb
  #往往可以看某个individual的predict value的随bootstrap次数增加而变化的稳定性来决定Nb是否选取的足够合适大了！
  #最后的RFQT predict结果可以取矩阵的最后一列结果(不用再cumsum操作了！)
}

#getPredict returns a matrix (individual number * Nb): 表示随着tree的增长目前的forest下的各个individual的predicted effects (即HTE)


# ###examples:
# #RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  )
# 
# predict_matrix<-getPredict(RES,indicator=1 )
# predict_matrix<-getPredict(RES,indicator=2 )
# 
# dim(predict_matrix)
# View(predict_matrix)
