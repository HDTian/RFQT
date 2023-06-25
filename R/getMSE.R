#RES为parSapply(   cl ,  1:100, BootstrapTreeFitting  )的结果

#getMSE: 基于RES，计算MSE (only applicable for simulatin case where the label is known)

getMSE<-function(RES,#RES为parSapply(   cl ,  1:100, BootstrapTreeFitting  )的结果
                 indicator=1#1: for OOB error    2: for test error
                 ){
  if(is.null(  dim(RES)  )){ #即， RES只是一次输出结果，不是Sapply的结果
      if(  indicator==2    ){
        #先判断是否有vdat
        if(  is.na(  RES$v_predict )){stop('No testing set used, consider OOB error only')}
        #再判断vdat是否有STE
        if(is.null(vdat$true_STE)){stop('Note: no individual treatment effect was defined. MSE values are not applicable.')}
        MSE<-mean( (RES$v_predict   - vdat$true_STE)^2  )
      }else{MSE<-mean(   ( (RES$OOB_predict  - odat$true_STE)[RES$OOB_predict!=0  ] )^2  )}
    }else{
      BN<-dim(RES)[2]  #dim(RES): 7 100
      #(7为输出的7种$结果: $end_node_information $OOB_predict $v_predict $vi1 $vi2 $ts1 $ts2;
      #100为bootstrapped次数/即RFQT的size/即Nb)
      if(  indicator==2    ){#for test error
        #先判断是否有vdat
        if(  is.na(  RES[3,1]$v_predict[1] )){stop('No testing set used, consider OOB error only')}
        #再判断vdat是否有STE
        if(is.null(vdat$true_STE)){stop('Note: no individual treatment effect was defined. MSE values are not applicable.')}
        v_predict<-c()
        for(i in 1:BN){
          v_predict<-cbind( v_predict ,  RES[3,i]$v_predict )  #v_predict
        }
        total_v_predict<-apply(v_predict, 1  , cumsum )
        v_MSE<-c()
        for(i in 1:BN){
          v_MSE<-c(     v_MSE,   mean((total_v_predict[i,]/i  - vdat$true_STE )^2)  )  #注意apply(v_predict, 1  , cumsum )的结果会自动变成转置矩阵！
        }
        MSE<-v_MSE
      }else{
        oob_predict<-c()
        oob_times<-c()
        for(i in 1:BN){
          oob_predict<-cbind( oob_predict ,  RES[2,i]$OOB_predict )  #OOB_predict
          oob_times<-cbind( oob_times ,    as.numeric(  abs( RES[2,i]$OOB_predict-0) > .Machine$double.eps^0.5 )  )
        }
        total_oob_predict<-apply(oob_predict, 1  , cumsum )
        total_oob_times<-apply(oob_times, 1  , cumsum )
        oob_MSE<-c()
        for(i in 1:BN){
          oob_MSE<-c(     oob_MSE,   mean(   (           (total_oob_predict[i,]/total_oob_times[i,]  - odat$true_STE)[total_oob_times[i,]!=0  ]     )^2   )  )
        }
        MSE<-oob_MSE
      }
  }

  return(MSE)
}
#MSE: 一个vector; length = Nb; 表示不同Nb下(即，随着tree的数量的增长)的MSE (iNdicator either for OOB or testing)

###examples:
#RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  )

getMSE(RES,indicator=1 )
getMSE(RES,indicator=2 )
