#RES为parSapply(   cl ,  1:100, BootstrapTreeFitting  )的结果

getMSE<-function(RES,
                 indicator=1){  #RES为parSapply(   cl ,  1:100, BootstrapTreeFitting  )的结果  ;
  #1: for OOB error    2: for test error
  if(is.null(vdat$true_STE)){print('Note: no individual treatment effect was defined. MSE values are meaningless.')}
  if(is.null(  dim(RES)  )){ #即， RES只是一次输出结果，不是Sapply的结果
    if(  indicator==2    ){
      if(  is.na(  RES$v_predict )){stop('No testing set used, consider OOB error only')}
      MSE<-mean( (RES$v_predict   - vdat$true_STE)^2  )
    }else{MSE<-mean(   ( (RES$OOB_predict  - odat$true_STE)[RES$OOB_predict!=0  ] )^2  )}
    }else{
      BN<-dim(RES)[2]  #dim(RES): 4 100
      #(4为输出的4种$结果: $OOB_predict $v_predict $vi1 $vi2 ; 100为bootstrapped次数/即RFQT的size/即BN)
      if(  indicator==2    ){
        if(  is.na(  RES[2,1]$v_predict[1] )){stop('No testing set used, consider OOB error only')}
        v_predict<-c()
        for(i in 1:BN){
          v_predict<-cbind( v_predict ,  RES[2,i]$v_predict )  #v_predict
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
          oob_predict<-cbind( oob_predict ,  RES[1,i]$OOB_predict )  #OOB_predict
          oob_times<-cbind( oob_times ,    as.numeric(  abs( RES[1,i]$OOB_predict-0) > .Machine$double.eps^0.5 )  )
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
#MSE: 一个vector; length = Nb; 表示不同Nb下的MSE (iNdicator either for OOB or testing)
