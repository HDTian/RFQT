#RES为parSapply(   cl ,  1:100, BootstrapTreeFitting  )的结果

#getMSE: 基于RES，计算MSE (only applicable for simulatin case where the label is known)
#如果是real data等没有HTE label的情况，则return ’no_label‘

getMSE<-function(RES,#RES为parSapply(   cl ,  1:100, BootstrapTreeFitting  )的结果
                 indicator=1#1: for OOB error    2: for test error
){
  if(is.null(  dim(RES)  )){ #即， RES只是一次输出结果，不是Sapply的结果
    
    if(indicator==1){#OOB MSE case
      #OOB sample一定存在
      #判断odat是否有STE
      if(is.null(odat$true_STE)){
        print('Note: no individual treatment effect for odat was defined. MSE values are not applicable.')
        MSE<-'no_label'
      }else{ MSE<-mean(   ( (RES$OOB_predict  - odat$true_STE)[RES$OOB_predict!=0  ] )^2  )  }
    }
    
    if(indicator==2){#test MSE case
      #先判断是否有vdat
      if(  is.na(  RES$v_predict[1] )){stop('No testing set used, consider OOB error only')}
      #再判断vdat是否有STE(此时一定存在vdat)
      if(is.null(vdat$true_STE)){
        print('Note: no individual treatment effect for vdat was defined. MSE values are not applicable.')
        MSE<-'no_label'
      }else{ MSE<-mean( (RES$v_predict   - vdat$true_STE)^2  ) }
    }
    
    
  }else{  #即，此时RES是Sapply的结果
    BN<-dim(RES)[2]  #dim(RES): 7 100
    #(7为输出的7种$结果: $end_node_information $OOB_predict $v_predict $vi1 $vi2 $ts1 $ts2;
    #100为bootstrapped次数/即RFQT的size/即Nb)
    
    if(indicator==1){#OOB MSE case
      #OOB sample一定存在
      #判断odat是否有STE
      if(is.null(odat$true_STE)){
        print('Note: no individual treatment effect for odat was defined. MSE values are not applicable.')
        MSE<-'no_label'
      }else{
        ### OOB MSE --------------------------------------
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
        ### -----------------------------------------------
      }
    }
    
    if(indicator==2){#test MSE case
      #先判断是否有vdat
      if(  is.na(  RES[3,1]$v_predict[1] )){stop('No testing set used, consider OOB error only')}
      #再判断vdat是否有STE(此时一定存在vdat)
      if(is.null(vdat$true_STE)){
        print('Note: no individual treatment effect for vdat was defined. MSE values are not applicable.')
        MSE<-'no_label'
      }else{
        ### test MSE --------------------------------------
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
        ### -----------------------------------------------
      }
    }
    
  }
  
  return(MSE)
}
#在存在labEl时,MSE是一个vector; length = Nb; 表示不同Nb下(即，随着tree的数量的增长)的MSE (iNdicator either for OOB or testing)
#不存在label时,MSE是一个string : 'no_label'

###examples:

#RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  )
#or
#RES<-ALLRES$RES
set.seed(60)
res<-getDat() #simulated data  #the deflaut setting: scenario='A' and SoM=0.5
res<-getDat(label=FALSE) 
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set

#When running RFQT with mutiple Q trees/bootstrap - use parallel computation
Nb<-5 #how many trees in the forest? e.g. 5 trees
cl<-makeCluster(detectCores()-1)#规定用多少并行clusters/线程
clusterEvalQ(cl=cl , expr=library(dplyr))  #给各个cluster中的运行一些表达式expression (例如导入一些包；做一些基础操作)
clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
clusterExport(  cl=cl ,  varlist=c('odat','vdat',#只有一个环境：全局环境
                                   'GetTree', 'GetNindex', 'GetIndex' )  )#or any other arguments
RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  ) #本脚本里parSapply的结果常用RES表示
stopCluster(cl)
dim(RES)#7 Nb


getMSE(RES,indicator=1 )
getMSE(RES,indicator=2 )



