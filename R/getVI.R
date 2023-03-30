

#只返回各个candidate covariate的VI measurements;不用表示成order顺序了

getVI<-function(RES,
                VItype=2){  #get平均的VI measurement
  if(is.null(  dim(RES)  )){ #即， RES只是一次输出结果，不是Sapply的结果
    stop('Only one tree; no need to use getVI, simply draw >RES$vi1 or >RES$vi2, where RES is your BootstrapTreeFitting result')
  }
  
  BN<-dim(RES)[2]  #dim(RES): 6 100 (6为输出的6种$结果: $OOB_predict  $v_predict $vi1   $vi2 $ts1  $ts2 ; 
  #100为boostrapped次数/即RFQT的size/即BN)
  VI<-c()
  if(VItype==2){   #single Q-tree 下的各个candidate variable的VI measurement
    for(i in 1:BN){   VI<-rbind( VI ,  RES[4,i]$vi2 )  }
  }else{
    for(i in 1:BN){   VI<-rbind( VI ,  RES[3,i]$vi1 )  }
  }
  
  VI_means<-apply(VI,2,mean)
  vires<-VI_means*(VI_means >0  ) 
  names(vires)<-paste0( 'Covariate', 1:length(vires)   )
 return( vires )     #这里让负数值全部变成0
}
