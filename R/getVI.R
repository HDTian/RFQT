
#getVI: get the average value of the variable importance (VI) measurement, base on OOB samples
#只返回各个candidate covariate的average VI measurements;不用表示成order顺序了


#get平均的VI measurement
getVI<-function(RES,#RES为parSapply(   cl ,  1:100, BootstrapTreeFitting  )的结果
                VItype=2#VItype=1 for vi1(with label known) VItype=2 for vi2(no label)
                ){
  if(is.null(  dim(RES)  )){ #即， RES只是一次输出结果，不是Sapply的结果
    stop('Only one tree; no need to use getVI, simply draw >RES$vi1 or >RES$vi2, where RES is your BootstrapTreeFitting result')
  }

  BN<-dim(RES)[2]  #dim(RES): 7 100 (7为输出的7种$结果:  $end_node_information $OOB_predict $v_predict $vi1 $vi2 $ts1 $ts2;
  #100为boostrapped次数/即RFQT的size/即BN)
  VI<-c()
  if(VItype==2){   #single Q-tree 下的各个candidate variable的VI measurement
    for(i in 1:BN){   VI<-rbind( VI ,  RES[5,i]$vi2 )  }
  }else{
    for(i in 1:BN){   VI<-rbind( VI ,  RES[4,i]$vi1 )  }
  }

  VI_means<-apply(VI,2,mean)
  vires<-VI_means*(VI_means >0  )#这里让负数值全部变成0
  names(vires)<-paste0( 'Covariate', 1:length(vires)   )
 return( vires )
}



###example:
#RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  )

getVI(RES,VItype=1 )
getVI(RES,VItype=2 )

