

#ֻ���ظ���candidate covariate��VI measurements;���ñ�ʾ��order˳����

getVI<-function(RES,
                VItype=2){  #getƽ����VI measurement
  if(is.null(  dim(RES)  )){ #���� RESֻ��һ��������������Sapply�Ľ��
    stop('Only one tree; no need to use getVI, simply draw >RES$vi1 or >RES$vi2, where RES is your BootstrapTreeFitting result')
  }
  
  BN<-dim(RES)[2]  #dim(RES): 6 100 (6Ϊ�����6��$���: $OOB_predict  $v_predict $vi1   $vi2 $ts1  $ts2 ; 
  #100Ϊboostrapped����/��RFQT��size/��BN)
  VI<-c()
  if(VItype==2){   #single Q-tree �µĸ���candidate variable��VI measurement
    for(i in 1:BN){   VI<-rbind( VI ,  RES[4,i]$vi2 )  }
  }else{
    for(i in 1:BN){   VI<-rbind( VI ,  RES[3,i]$vi1 )  }
  }
  
  VI_means<-apply(VI,2,mean)
  vires<-VI_means*(VI_means >0  ) 
  names(vires)<-paste0( 'Covariate', 1:length(vires)   )
 return( vires )     #�����ø���ֵȫ�����0
}