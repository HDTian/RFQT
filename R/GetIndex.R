
#dat_current : usually the training data

#NDR就是SoP

GetIndex<-function(dat_current,
                   JJ,
                   rate=1,
                   method='DR',
                   SoP=10,
                   howGX='SpecificGX',
                   const=NA){#dat_current: current data; rate: the proportion of M considered
  N<-dim( dat_current )[1]
  #getIndex这个函数不适合通过dat_current来计算JJ，因为在tree growth中，dat_current会不断增广
  #JJ是GetIndex这个函数外已经被定义的量
  if((rate>1)|(rate<=0)){ stop('PLease use a correct rate (0<= rate <1)')}
  J<-round(rate*JJ)
  if(J<1){  stop('No candidiate cvoariate to be determine, use a larger rate value')     }
  Mcandidates<-sample(1:JJ,J)
  QQs<-c()
  Ns<-2 #在RFQT中，一直只用2个strata
  for(j in Mcandidates){#J: the number of total variables Ms (注意，每一次search都要用随机的partial M candidates，为了de-correlation trees in random forest)
    if(method=='DR'){
      ###DR---
      dat_current$Mobj<-dat_current[,4+j] #J-th M的值； 注意前4项是I Z X Y
      dat_order<-dat_current[ order(dat_current[,2]  ),  ]  #ordered by Z
      
      dat_order$pre_stratum<-    rep(1:(floor(N/SoP)+1), each=SoP,length.out=N)     
      #(floor(N/SoP)+1)*SoP >= N 保证能超过就行
      #rank twice (ie doubly-ranked)
      temp<-arrange(  dat_order, Mobj )  #按照Mobj升序排一下  #arrange()应该也没有随机性
      dat_order<-arrange(  temp ,pre_stratum ) #即，保证pre_strata按顺序排列，并且每个pre_strata中的目标量都是升序
      dat_order$strata<-as.vector( unlist(sapply( as.numeric(table( dat_order$pre_stratum )) , function(x) sort(rep(1:Ns,length.out=x)) )   ) )
      
      # No<-round(N/SoP)  #NDR是我们想控制的！ 含义和N/No类似：代表着每个pre-strata中individual的数量
      # dat_order$pre_strata<- sort(    rep(1:No, length.out=N)  )      #No: # of pre-strata ; 不用纠结N/No是否为整数或者双数
      # temp<-arrange(  dat_order, Mobj )  #此时Mobj为目标量！
      # dat_order<-arrange(  temp ,pre_strata ) #即，保证pre_strata按顺序排列，并且每个pre_strata中的目标量都是升序
      # dat_order$strata<-as.vector( unlist(    sapply( as.numeric(table( dat_order$pre_strata ))    ,    function(x) sort(rep(  1:2, length.out=x ))    )   ) )
    }else{
      if(method=='Residual'){
        ###Residual---
        dat_current$Mobj<-dat_current[,4+j] #J-th M的值； 注意前4项是I Z X Y
        dat_order<-dat_current#[ order(dat_current[,2]  ),  ]  #ordered by Z
        dat_order$residual<-dat_order$Mobj-lm(  dat_order$Mobj~dat_order$Z   )$fitted
        
        dat_order<-dat_order[  order(dat_order$residual  )  ,  ] #ordered by residuals
        dat_order$strata<- sort(    rep(1:Ns, length.out=N)  )#一定可以控制Ns
        
        #dat_order$strata<-  floor( (rank( dat_order$residual,ties.method ='random' )/((N/2)+0.000000001) ) )+1 #N为奇数也不要紧啊
      }else{
        ###Naive rank---
        dat_current$Mobj<-dat_current[,4+j] #J-th M的值； 注意前4项是I Z X Y
        dat_order<-dat_current#[ order(dat_current[,2]  ),  ]  #ordered by Z
        
        dat_order<-dat_order[  order(dat_order$Mobj  )  ,  ] #ordered by Mobj
        dat_order$strata<- sort(    rep(1:Ns, length.out=N)  )#一定可以控制Ns
        
        #dat_order$strata<-  floor( (rank( dat_order$Mobj,ties.method ='random' )/((N/2)+0.000000001) ) )+1
      }
    }

    ###IV estimates---
    RES<-c();Means<-c()
    Bx1<-c(); Bxse1<-c(); By1<-c() ; Byse1<-c()  #这三个是用来算Q的
    for(i in 1:Ns){
      dat_sub<-dat_order[dat_order$strata==i,]
      Means<-rbind(Means,summary(dat_sub$Mobj  ) )
      if(howGX=='const'){
        if(is.na(const)){stop('You use a constant GX association, please tell me the constant value by const=my.const ' ) }
        bx<-const#population估计出来的G-X effect
        bxse<-0
      }else{
        fitGX<-lm(    dat_sub[,3]~  dat_sub[,2]  )  ; bx<-as.numeric( summary(fitGX)$coef[-1,1]  ); bxse<-as.numeric(  summary(fitGX)$coef[-1,2])
      }
      fitGY<-lm(    dat_sub[,4]~  dat_sub[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
      Bx1<-c( Bx1, bx  ); Bxse1<-c( Bxse1, bxse  ); By1<-c(By1, by  ) ; Byse1<-c(  Byse1,  byse)
      MRres<-mr_ivw(mr_input(bx, bxse, by, byse)); res1<-c(MRres@Estimate , MRres@StdError)
      RES<-rbind( RES,   c( res1 )   )
    }
    RES;Means

    ###Q statistic---
    Sr<-mr_ivw(mr_input(Bx1, (1:2) , By1, Byse1))#先算effect via naive IVW under null (no effect difference)
    delta<-Sr@Estimate #一维的
    #Cochran's Q statistic
    ssigma_square <- Byse1^2 + delta^2*Bxse1^2
    QQ<-sum( ( By1 - delta*Bx1     )^2 /ssigma_square      )#critical value:  qchisq(0.95, NC-1)==3.841459

    ###Storage
    QQs<-c(QQs,QQ )
  }
  results<-c(  Mcandidates[which.max(QQs )],  max(QQs), N  )
  names(results  )<-c(  'Candidate.index', 'Q.value' , 'node.size' )
  return( results) #返回M编号 和 最大的Q值 和此时date (没分之前的)sample size [这两个用在stoping rule]
}
#返回M编号 和 最大的Q值 和此时date (没分之前的)sample size [这两个用在stoping rule]

