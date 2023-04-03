###
###GetTree    #get a single tree fitting
###

GetTree<-function(dat,#input: the training set; must be data.frame ; return: rdat
                  S=5, #only makes sense for GetTree
                  Qthreshold=3.0,#only makes sense for GetTree
                  rate=1,
                  method='DR',
                  SoP=10,
                  howGX='SpecificGX',
                  const=NA,
                  endsize=5000 #only makes sense for GetTree
                  ){
  #先定义一个JJ；这个JJ是给Getindex用的
  if( is.null(dat$true_STE[1]) ){JJ<-ncol( dat )-4}else{JJ<-ncol( dat )-5}
  dat$Nindex<-0 #最开始先提供统一的Nindex
  NDR<-SoP
  NNN<-nrow(dat)
  Ns<-2 #在RFQT中，一直只用2个strata
  for(s in 1:S){  #S: total  split times
    infvec<-rep(NA,NNN)  # information vector; 最后加到dat上，作为reference table
    #很棒，这解决了我的一个担忧：是每次S生长时，重新全部设空一次infvec
    #但还有一个问题：已经停止的node是否还需要再生长？
    Ls<-levels(factor(    dat$Nindex )    ) #一定要先存好，否则loop中会不断更新dat
    for(k in 1:length(Ls)  ){ #Nindex: Node index
      dat_current<-dat[ abs( dat$Nindex-as.numeric(  Ls  )[k])<.Machine$double.eps^0.5 ,  ]
      GIres<-GetIndex(dat_current,
                      JJ=JJ,
                      rate=rate,
                      method=method,
                      SoP=SoP,
                      howGX=howGX,
                      const=const)  #GetIndex function result       #这行代码占了一半！
      J_<-GIres[1] #每次用之前记得指定处理这个dat中M的数量;GetIndex已经自动随机patrial选择M candidates了
      #这里主要是告诉我们用哪个M的
      if( (GIres[2]>Qthreshold)&(  GIres[3]>endsize ) ){  #只有此时才会接着细分，否则不用管  (注意 final end node size >= endsize/2)
        N<-dim( dat_current )[1] #其实奇数也不要紧；rank无视单数奇数
        if(method=='DR'){
          ###DR---
          dat_current$Mobj<-dat_current[,4+J_] #J-th M的值； 注意前4项是I Z X Y
          dat_order<-dat_current[ order(dat_current[,2]  ),  ]  #ordered by Z
          
          dat_order$pre_stratum<-    rep(1:(floor(N/SoP)+1), each=SoP,length.out=N)     
          #(floor(N/SoP)+1)*SoP >= N 保证能超过就行
          #rank twice (ie doubly-ranked)
          temp<-arrange(  dat_order, Mobj )  #按照Mobj升序排一下  #arrange()应该也没有随机性
          dat_order<-arrange(  temp ,pre_stratum ) #即，保证pre_strata按顺序排列，并且每个pre_strata中的目标量都是升序
          dat_order$strata<-as.vector( unlist(sapply( as.numeric(table( dat_order$pre_stratum )) , function(x) sort(rep(1:Ns,length.out=x)) )   ) )
          
          # No<-round(N/NDR)  #NDR是我们想控制的！ 含义和N/No类似：代表着每个pre-strata中individual的数量
          # dat_order$pre_strata<- sort(    rep(1:No, length.out=N)  )      #No: # of pre-strata ; 不用纠结N/No是否为整数或者双数
          # temp<-arrange(  dat_order, Mobj )  #此时Mobj为目标量！
          # dat_order<-arrange(  temp ,pre_strata ) #即，保证pre_strata按顺序排列，并且每个pre_strata中的目标量都是升序
          # dat_order$strata<-as.vector( unlist(    sapply( as.numeric(table( dat_order$pre_strata ))    ,    function(x) sort(rep(  1:2, length.out=x ))    )   ) )
          
          INFvect<-rep(NA,N)
          INFvect[dat_order$strata==1  ]<-paste0( J_ , '_' , 1, '_', mean( dat_order$Mobj[dat_order$strata==1] ) )
          INFvect[dat_order$strata==2  ]<-paste0( J_ , '_' , 2, '_', mean( dat_order$Mobj[dat_order$strata==2] ) )
          dat_order$inf<-INFvect
          cuttingvalue<-  (mean( dat_order$Mobj[dat_order$strata==1] ) +  mean( dat_order$Mobj[dat_order$strata==2] ) )/2
        }else{
          if(method=='Residual'){
            ###Residual---
            dat_current$Mobj<-dat_current[,4+J_] #J-th M的值； 注意前4项是I Z X Y
            dat_order<-dat_current#[ order(dat_current[,2]  ),  ]  #ordered by Z
            dat_order$residual<-dat_order$Mobj-lm(  dat_order$Mobj~dat_order$Z   )$fitted
            
            dat_order<-dat_order[  order(dat_order$residual  )  ,  ] #ordered by residuals
            dat_order$strata<- sort(    rep(1:Ns, length.out=N)  )#一定可以控制Ns
            
            #dat_order$strata<-  floor( (rank( dat_order$residual,ties.method ='random' )/((N/2)+0.000000001) ) )+1 #N为奇数也不要紧啊
            
            INFvect<-rep(NA,N)
            INFvect[dat_order$strata==1  ]<-paste0( J_ , '_' , 1, '_', mean( dat_order$Mobj[dat_order$strata==1] ) )
            INFvect[dat_order$strata==2  ]<-paste0( J_ , '_' , 2, '_', mean( dat_order$Mobj[dat_order$strata==2] ) )
            dat_order$inf<-INFvect
            cuttingvalue<-  (mean( dat_order$Mobj[dat_order$strata==1] ) +  mean( dat_order$Mobj[dat_order$strata==2] ) )/2
          }else{
            ###Naive rank---
            dat_current$Mobj<-dat_current[,4+J_] #J-th M的值； 注意前4项是I Z X Y
            dat_order<-dat_current#[ order(dat_current[,2]  ),  ]  #ordered by Z
            
            dat_order<-dat_order[  order(dat_order$Mobj  )  ,  ] #ordered by Mobj
            dat_order$strata<- sort(    rep(1:Ns, length.out=N)  )#一定可以控制Ns
            
            #dat_order$strata<-  floor( (rank( dat_order$Mobj,ties.method ='random' )/((N/2)+0.000000001) ) )+1
            
            INFvect<-rep(NA,N)
            INFvect[dat_order$strata==1  ]<-paste0( J_ , '_' , 1, '_', mean( dat_order$Mobj[dat_order$strata==1] ) )
            INFvect[dat_order$strata==2  ]<-paste0( J_ , '_' , 2, '_', mean( dat_order$Mobj[dat_order$strata==2] ) )
            dat_order$inf<-INFvect
            cuttingvalue<-  (mean( dat_order$Mobj[dat_order$strata==1] ) +  mean( dat_order$Mobj[dat_order$strata==2] ) )/2
          }
        }

        ###storage into dat---
        #update the Nindex
        dat_current<-arrange(  dat_order, I )  #最后按照I排序，方便融入dat中; #dat_current$I 中肯定不包含Iindex
        vect<-rep(0,NNN)  #NNN为dat (total) sample size
        vect[(dat$I)%in%(dat_current$I )]<-dat_current$strata*0.1^s
        #update the Nindex
        dat$Nindex<-dat$Nindex+vect
        #update infvec
        infvec[  (dat$I)%in%(dat_current$I )  ]<- dat_current$inf  #其实和Nindex应该要匹配
      }
    }
    dat<-cbind( dat , infvec )  #不用管NA项目，反正都是 给定Nindex具体值再做处理
    names(dat)[dim(dat)[2]    ]<-paste0( 'infvec_',s  ) #最后一项换个名字，防止重名导致无法使用arrange
  }
  #Reference table
  rdat<-dat #rdat和dat一致 #rdat是dat的增广矩阵  #不需要用subdata: dat[, ( (dim(dat)[2]-S+1)  :   (dim(dat)[2]) )]

  return(rdat)
}


#返回一个和inputed dat增广的dat set
