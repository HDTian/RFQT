###
###GetTree    #construct a single Q-tree (e.g. for training dat or treedat) and store the decision rule
###return: rdat (the rdat contains the decision rule

GetTree<-function(dat,#input: the data set (either training dats or tree data);  must be data.frame ;the I index must be strictly increasing
                  S=5, #the largest depth #only makes sense for GetTree
                  Qthreshold=3.84,#the threshold for Q heterogneity assessment #only makes sense for GetTree
                  rate=1,# the proportion of candidate variables Ms considered
                  method='DR',
                  SoP=10, #size of pre-stratum #only make sense to DR stratification
                  howGX='SpecificGX',#'const' means use extra constant; otherwise estimated by stratum data (stratum-specific GXeffect)
                  const=NA, #only make sense when howGX='const'
                  endsize=5000 #the minimual size of the node of Q-tree allowed to exist. S and endsize work in a similar way. #only makes sense for GetTree
                  ){
  #先定义一个JJ；这个JJ是给Getindex用的
  if( is.null(dat$true_STE[1]) ){JJ<-ncol( dat )-4}else{JJ<-ncol( dat )-5}#如果是sim data，就减去5，否则减去4

  #check the IDindex in dat is strictly increasing (otherwise induce problems matching the stratification results)
  if(  sum((dat$I[-1]-dat$I[-length(dat$I)])<=0)>0  ){
    stop('the ID index of the inputed data is not strictly increasing: please first arrange them by yourself')
  }

  dat$Nindex<-0 #最开始先提供统一的Nindex: 0 #Nindex: 目前的Node index (对应着所有的current end node)
  NDR<-SoP
  NNN<-nrow(dat)
  #Ns<-2 #在RFQT中，一直只用2个strata
  for(s in 1:S){  #S:the largest depth #s: the S round
    infvec<-rep(NA,NNN)  # information vector; 最后加到dat上，作为reference table
    #很棒，这解决了我的一个担忧：每次S生长时，会重新全部设空一次infvec
    #但还有一个问题：已经停止的node是否还需要再生长？
    Ls<-levels(factor(    dat$Nindex )    ) #一定要先存好，否则loop中会不断更新dat
    for(k in 1:length(Ls)  ){ #Nindex: Node index
      #go to the specific Nindex subgroup/node
      dat_current<-dat[ abs( dat$Nindex-as.numeric(  Ls  )[k])<.Machine$double.eps^0.5 ,  ]
      #get the best candidate M for this node
      if(nrow(dat_current)>=endsize*2){#如果当前data小于endsize*2 那么无论怎么分都会出现小于endsize的node;没必要花时间算GetIndex
        GIres<-GetIndex(dat_current,
                        JJ=JJ,
                        rate=rate,
                        method=method,
                        SoP=SoP,
                        howGX=howGX,
                        const=const)  #GetIndex function result       #这行代码占了一半！
        J_<-GIres[1] #每次用之前记得指定处理这个dat中M的数量;GetIndex已经自动随机patrial选择M candidates了
        #这里主要是告诉我们用哪个M的
        split_style_used<-GIres[4]
        #以及用什么split_style

        if( (GIres[2]>Qthreshold)&(  GIres[3]>endsize ) ){  #只有同时满足 此时才会接着细分，否则不用管
          N<-dim( dat_current )[1] #其实奇数也不要紧；rank无视单数奇数
          #根据split_style_used设置the vector for rep
          repvector_<-!(c(1,1,1,2,2,3,3,4,4,4)%in%(1:split_style_used) )
          repvector<-repvector_+1
          #split style 1: 111   2222222
          #split style 2: 11111   22222
          #split style 3: 1111111   222

          #stratification; i.e. node split
          if(method=='DR'){
            ###DR---
            dat_current$Mobj<-dat_current[,4+J_] #J-th M的值； 注意前4项是I Z X Y
            dat_order<-dat_current[ order(dat_current[,2]  ),  ]  #ordered by Z

            dat_order$pre_stratum<-    rep(1:(floor(N/SoP)+1), each=SoP,length.out=N)#(floor(N/SoP)+1)*SoP >= N 保证能超过就行
            #rank twice (ie doubly-ranked)
            temp<-arrange(  dat_order, Mobj )  #按照Mobj升序排一下  #arrange()应该也没有随机性
            dat_order<-arrange(  temp ,pre_stratum ) #即，保证pre_strata按顺序排列，并且每个pre_strata中的目标量都是升序
            dat_order$strata<-as.vector( unlist(sapply( as.numeric(table( dat_order$pre_stratum )) , function(x) sort(rep(repvector,length.out=x)) )   ) )

          }else{
            if(method=='Residual'){
              ###Residual---
              dat_current$Mobj<-dat_current[,4+J_] #J-th M的值； 注意前4项是I Z X Y
              dat_order<-dat_current#[ order(dat_current[,2]  ),  ]  #ordered by Z
              dat_order$residual<-dat_order$Mobj-lm(  dat_order$Mobj~dat_order$Z   )$fitted

              dat_order<-dat_order[  order(dat_order$residual  )  ,  ] #ordered by residuals
              dat_order$strata<- sort(    rep(repvector, length.out=N)  )#一定可以控制Ns

            }else{
              ###Naive rank---
              dat_current$Mobj<-dat_current[,4+J_] #J-th M的值； 注意前4项是I Z X Y
              dat_order<-dat_current#[ order(dat_current[,2]  ),  ]  #ordered by Z

              dat_order<-dat_order[  order(dat_order$Mobj  )  ,  ] #ordered by Mobj
              dat_order$strata<- sort(    rep(repvector, length.out=N)  )#一定可以控制Ns

            }
          }

          ##split/tree decision information storage: Mindex + left-lower/right-upper go + mean M values + split style (后两者决定cuttingvalue)
          INFvect<-rep(NA,N)
          INFvect[dat_order$strata==1  ]<-paste0( J_ , '_' , 1, '_', mean( dat_order$Mobj[dat_order$strata==1] ), '_' ,split_style_used  )#1: lower M, left node
          INFvect[dat_order$strata==2  ]<-paste0( J_ , '_' , 2, '_', mean( dat_order$Mobj[dat_order$strata==2] ), '_' ,split_style_used  )#2: upper M, right node
          dat_order$inf<-INFvect
          #cuttingvalue<-  (mean( dat_order$Mobj[dat_order$strata==1] ) +  mean( dat_order$Mobj[dat_order$strata==2] ) )/2
          #this cuttingvalue means new input will go to the node that has more closer mean M value

          ###storage/merge into dat---
          #update the Nindex
          dat_current<-arrange(  dat_order, I )  #最后按照I排序，方便融入dat中; #dat_current$I 中肯定不包含Iindex
          vect<-rep(0,NNN)  #NNN为dat (total) sample size
          vect[(dat$I)%in%(dat_current$I )]<-dat_current$strata*0.1^s #times 0.1^s 相当于在Nindex后面添上新的tree information
          #update the Nindex
          dat$Nindex<-dat$Nindex+vect#dat是total dat
          #update infvec (infvec可能是有些才会update，否则就是NA)
          infvec[  (dat$I)%in%(dat_current$I )  ]<- dat_current$inf  #其实和Nindex应该要匹配
        }
      } #end of single node split in the current end node situation
    } #end of the node split for all current Nindex in the present S round

    ###added the infvec representing the current s round
    dat<-cbind( dat , infvec )  #不用管NA项目，反正都是 给定Nindex具体值再做处理
    names(dat)[dim(dat)[2]    ]<-paste0( 'infvec_',s  ) #给最后一项换个名字，防止重名导致无法使用arrange
  } #end of all S round; the tree stop growing

  ###final Reference table
  rdat<-dat #rdat和dat一致 #rdat是dat的增广矩阵(包含每轮s生长的tree information / decision rule)  #不需要用subdata: dat[, ( (dim(dat)[2]-S+1)  :   (dim(dat)[2]) )]

  return(rdat)
}


#返回一个和inputed dat增广的dat set(包含每轮s生长的tree information)

###exmaples:
set.seed(60)
res<-getDat() #simulated data #the deflaut setting: scenario='A' and SoM=0.5
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set
rdat<-GetTree(odat)

#easy check
View(rdat)
table(rdat$infvec_3)

