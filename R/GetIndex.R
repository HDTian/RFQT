
#GetIndex:accodring to the current data, return the best candidate M (i.e. the most possible modifiers) and its splitting point (1:1? 1:2? 2:1?)


GetIndex<-function(dat_current,#dat_current: current data #must be data.frame with first four columns: $I $Z $X $Y
                   JJ,#the number of total candidate variables
                   rate=1,# rate: the proportion of M considered #if SpecificM is used, rate will not make a difference
                   SpecificM=NA, #vector: user specific M index
                   method='DR',#stratification method used : 'DR' 'Residual' others
                   SoP=10,#size of pre-stratum  #SoP=10: better for operation: (1,1,1,2,2,3,3,4,4,4)
                   howGX='SpecificGX', #how to caculate the GX effect?  'const' means use extra constant; otherwise estimated by stratum data (stratum-specific GXeffect)
                   Halve=FALSE, #if only use half splitting? Default is FALSE, i.e. use three splitting possiblies: (3:7) (5:5) (7:3)
                   const=NA){
  N<-dim( dat_current )[1]
  ###getIndex这个函数不适合通过dat_current来计算JJ，因为在tree growth中，dat_current会不断增广
  #JJ是GetIndex这个函数外已经被定义的量
  if((rate>1)|(rate<=0)){ stop('PLease use a correct rate (0<= rate <1)')}
  J<-round(rate*JJ)
  if(J<1){  stop('No candidiate cvoariate to be determine, use a larger rate value')     }
  if(is.na(SpecificM)){
    Mcandidates<-sample(1:JJ,J)
    }else{
      if(  sum(!SpecificM%in%(1:JJ))>0  ){stop('SpecificM is not correctly indexed!')}
      Mcandidates<-SpecificM
  }
  #print(Mcandidates)
  QQs_allM<-c()
  Size_allM<-c()#record the mininal node size after the specific split for specific M
  #Ns<-2 #在RFQT中，一直只用2个strata
  for(j in Mcandidates){#J: the number of variables Ms considered (注意，每一次search都要用随机的partial M candidates，为了de-correlation trees in random forest)
    ###stratification
    if(method=='DR'){
      ###DR---
      dat_current$Mobj<-dat_current[,4+j] #J-th M的值； 注意前4项是I Z X Y
      dat_order<-dat_current[ order(dat_current[,2]  ),  ]  #ordered by Z

      dat_order$pre_stratum<-    rep(1:(floor(N/SoP)+1), each=SoP,length.out=N)#(floor(N/SoP)+1)*SoP >= N 保证能超过就行

      #rank twice (ie doubly-ranked)
      temp<-arrange(  dat_order, Mobj )  #按照Mobj升序排一下  #arrange()应该也没有随机性
      dat_order<-arrange(  temp ,pre_stratum ) #即，保证pre_strata按顺序排列，并且每个pre_strata中的目标量都是升序
      dat_order$strata<-as.vector( unlist(sapply( as.numeric(table( dat_order$pre_stratum )) , function(x) sort(rep(c(1,1,1,2,2,3,3,4,4,4),length.out=x)) )   ) )

    }else{
      if(method=='Residual'){
        ###Residual---
        dat_current$Mobj<-dat_current[,4+j] #J-th M的值； 注意前4项是I Z X Y
        dat_order<-dat_current#[ order(dat_current[,2]  ),  ]  #ordered by Z
        dat_order$residual<-dat_order$Mobj-lm(  dat_order$Mobj~dat_order$Z   )$fitted

        dat_order<-dat_order[  order(dat_order$residual  )  ,  ] #ordered by residuals
        dat_order$strata<- sort(    rep(c(1,1,1,2,2,3,3,4,4,4), length.out=N)  )#一定可以控制Ns

      }else{
        ###Naive rank---
        dat_current$Mobj<-dat_current[,4+j] #J-th M的值； 注意前4项是I Z X Y
        dat_order<-dat_current#[ order(dat_current[,2]  ),  ]  #ordered by Z

        dat_order<-dat_order[  order(dat_order$Mobj  )  ,  ] #ordered by Mobj
        dat_order$strata<- sort(    rep(c(1,1,1,2,2,3,3,4,4,4), length.out=N)  )#一定可以控制Ns

      }
    }#end of stratification

    ###IV estimates---
    #RES<-c();Means<-c()
    ##3种dochotomising scenarios: c(1,1,1,2,2,3,3,4,4,4)
    #split style 1: 111   2222222
    #split style 2: 11111   22222
    #split style 3: 1111111   222
    QQs<-c() #final length of QQs: 3
    Size<-c() #store the minial node size for each split
    
    if(Halve){ dd_candidates<-2    }else{ dd_candidates<-1:3  }
    for(dd in dd_candidates){ #原本是 dd in 1:3； 如果只用Halve 则dd in 2
      strata_range<-list()
      strata_range[[1]]<-(1:4)[1:dd]
      strata_range[[2]]<-(1:4)[-(1:dd)]
      Bx1<-c(); Bxse1<-c(); By1<-c() ; Byse1<-c()  #这些来算Q足矣
      ssize<-c()
      for(i in 1:2){
        dat_sub<-dat_order[(dat_order$strata)%in%strata_range[[i]],]
        ssize<-c(ssize,nrow(dat_sub))
        #Means<-rbind(Means,summary(dat_sub$Mobj  ) )
        if(howGX=='const'){
          if(is.na(const)){stop('You use a constant GX association, please tell me the constant value by const=my.const ' ) }
          bx<-const#population估计出来的G-X effect
          bxse<-0
        }else{
          fitGX<-lm(    dat_sub[,3]~  dat_sub[,2]  )  ; bx<-as.numeric( summary(fitGX)$coef[-1,1]  ); bxse<-as.numeric(  summary(fitGX)$coef[-1,2])
        }
        fitGY<-lm(    dat_sub[,4]~  dat_sub[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
        Bx1<-c( Bx1, bx  ); Bxse1<-c( Bxse1, bxse  ); By1<-c(By1, by  ) ; Byse1<-c(  Byse1,  byse)
        #MRres<-mr_ivw(mr_input(bx, bxse, by, byse)); res1<-c(MRres@Estimate , MRres@StdError)
        #RES<-rbind( RES,   c( res1 )   )
      }
      #RES;Means

      ###Q statistic---
      Sr<-mr_ivw(mr_input(Bx1, (1:2) , By1, Byse1))#先算effect via naive IVW under null (no effect difference)
      delta<-Sr@Estimate #一维的
      #Cochran's Q statistic
      ssigma_square <- Byse1^2 + delta^2*Bxse1^2
      QQ<-sum( ( By1 - delta*Bx1     )^2 /ssigma_square      )#critical value:  qchisq(0.95, NC-1)==3.841459

      ###Storage
      QQs<-c(QQs,QQ )
      Size<-c(Size,min(ssize))
    }#end of the 3 dichotomization for the current M

    QQs_allM<-rbind(QQs_allM,QQs) # #dim(QQs_allM) #JJ 3
    Size_allM<-rbind(Size_allM,Size) ##dim(Size_allM) #JJ 3
  } #end of the Q calculation for all Ms     #dim(QQs_allM) #JJ 3

  ##result storage
  row_col<-which(QQs_allM == max(QQs_allM), arr.ind = TRUE) #vector of 2 elements: the row and col index
  results<-c(  Mcandidates[row_col[1]],  max(QQs_allM), Size_allM[row_col[1],row_col[2]] , row_col[2]  )
  
  if(Halve){ results[4]<-results[4]+1   }#为了让split.style的数字为2，此时对应着5:5
  
  names(results  )<-c(  'Candidate.index', 'Q.value' , 'minimal.node.size.after.split','split.style' )
  return( results) #返回M编号 和 最大的Q值 和此时date (没分之前的)sample size [这两个用在stoping rule]
}
#返回 1.M编号 和 2.对应的最大的Q值 和 3.此时date分之后的最小的node size [这两个用在stoping rule] 4.split.style


#splitting style:
#1为 3：7   2为5：5    3为7：3

###exmaples:
#set.seed(60)
#res<-getDat() #simulated data #the deflaut setting: scenario='A' and SoM=0.5
#odat<-res$traning.set  #training set
#vdat<-res$testing.set  #testing set
#GetIndex(odat,JJ=20)
#Candidate.index         Q.value       node.size
#         3.0000        324.7219     100000.0000


#Candidate.index                      Q.value     minimal.node.size.after.split    split.style
#         3.0000                      382.1728    30000.0000                       1.0000


#split.style: (lower node : upper node)
#1: 3:7
#2: 5:5
#3: 7:3



