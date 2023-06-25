
###BootstrapTreeFitting: Q-Tree fit with seed return the relevant results
#BootstrapTreeFitting is particularly suitable for parallel computation

#把sim data和real data整合到一个共有的tree fitting函数好了

#real data的区别：
#(1)存在evidence_effect (now has been removed)
#(2)不存在true effect的概念，也没有MSE，更没有VI1
#(3)real中新加了一个ts1 ts2的概念；sim当然也可以加进去


BootstrapTreeFitting<-function(seed=1,
                               Odat=odat,  #inputed training data
                               Vdat=vdat,  #inputed testing data
                               honest=TRUE, #use honiest estimation or not?
                               S=5,  #maximal depth
                               rate=0.4, # the proportion of candidate variables Ms considered
                               Qthreshold=3.0, ##the threshold for Q heterogneity assessment
                               method='DR',
                               SoP=10, ##size of pre-stratum #only make sense to DR stratification
                               howGX='SpecificGX',##'const' means use extra constant; otherwise estimated by stratum data (stratum-specific GXeffect)
                               endsize=5000,##the two-times of the minimual size of the node of Q-tree. S and endsize work in a similar way. #only makes sense for GetTree
                               const=NA){ #boostrap一次odat得到OOB samples effect prediction;顺便vdat的prediction也做了
  if( is.null(Odat$true_STE[1]) ){JJ<-ncol( Odat )-4}else{JJ<-ncol( Odat )-5}
  if( JJ<1 ){ stop('No candidate covariate, or the dat is not regonized')  }
  #这个JJ主要是给后面的RES$vi用的，
  #至于treefitting中的splitting variable index选择，GetTree自己会再定义(但这两个定义的JJ当然得一样)
  ###Bootstrap一次------------------------------------------------------------------------
  set.seed(seed)
  NNN<-nrow(Odat)
  Bindex<-sample( 1:NNN, NNN , replace=TRUE)
  dat<-Odat[Bindex, ]             #一次simulation中： 分为odat 和 vdat： odat就是training data！
  dat$I<-1:NNN  #当成new data  reindexing ID   #dat就是odat bootstrapped后的for-tree-built的data
  OOBdat<-Odat[  -as.numeric( levels( factor(Bindex  ) ) ), ]  #OOB error和Variable Importance的关键


  ###single tree fitting一次---------------------------------------------------------------
  rdat<-GetTree( dat,  #把dat作为GetTree的输入  #dat可以为bootstrap后的trdat or直接odat
                 S=S,
                 rate=rate,
                 Qthreshold=Qthreshold,
                 method=method,
                 SoP=SoP,
                 howGX=howGX,
                 const=const,
                 endsize=endsize)  #返回rdat这个data.frame  ;  GetTree是最耗时的function

  ###get the MR est for each Nindex based on this single tree fitting----------------------
  NNindex<-as.numeric(  levels(     factor(   rdat$Nindex  )   )   )  #所有NNindex的种类
  MRest<-rep( NA, length( NNindex ) ) #NNindex 和MRest的组合构成了tree的所有结果
  Bx<-c();Bxse<-c();By<-c(); Byse<-c() #为了分析stratum-specific G-X effect
  Sizeratio<-c()
  #此处注意下estdata；如果是odat，那么rdat是自带Nindex；如果是用estdata，那么先算一下estdata$Nindex的Nindex，然后替换rdat为estdat后续代码是一样的
  if(howGX=='const'){
    if(is.na(const)){stop('You used a constant GX association, please tell me the constant value by const=my.const ' ) }
    bx<-const;bxse<-0
    for(nn in 1:length(NNindex)){
      dat_sub<-rdat[abs(rdat$Nindex- NNindex[nn])< .Machine$double.eps^0.5,]#focus on the Nindex-specific subgroup
      fitGY<-lm(    dat_sub[,4]~  dat_sub[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
      Bx<-c(Bx  , bx   ) ;  Bxse<-c(Bxse,  bxse  ) ; By<- c(  By , by ) ; Byse<-c(  Byse, byse )
      Sizeratio<-c( Sizeratio ,nrow( dat_sub) /nrow(rdat)   )  #for ts1  #注意honest style则把rdat换成estdat即可
      MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
      MRest[nn]<-MRres@Estimate
    }
  }else{
    for(nn in 1:length(NNindex)){
      dat_sub<-rdat[abs(rdat$Nindex- NNindex[nn])< .Machine$double.eps^0.5,]
      fitGX<-lm(    dat_sub[,3]~  dat_sub[,2]  )  ; bx<-as.numeric( summary(fitGX)$coef[-1,1]  ); bxse<-as.numeric(  summary(fitGX)$coef[-1,2])
      fitGY<-lm(    dat_sub[,4]~  dat_sub[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
      Bx<-c(Bx  , bx   ) ;  Bxse<-c(Bxse,  bxse  ) ; By<- c(  By , by ) ; Byse<-c(  Byse, byse )
      Sizeratio<-c( Sizeratio ,nrow( dat_sub) /nrow(rdat)   )  #for ts1
      MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
      MRest[nn]<-MRres@Estimate
    }
  }


  ###results calculation---------------------------------------------------------------------
  ###----------------------------------------------------------------------------------------
  RES<-list()

  ###RESULT 1: OOB and testing set predicts
  #OOB individual predicted data
  theMRest<- MRest[match(  as.character(GetNindex(  OOBdat[,5:( JJ+4 ) ] ,rdat,S=S ))  ,  as.character(NNindex)    )  ] #as.character  很重要的操作！ 防止数值误骗
  #注意：NNindex和MRest构成了重要的tree end-node information 即使MRest可能是由estdat估计出来的
  vect<-rep(0,NNN) #NNN 为odat size (无论是否用honest style，这么做都是没问题的)
  vect[  (Odat$I)%in%(OOBdat$I ) ]<-theMRest
  RES$OOB_predict<-vect

  #Validation individual predicted predicted values
  #(如果没有testing set，那么就把v_predict设空)
  if(!exists('vdat')){ RES$v_predict<-NA  }else{
    theMRest<- MRest[match(  as.character(GetNindex(  Vdat[,5:( JJ+4 ) ] ,rdat,S=S ))  ,  as.character(NNindex)    )  ]
    RES$v_predict<-theMRest
  }


  ###RESULT 2: Variable Importance (VI) - for simulation data with label known (vi1) or unpermuted estimates as the true treatment effect (vi2)
  #(VIs 只在training set 和OOB set里进行，和testing set没关系)
  vi1<-'N/A'
  vi2<-c()
  theMRest<-MRest[match(  as.character(GetNindex(  OOBdat[,5:( JJ+4 ) ] ,rdat,S=S ))  ,  as.character(NNindex)    )  ] #这里其实重复计算了一下
  if(  !is.null(OOBdat$true_STE)   ){#即，如果是simulate case
    vi1<-c()
    present_MSE<-mean((OOBdat$true_STE-theMRest)^2)  #当前的single tree MSE
  }
  for(jj in 1:JJ){   #based on OOB set
    ##permuting the jj-th M for OOBdat
    OOBdat_p<-OOBdat[,5:( JJ+4 ) ]  #只保留M信息即可
    OOBdat_p[,jj]<-sample(OOBdat_p[,jj])  #no replacement
    ##get the predicted values and MSE (for sim data)
    theMRest_p<- MRest[match(  as.character(    GetNindex(  OOBdat_p[,1: JJ ] ,rdat,S=S  )  )  ,  as.character(NNindex)    )  ]
    if(  !is.null(OOBdat$true_STE)   ){#即，如果是simulate case
      permuted_MSE<- mean((OOBdat$true_STE-theMRest_p)^2)  #permuted后的MSE ； 一般会更大/越大代表该variable越重要
      vi1<-c(vi1,   permuted_MSE- present_MSE   )
    }
    vi2<-c(vi2,   mean(    (theMRest -theMRest_p)^2   )    )#最开始的unpermuted estimate as the true STE (sharp treatment effect)
  }
  #本次Q-tree内的各个candidate variable的VI measurement
  RES$vi1<-vi1   #如果是用real data,那么vi1就是c() #vi1没法用在real data中
  RES$vi2<-vi2

  ###RESULT 3: Permutation test statistics ts1 and ts2
  #sum(Sizeratio) == 1 #checked!
  ts1<-  sum(Sizeratio*(MRest-  sum(  Sizeratio*MRest  ))^2)
  ##ts2 -> Q statistic:
  Sr<-mr_ivw(mr_input(Bx, Bxse , By, Byse))#先算effect via naive IVW under null (no effect difference)
  delta<-Sr@Estimate #一维的; 已验证！就是y_<-(By/Byse) ;  x_<-1/Byse; lm( y_ ~-1+x_)的结果
  #Cochran's Q statistic
  ssigma_square <- Byse^2 + delta^2*Bxse^2
  ts2<-sum( ( By - delta*Bx     )^2 /ssigma_square      )#critical value:  qchisq(0.95, NC-1)==3.841459
  RES$ts1<-ts1
  RES$ts2<-ts2

  return( RES )
}

#RES list:  $OOB_predict $v_predict $vi1 $vi2 $ts1 $ts2



###examples:

set.seed(60)
res<-getDat() #simulated data
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set

#When running RFQT with mutiple Q trees/bootstrap - use parallel computation
n_cores_used<-detectCores()-1
Nb<-100 #how many trees in the forest?
cl<-makeCluster(n_cores_used)#规定用多少并行clusters/线程
clusterEvalQ(cl=cl , expr=library(dplyr))  #给各个cluster中的运行一些表达式expression (例如导入一些包；做一些基础操作)
clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat',
                                    'GetTree', 'GetNindex', 'GetIndex' )  )#or any other arguments
RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  ) #本脚本里parSapply的结果常用RES表示
stopCluster(cl)

##If you wish to use your own parameters rather than the default parameters, try:
#general exmaple
user_BootstrapTreeFitting<-function(seed){
  RES<-BootstrapTreeFitting(seed,
                            S=my.S,
                            JJ=my.JJ,   #or any partial of the arguments
                            rate=my.rate,
                            Qthreshold=my.Qthreshold,
                            method=my.method,
                            SoP=my.SoP,
                            howGX=my.howGX,
                            endsize=my.endsize)
  return(RES)
}
#specific exmaple
user_BootstrapTreeFitting<-function(seed){
  RES<-BootstrapTreeFitting(seed,SoP=20)
  return(RES)
}

###parallel computation  -> RFQT
n_cores_used<-detectCores()-1
Nb<-100 #how many trees in the forest?
cl<-makeCluster(n_cores_used)#规定用多少并行clusters/线程
clusterEvalQ(cl=cl , expr=library(dplyr))  #给各个cluster中的运行一些表达式expression (例如导入一些包；做一些基础操作)
clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat',
                                    'GetTree', 'GetNindex', 'GetIndex' , 'BootstrapTreeFitting')  )
RES<-parSapply(   cl ,  1:Nb, user_BootstrapTreeFitting  ) #本脚本里parSapply的结果常用RES表示
stopCluster(cl)

