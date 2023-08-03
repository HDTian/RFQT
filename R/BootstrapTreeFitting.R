
###BootstrapTreeFitting: Q-Tree fit with seed return the relevant results
#BootstrapTreeFitting is particularly suitable for parallel computation

#把sim data和real data整合到一个共有的tree fitting函数好了

#real data的区别：
#(1)存在evidence_effect (now has been removed)
#(2)不存在true effect的概念，也没有MSE，更没有VI1
#(3)real中新加了一个ts1 ts2的概念；sim当然也可以加进去

#注意：和DTfit相比，区别在于： 
#没有bootstrap, 因此没有OOB，也没有VI1 VI2；只有$end_node_information $v_predict $ts1 $ts2
#强制默认rate=1


BootstrapTreeFitting<-function(seed=1,
                               Odat=odat,  #inputted training data
                               Vdat=vdat,  #inputted testing data  #Vdat can be empty or nonexistent, not defined
                               honest=FALSE, #use honest estimation or not?
                               S=5,  #maximal depth
                               rate=0.4, # the proportion of candidate variables Ms considered #if SingleM=TRUE, rate will be ignored
                               SingleM=FALSE, #whether to ue single stratification style?(i,e, use a fixed M determined in the begining)
                               Qthreshold=3.84, ##the threshold for Q heterogneity assessment
                               method='DR',#stratification method used: 'DR' 'Residual' others
                               SoP=10, ##size of pre-stratum #only make sense to DR stratification
                               howGX='SpecificGX',##'const' means use extra constant; otherwise estimated by stratum data (stratum-specific GXeffect)
                               Halve=FALSE, #if only use half splitting? Default is FALSE, i.e. use three splitting possiblies: (3:7) (5:5) (7:3)
                               endsize=1000,##the minimal size of the node of Q-tree allowed to exist
                               const=NA){ #boostrap一次odat得到OOB samples effect prediction;顺便vdat的prediction也做了
  print(paste0( 'seed:',seed )  )
  if( is.null(Odat$true_STE[1]) ){JJ<-ncol( Odat )-4}else{JJ<-ncol( Odat )-5}
  if( JJ<1 ){ stop('No candidate covariate, or the dat is not regonized')  }
  #这个JJ主要是给后面的RES$vi用的，
  #至于treefitting中的splitting variable index选择，GetTree自己会再定义(但这两个定义的JJ当然得一样)
  
  
  #check the IDindex in Odat is strictly increasing
  #(因为最后返回的结果OOBpredict是和odat index顺序匹配的)
  if(  sum((Odat$I[-1]-Odat$I[-length(Odat$I)])<=0)>0  ){
    stop('the ID index of the inputed Odat is not strictly increasing: please first arrange them by yourself')
  }
  ###Bootstrap一次------------------------------------------------------------------------
  set.seed(seed)
  NNN<-nrow(Odat)
  Bindex<-sample( 1:NNN, NNN , replace=TRUE)
  bdat<-Odat[Bindex, ]  #bootstrapped data         #一次simulation中： 分为odat 和 vdat： odat就是training data！
  bdat$I<-1:NNN  #当成new data and !reindexing! ID  #ID严格升序
  #bdat就是odat bootstrapped后的for-tree-built的data (即可当作training dat or divide into treedat + estdat)
  OOBdat<-Odat[  -as.numeric( levels( factor(Bindex  ) ) ), ]  #OOBdat的I是严格升序的
  #OOB error和Variable Importance的关键
  dat<-bdat
  
  if(honest){ #honest estimation
    #random split (注意：必须要random split；但是由于bootstrap index是随机的且没有sort,可以直接按照bdat的I顺序来分)
    treedat<-bdat[(1:round(nrow(bdat)/2)),]#bdat的I是严格升序的
    estdat<-bdat[-(1:round(nrow(bdat)/2)),]
    dat<-treedat
  }
  
  
  ###single tree fitting一次---------------------------------------------------------------
  #如果使用single stratification；需要先依据dat判断出one single best candidate M,然后使用SpecificM
  if(SingleM){
    GetIndex_res<-GetIndex(dat,
                           JJ=JJ, #JJ has been decided in the beginning of the BootstrapTreeFitting function
                           rate=1, #must be  1
                           SpecificM=NA, #must be NA #vector: user specific M index
                           method=method,#stratification method used
                           SoP=SoP,#size of pre-stratum  #SoP=10: better for operation: (1,1,1,2,2,3,3,4,4,4)
                           howGX=howGX, #how to calculate the GX effect?  'const' means use extra constant; otherwise estimated by stratum data (stratum-specific GXeffect)
                           Halve=Halve, #if only use half splitting? Default is FALSE, i.e. use three splitting possiblies: (3:7) (5:5) (7:3)
                           const=const)
    SpecificM_used<-GetIndex_res[1]
  }else{
    SpecificM_used<-NA
  }
  
  rdat<-GetTree( dat,  #把dat作为GetTree的输入  #dat可以为bdat或者的treedat或者直接odat
                 S=S,
                 rate=rate,
                 SpecificM=SpecificM_used,#一旦SpecificM被启用，那么rate会被忽视，不用担心默认定义了rate=0.4
                 Qthreshold=Qthreshold,
                 method=method,
                 SoP=SoP,
                 howGX=howGX,
                 Halve=Halve,
                 const=const,
                 endsize=endsize)  #返回rdat这个data.frame  ;  GetTree是最耗时的function
  
  ###get the MR est for each Nindex based on this single tree fitting----------------------
  ddat<-rdat #ddat是用来做Node/stratum-specific MR estimates
  NNindex<-as.numeric(  levels(     factor(   rdat$Nindex  )   )   )  #所有NNindex的种类
  #MRest vector  (same length and order as NNindex)
  MRest<-rep( NA, length( NNindex ) ) #NNindex 和MRest的组合构成了tree的end node所有结果(decision rule由rdat提供)
  Bx<-c();Bxse<-c();By<-c(); Byse<-c() #为了分析stratum-specific G-X effect
  Sizeratio<-c()
  #此处注意下estdata；如果是odat，那么rdat是自带Nindex；如果是用estdata，那么先算一下estdata$Nindex的Nindex，然后替换rdat为estdat后续代码是一样的
  if(honest){
    #estdat Nindex
    estdat$Nindex<-GetNindex(estdat[,5:(5+JJ-1)] ,rdat )#基于刚fit好的rdat来给estdat赋Nindex
    ddat<-estdat
  }
  
  
  ###node/stratum-specific MR estimates----------------------------------------------------------------
  
  #注意，如果时ddat 为rdat本身，那么肯定不会出现empty node for any Nindex,
  #而如果是estdat作为ddat，则根据decision rule的差异性会导致一些Nindex不存在样本，即无法估计出this Nindex-specific MR estimates
  if(howGX=='const'){
    if(is.na(const)){stop('You used a constant GX association, please tell me the constant value by const=my.const ' ) }
    bx<-const;bxse<-0
    for(nn in 1:length(NNindex)){
      dat_sub<-ddat[abs(ddat$Nindex- NNindex[nn])< .Machine$double.eps^0.5,]#focus on the Nindex-specific subgroup
      if(nrow(dat_sub)!=0){
        fitGY<-lm(    dat_sub[,4]~  dat_sub[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
        Bx<-c(Bx  , bx   ) ;  Bxse<-c(Bxse,  bxse  ) ; By<- c(  By , by ) ; Byse<-c(  Byse, byse )
        Sizeratio<-c( Sizeratio ,nrow( dat_sub) /nrow(ddat)   )  #for ts1  #注意honest style不需要做变动
        MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
        MRest[nn]<-MRres@Estimate
      }else{
        Sizeratio<-c( Sizeratio ,0   )
        MRest[nn]<-NA  #如果当前Nindex没有samples (只可能对于estdata出现)，那么就先设置为NA，最后再借其他Nindex的结果
      }
    }
  }else{
    for(nn in 1:length(NNindex)){
      dat_sub<-ddat[abs(ddat$Nindex- NNindex[nn])< .Machine$double.eps^0.5,]
      if(nrow(dat_sub)!=0){
        fitGX<-lm(    dat_sub[,3]~  dat_sub[,2]  )  ; bx<-as.numeric( summary(fitGX)$coef[-1,1]  ); bxse<-as.numeric(  summary(fitGX)$coef[-1,2])
        fitGY<-lm(    dat_sub[,4]~  dat_sub[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
        Bx<-c(Bx  , bx   ) ;  Bxse<-c(Bxse,  bxse  ) ; By<- c(  By , by ) ; Byse<-c(  Byse, byse )
        Sizeratio<-c( Sizeratio ,nrow( dat_sub) /nrow(ddat)   )  #for ts1
        MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
        MRest[nn]<-MRres@Estimate
      }else{
        Sizeratio<-c( Sizeratio ,0   )
        MRest[nn]<-NA  #如果当前Nindex没有samples (只可能对于estdata出现)，那么就先设置为NA，最后再借其他Nindex的结果
      }
    }
  }#end of N-index specific MR estimAte results -> MRest
  
  
  ###results calculation---------------------------------------------------------------------
  ###----------------------------------------------------------------------------------------
  RES<-list()
  
  ###RESULT 0: NNindex and MRest
  #如果MRest存在NA，那么就得借用相邻的end node (Nindex)  #注意，NNindex是严格升序的
  if( sum(is.na(MRest))>0  ){
    na_position<-which( is.na(MRest) )
    for(ppp in 1:length(na_position)  ){
      NNindex_<-NNindex
      NNindex_[na_position]<-1123#give  large value  #注意应该让所有NA的位置都赋很大的值；否则很容易出现两个NA相互借力的情况
      MRest[na_position[ppp]]<-MRest[ which.min( abs(NNindex_-NNindex[na_position[ppp]])  ) ]
    }
  }
  
  RES$end_node_information<-rbind(NNindex,MRest,Sizeratio)
  ##注意：NNindex和MRest构成了重要的tree end-node information 即使MRest可能是由estdat估计出来的
  
  ###RESULT 1: OOB and testing set predicts
  #OOB individual predicted data
  theMRest<- MRest[match(  as.character(GetNindex(  OOBdat[,5:( JJ+4 ) ] ,rdat,S=S ))  ,  as.character(NNindex)    )  ] #as.character  很重要的操作！ 防止数值误骗
  
  vect<-rep(0,NNN) #NNN 为odat size (无论是否用honest style，这么做都是没问题的)
  vect[  (Odat$I)%in%(OOBdat$I ) ]<-theMRest
  RES$OOB_predict<-vect
  
  #testing set individual predicted predicted values
  #(如果没有testing set，那么就把v_predict设空)
  #Vdat=vdat; vdat来自于当前环境(不一定是全局环境)
  if(as.matrix(is.na(Vdat))[1,1]){  #exists('vdat')
    RES$v_predict<-NA
  }else{
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
  ts1<-  sum(Sizeratio*(MRest-  sum(  Sizeratio*MRest  ))^2)#MRest不可能存在NA
  ##ts2 -> Q statistic:
  Sr<-mr_ivw(mr_input(Bx, Bxse , By, Byse))#先算effect via naive IVW under null (no effect difference)
  delta<-Sr@Estimate #一维的; 已验证！就是y_<-(By/Byse) ;  x_<-1/Byse; lm( y_ ~-1+x_)的结果
  #Cochran's Q statistic (more naive version - as we ignore the correlation of the GX estimate and GY estimate)
  ssigma_square <- Byse^2 + delta^2*Bxse^2
  ts2<-sum( ( By - delta*Bx     )^2 /ssigma_square      )#critical value:  qchisq(0.95, NC-1)==3.841459
  RES$ts1<-ts1
  RES$ts2<-ts2
  
  return( RES )
}

#RES list:  $end_node_information $OOB_predict $v_predict $vi1 $vi2 $ts1 $ts2
#$OOB_predict: index顺序和odat一致
#$v_predict: index顺序和vdat一致

#$vi1 和 $vi2 是按照1：JJ即covariate的顺序


# 
# ###examples:
# 
# set.seed(60)
# res<-getDat() #simulated data  #the deflaut setting: scenario='A' and SoM=0.5
# odat<-res$traning.set  #training set
# vdat<-res$testing.set  #testing set
# 
# #When running RFQT with mutiple Q trees/bootstrap - use parallel computation
# Nb<-5 #how many trees in the forest? e.g. 5 trees
# cl<-makeCluster(detectCores()-1)#规定用多少并行clusters/线程
# clusterEvalQ(cl=cl , expr=library(dplyr))  #给各个cluster中的运行一些表达式expression (例如导入一些包；做一些基础操作)
# clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
# clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat',
#                                     'GetTree', 'GetNindex', 'GetIndex' )  )#or any other arguments
# RES<-parSapply(   cl ,  1:Nb, BootstrapTreeFitting  ) #本脚本里parSapply的结果常用RES表示
# stopCluster(cl)
# dim(RES)#7 Nb
# 
# ##If you wish to use your own parameters rather than the default parameters, try:
# #general exmaple
# user_BootstrapTreeFitting<-function(seed){
#   RES<-BootstrapTreeFitting(seed,
#                             honest=my.honest,
#                             S=my.S,
#                             JJ=my.JJ,   #or any partial of the arguments
#                             rate=my.rate,
#                             Qthreshold=my.Qthreshold,
#                             method=my.method,
#                             SoP=my.SoP,
#                             howGX=my.howGX,
#                             endsize=my.endsize)
#   return(RES)
# }
# 
# 
# #specific exmaple
# user_BootstrapTreeFitting<-function(seed){
#   RES<-BootstrapTreeFitting(seed,SoP=20)
#   return(RES)
# }
# 
# ###parallel computation  -> RFQT
# Nb<-5 #how many trees in the forest?
# cl<-makeCluster(detectCores()-1)#规定用多少并行clusters/线程
# clusterEvalQ(cl=cl , expr=library(dplyr))  #给各个cluster中的运行一些表达式expression (例如导入一些包；做一些基础操作)
# clusterEvalQ(cl=cl , expr=library(MendelianRandomization) )
# clusterExport(  cl=cl ,  varlist=c( 'odat', 'vdat',
#                                     'GetTree', 'GetNindex', 'GetIndex' , 'BootstrapTreeFitting')  )
# RES<-parSapply(   cl ,  1:Nb, user_BootstrapTreeFitting  ) #本脚本里parSapply的结果常用RES表示
# stopCluster(cl)
# 
