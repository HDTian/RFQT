##DTfit: Decision Tree fit and its results

#similar BootstrapTreeFitting but without bootstrap (all related metrics will be removed) and use all candidate covariates



DTfit<-function(Odat=odat,  #inputted training data
                Vdat=NA,  #inputted testing data  #Vdat can be empty or nonexistent, not defined
                honest=FALSE, #use honest estimation or not?
                S=5,  #maximal depth
                rate=1, # the proportion of candidate variables Ms considered #if SingleM=TRUE, rate will be ignored
                SingleM=FALSE, #whether to ue single stratification style?(i,e, use a fixed M determined in the begining)
                Qthreshold=3.84, ##the threshold for Q heterogneity assessment
                method='DR',#stratification method used: 'DR' 'Residual' others
                SoP=10, ##size of pre-stratum #only make sense to DR stratification
                howGX='SpecificGX',##'const' means use extra constant; otherwise estimated by stratum data (stratum-specific GXeffect)
                Halve=FALSE, #if only use half splitting? Default is FALSE, i.e. use three splitting possiblies: (3:7) (5:5) (7:3)
                endsize=1000,##the minimal size of the node of Q-tree allowed to exist
                const=NA){
  if( is.null(Odat$true_STE[1]) ){JJ<-ncol( Odat )-4}else{JJ<-ncol( Odat )-5}
  if( JJ<1 ){ stop('No candidate covariate, or the dat is not regonized')  }
  
  #check the IDindex in Odat is strictly increasing
  if(  sum((Odat$I[-1]-Odat$I[-length(Odat$I)])<=0)>0  ){
    stop('the ID index of the inputed Odat is not strictly increasing: please first arrange them by yourself')
  }
  
  bdat<-Odat #不需要bootstrap
  dat<-bdat  #dat 是最终喂给getTree的data
  
  
  if(honest){ #honest estimation
    #random split (注意：必须要random split)
    rm_index<-sort(   sample( 1:nrow(bdat) , round(nrow(bdat)/2)   )   )
    treedat<-bdat[rm_index,]#bdat的I是严格升序的
    estdat<-bdat[-rm_index,]
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
  
  rdat<-GetTree( dat,  #把dat作为GetTree的输入  #dat可以为Odat或者的treedat
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
  
  
  ###RESULT 1:testing set predicts (没有OOB set predicts) 
  #testing set individual predicted predicted values
  #(如果没有testing set，那么就把v_predict设空)
  if(as.matrix(is.na(Vdat))[1,1]){
    RES$v_predict<-NA
  }else{
    theMRest<- MRest[match(  as.character(GetNindex(  Vdat[,5:( JJ+4 ) ] ,rdat,S=S ))  ,  as.character(NNindex)    )  ]
    RES$v_predict<-theMRest
  }
  
  ###RESULT 2: Variable Importance (VI) - 目前只定义于OOB data, 所以暂时不需要
  
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
  
  
  ###RESULT 4: testing MSE
  #(如果没有testing set，那么就把MSE设空)
  if(as.matrix(is.na(Vdat))[1,1]){
      RES$MSE<-NA   }else{
        if(  is.null(  Vdat$true_STE   ) ){  #p判断是否有label   
          RES$MSE<-NA    }else{
            RES$MSE<-  mean( ( theMRest- Vdat$true_STE   )^2   )
          }
      }
  
  
  
  return( RES )
}


#RES包括：$end_node_information  $v_predict  $ts1  $ts2 $MSE
#不需要getMSE算MSE，直接自己算即可


###examples:
set.seed(60)
res<-getDat() #simulated data  #the default setting: scenario='A' and SoM=0.5
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set

DTRES<-DTfit(odat,vdat)

DTRES

#get (testing set) MSE manually
mean( ( DTRES$v_predict - vdat$true_STE   )^2   )



