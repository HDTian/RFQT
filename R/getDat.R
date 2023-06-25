#getDat


#按照论文里的Model生成一组well-regonized standard data

#就是原始代码里的scenarioA,即论文里的scenarioB
getDat<-function(N=150000,
                Nt=100000,
                Nc=20, #Nc: number of candidate variables
                scenario='A', #'A' 'B' or other character (i.e. 'C')
                SoM=0.5, #SoM: strength of modification (i.e. gamma)
                ZXeffect=0.5,  #instrument effect/strength
                split=TRUE  #traning-testing split or not? if not, there is only one complete data
){
  gamma<-SoM
  NNNN<-N; NNN<-Nt
  ##Z: instrument vector
  Z<-rnorm(NNNN,0,1)
  e1<-rnorm(NNNN,0,1);e2<-rnorm(NNNN,0,1) #one for exposure, one for outcome
  JJ<-Nc   #JJ: candidate variable的数量
  ##UU: confounder matrix
  UU<-matrix(   rnorm( NNNN*JJ,0,1  ), NNNN ,JJ    )
  ##X: exposure vector
  X<-ZXeffect*Z+UU%*%rep(0.5,JJ)+e1

  ##MM: candidate variable (Ms) matrix
  if(scenario=='A'){   #no collider
    #scenarioA
    MM<-X%*%t(    rep(  c(0,0)  , JJ/2 )     )+UU
  }else{
    if(scenario=='B'){  #half collider
      #scenarioB
      MM<-X%*%t(    rep(  c(0,0.5)  , JJ/2 )     )+UU   #这里主要区别是是否有collider (-Scenario B)  ; 注意后面的代码也要相应匹配！
    }else{   #all collider and also complicated collider
        #scenarioC
        MM<-as.vector(X)*(0.1+0.5*UU)+UU #直接用*即可，diag()对于150000这个量维度太大了！  注意之前的X是一个Matrix而不是vector!
    }
  }


  modifier_vec<-c(  rnorm(5,gamma,0.1) , rep( 0,JJ-5 ) )    #前5个存在effect-modification (modification strength = gamma)
  ##Y: the outcome vector
  Y<-(0.5+MM%*%modifier_vec   )*X+ UU%*%rep(0.5,JJ) +e2

  ##One data set
  Dat<-cbind( 1:NNNN , Z , X ,Y ,MM  )
  ##The_true_STE: HTE/STE vector (labels)
  The_true_STE<-(0.5+MM%*%modifier_vec   )#很关键的量！每次simulation都要根据当前的DAT新定义一遍

  #同一个Dat下，odat和vdat永远不会改变
  vdat<-Dat[(NNN+1):NNNN, ]  #validation (testing) dat
  odat<-Dat[1:NNN,]   #original (training) dat
  #[trainging dat may be further divided into tree dat and estimation dat - honest estimation]
  Dat<-as.data.frame(Dat)
  odat<-as.data.frame(odat);vdat<-as.data.frame(vdat)
  names(Dat)[1:4]<-c('I','Z','X','Y')
  names(odat)[1:4]<-c('I','Z','X','Y')
  names(vdat)[1:4]<-c('I','Z','X','Y')
  Dat$true_STE<-The_true_STE
  odat$true_STE<-The_true_STE[  1:NNN  ]
  vdat$true_STE<-The_true_STE[  (NNN+1):NNNN ]
  if(split == TRUE){
    return(   list( whole.data=Dat, traning.set=odat , testing.set = vdat   ))
  }else{
    return(Dat)
  }
}

##examples:

set.seed(60)
res<-getDat() #simulated data #the deflaut setting: scenario='A' and SoM=0.5
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set












