#getDat


#按照论文里的Model生成一组well-regonized standard data

#就是原始代码里的scenarioA,即论文里的scenarioB
getDat<-function(N=150000,
                Nt=100000,
                Nc=20, #Nc: number of candidate variables
                SoM=0.5, #SoM: strength of modification
                ZXeffect=0.5,  #instrument effect/strength
                XMeffect=0.5,  #the effect of X on covariates (for 2, 4, 6, ...)
                split=TRUE  #traning-testing split or not?

){
  gamma<-SoM
  NNNN<-N; NNN<-Nt
  Z<-rnorm(NNNN,0,1)
  e1<-rnorm(NNNN,0,1);e2<-rnorm(NNNN,0,1) #one for exposure, one for outcome
  JJ<-Nc #好的，现在我们直接考虑30个  #JJ: candidate variable的数量
  UU<-matrix(   rnorm( NNNN*JJ,0,1  ), NNNN ,JJ    )
  X<-ZXeffect*Z+UU%*%rep(0.5,JJ)+e1

  MM<-X%*%t(    rep(  c(0,XMeffect)  , JJ/2 )     )+UU

  modifier_vec<-c(  rnorm(5,gamma,0.1) , rep( 0,JJ-5 ) )    #前5个存在effect-modification
  Y<-(0.5+MM%*%modifier_vec   )*X+ UU%*%rep(0.5,JJ) +e2

  Dat<-cbind( 1:NNNN , Z , X ,Y ,MM  )
  The_true_STE<-(0.5+MM%*%modifier_vec   )#很关键的量！每次simulation都要根据当前的DAT新定义一遍

  #同一个Dat下，odat和vdat永远不会改变
  vdat<-Dat[(NNN+1):NNNN, ]  #validation dat
  odat<-Dat[1:NNN,]   #original (training) dat
  Dat<-as.data.frame(Dat)
  odat<-as.data.frame(odat);vdat<-as.data.frame(vdat)
  names(odat)[1:4]<-c('I','Z','X','Y')
  Dat$true_STE<-The_true_STE
  odat$true_STE<-The_true_STE[  1:NNN  ]
  vdat$true_STE<-The_true_STE[  (NNN+1):NNNN ]
  if(split == TRUE){
    return(   list( whole.data=Dat, traning.set=odat , testing.set = vdat   ))
  }else{
    return(Dat)
  }
}
