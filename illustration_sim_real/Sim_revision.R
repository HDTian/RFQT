
###supplememtary simualtion for revision


###getDat() for supplementary/revision; SoM as the argument; return Dat, which contans $traning.set and $testing.set 

getDat_supp<-function(N=150000,
                 Nt=100000,
                 Nc=20, #Nc: number of candidate variables
                 # scenario='A', #'A' 'B' or other character (i.e. 'C')
                 SoM=0.5, #SoM: strength of modification (i.e. gamma)
                 ZXeffect=0.5,  #instrument effect/strength
                 Random=TRUE, #if use the random value, N(0,0.1^2), for the strenght of modifcation
                 label=TRUE, #if to added the final column as the true HTE (i.e. label)
                 split=TRUE  #traning-testing split or not? and if store other information? if not, there is only one complete data and no other information
){
  gamma<-SoM
  NNNN<-N; NNN<-Nt
  ##Z: instrument vector
  Z<-rnorm(NNNN,0,1)
  e1<-rnorm(NNNN,0,1);e2<-rnorm(NNNN,0,1)#e1 for exposure, e2 for outcome, 
  JJ<-Nc   #JJ: candidate variable的数量
  
  ##UU: confounder matrix
  UU<-matrix(   rnorm( NNNN*JJ,0,1  ), NNNN ,JJ    ) #NNNN=150000 * JJ=20
  
  #ee: error term matrix for the covariate
  ee<-matrix(     rnorm( NNNN*JJ,0,1  )      ,NNNN,JJ   )
  
  
  ##first generate Mj for j=1:10
  MM1to10<- UU[,1:10] + ee[,1:10]
  
  
  ##X: exposure vector
  X<-ZXeffect*Z +  MM1to10%*%rep(0.5, 10)  + UU%*%rep(0.5,JJ)+e1 #length(X) #150000
  
  ##Y: the outcome vector
  if(Random){  
    modifier_vec<-c(  rnorm(5,gamma,0.1) , rep( 0,JJ-5 ) )    #前5个存在effect-modification (modification strength = gamma)
  }else{
    modifier_vec<-c(  rnorm(5,gamma,0) , rep( 0,JJ-5 ) )    #前5个存在effect-modification (modification strength = gamma)
  }
  
  Y<-( 0.5+MM1to10%*%modifier_vec[1:10]  )*X+ UU%*%rep(0.5,JJ) +e2
  
  ##generate Mj for j=11:20
  MM11to15<- UU[,11:15] + ee[,11:15] + as.vector(0.5*Y) + as.vector( 0.5*X)  #collider
  MM16to20<- UU[,16:20] + ee[,16:20] + as.vector(0.5*Y) + as.vector( 0.0*X)  #collider
  
  MM<-cbind(MM1to10, MM11to15,  MM16to20  )
  
  
  ##One data set
  Dat<-cbind( 1:NNNN , Z , X ,Y ,MM  )
  
  ##The_true_STE: HTE/STE vector (labels)
  The_true_STE<-( 0.5+MM1to10%*%modifier_vec[1:10]  )#很关键的量！每次simulation都要根据当前的DAT新定义一遍
  
  #同一个Dat下，odat和vdat永远不会改变
  vdat<-Dat[(NNN+1):NNNN, ]  #validation (testing) dat
  odat<-Dat[1:NNN,]   #original (training) dat
  #[trainging dat may be further divided into tree dat and estimation dat - honest estimation]
  Dat<-as.data.frame(Dat)
  odat<-as.data.frame(odat);vdat<-as.data.frame(vdat)
  names(Dat)[1:4]<-c('I','Z','X','Y')
  names(odat)[1:4]<-c('I','Z','X','Y')
  names(vdat)[1:4]<-c('I','Z','X','Y')
  if(label){
    Dat$true_STE<-The_true_STE
    odat$true_STE<-The_true_STE[  1:NNN  ]
    vdat$true_STE<-The_true_STE[  (NNN+1):NNNN ]
  }
  if(split == TRUE){
    return(   list( whole.data=Dat, traning.set=odat , testing.set = vdat, SoM=gamma , modifier_vec=modifier_vec    ))
  }else{
    return(Dat)
  }
}

#example
set.seed(60)
res<-getDat_supp() #simulated data #the deflaut setting: scenario='A' and SoM=0.5
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set



###get data
Gammas<-c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    )

ATEresults<-c()
Nresults<-c()
DRresults<-c()

for(i in 1:length(Gammas)){
  set.seed(10*i)
  gamma_used<-Gammas[i]

  ###get data
  Dat<-getDat_supp(SoM=gamma_used)  #SoM: strength of modification
  odat<-Dat$traning.set  #training set in 全局环境
  vdat<-Dat$testing.set  #testing set in 全局环境


  ###RFQT fitting and ATE fitting 

  ###ATE
  fitGX<-lm(    odat[,3]~  odat[,2]  )  ; bx<-as.numeric( summary(fitGX)$coef[-1,1]  ); bxse<-as.numeric(  summary(fitGX)$coef[-1,2])
  fitGY<-lm(    odat[,4]~  odat[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
  MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
    
  ATEresults<-c( ATEresults , mean(  (vdat$true_STE - MRres@Estimate )^2 )  )
  #3.057939


  
  ###Naive
  ALLRES<-RFQTfit(odat,vdat,Nb=200,method='N')  
  saveRDS(ALLRES,file=paste0('/Users/haodongtian/Documents/HTEdata/RFQT_N_','supp_',gamma_used,'.RData'))
  
  #MSE
  Nresults<-c( Nresults , ALLRES$MSE_test[length(ALLRES$MSE_test)] )
  #1.698725 for N

  ###DR
  ALLRES<-RFQTfit(odat,vdat,Nb=200,method='DR')
  saveRDS(ALLRES,file=paste0('/Users/haodongtian/Documents/HTEdata/RFQT_DR_','supp_',gamma_used,'.RData'))
  
  #MSE
  DRresults<-c( DRresults , ALLRES$MSE_test[length(ALLRES$MSE_test)] )
  #1.248169 for DR
 
}  


###visualization
ggdata<-data.frame(  SoM=c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    ),
                     results= c( Nresults ,  DRresults, ATEresults  ) ,
                     type=  rep( c('Naive stratification','Doubly-ranked stratification' ,'No stratification') , each=6     ) )


p<- ggplot(data = ggdata, mapping = aes(x = SoM, y = results))+
  geom_line(mapping = aes(color= type ),linewidth=1.0,alpha=1.0)+
  geom_point(mapping = aes(color= type ) ,size = 2.5,alpha=1.0)+
  labs(x = "Strength of modification", y = "MSE")+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  guides(  color=guide_legend(title='Stratification Method')  )  
p




###or use HPC (100 independent simulations; each contains the 6*3=18 results, representing ATEresults , Nresults ,  DRresults  )

#only consider ATEresults and DRresults

RESmatrix<-c()
for(iii in 1:100){
  test <- read.table( paste0(  "/Users/haodongtian/Documents/HTEdata/HPCrevision/vector_",iii,".txt"  ), header = FALSE, sep = "\t")
  RESmatrix<-cbind( RESmatrix , test$V1  )
}

#get the median value
result_vector<-apply( RESmatrix, 1, median )

ggdata<-data.frame(  SoM=c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    ),
                     results= result_vector[-c(7:12)],
                     type=  rep( c('ATE estimation without stratification','RFQT with the doubly-ranked stratification' ) , each=6     ) )


p<- ggplot(data = ggdata, mapping = aes(x = SoM, y = results))+
  geom_line(mapping = aes(color= type ),linewidth=1.0,alpha=1.0)+
  geom_point(mapping = aes(color= type ) ,size = 2.5,alpha=1.0)+
  scale_color_manual(values = c("ATE estimation without stratification" = "black", 
                                "RFQT with the doubly-ranked stratification" = "red"))+
  labs(x = "Strength of modification", y = "MSE")+
  theme_bw()+theme(legend.position = "bottom")+
  guides(  color=guide_legend(title='Method')  )  
p
ggsave(paste0('Sim_MSE_supp.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('/Users/haodongtian/Documents/HTE_revision'), 
       height = 4, width = 7, units = "in",limitsize=TRUE)




###Variable Importance

#taking SoM=0.5 as the example to present the variable importance
ALLRES<-readRDS(paste0('/Users/haodongtian/Documents/HTEdata/RFQT_DR_','supp_',0.5,'.RData'))


ALLRES$VI1; ALLRES$VI2


names(odat)[5:24]<-paste0('M',1:20)

#left plot
VIindex<-order(ALLRES$VI1)#得出一组position index用以表示最小到最大unit的各自的位置 
ggdata_VI<-data.frame(candidate=factor(names(odat)[VIindex+4],levels =(names(odat)[VIindex+4])[length(VIindex):1] ),  #逆向排序
                      VI=( ALLRES$VI1 )[VIindex]   )  
p<-ggplot(ggdata_VI,aes(x=candidate,y=VI))+geom_col( )+
  xlab('Candidate covariate')+ylab('Variable Importance')+#labs(title=c('Doubly-ranked method', 'Residual method' , 'Naive method'   )[RESindex] )+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  coord_flip()
p
ggsave(paste0('sim_supp_VI_left.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('/Users/haodongtian/Documents/HTE_revision'), 
       height = 4, width = 4, units = "in",limitsize=TRUE)

#right plot
VIindex<-order(ALLRES$VI2)#得出一组position index用以表示最小到最大unit的各自的位置 
ggdata_VI<-data.frame(candidate=factor(names(odat)[VIindex+4],levels =(names(odat)[VIindex+4])[length(VIindex):1] ),  #逆向排序
                      VI=( ALLRES$VI2 )[VIindex]   )  
p<-ggplot(ggdata_VI,aes(x=candidate,y=VI))+geom_col( )+
  xlab('Candidate covariate')+ylab('Variable Importance')+#labs(title=c('Doubly-ranked method', 'Residual method' , 'Naive method'   )[RESindex] )+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  coord_flip()
p
ggsave(paste0('sim_supp_VI_right.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('/Users/haodongtian/Documents/HTE_revision'), 
       height = 4, width = 4, units = "in",limitsize=TRUE)






