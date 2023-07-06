



###simulation code

###result1: (Figure 3)
#MSE results of the three X-M model Scenarios with different Strength of Modification (SoM) using RFQT and single-stratification


for(s in c('A','B','C')){
  scenario_used<-s
  Gammas<-c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    )#共6中SoM情况
  for(i in 1:length(Gammas)){
    set.seed(10*i)
    gamma_used<-Gammas[i]
    
    #get data
    Dat<-getDat(scenario=scenario_used, SoM=gamma_used)
    odat<-Dat$traning.set  #training set in 全局环境
    vdat<-Dat$testing.set  #testing set in 全局环境
    
    ##############ATE
    fitGX<-lm(    odat[,3]~  odat[,2]  )  ; bx<-as.numeric( summary(fitGX)$coef[-1,1]  ); bxse<-as.numeric(  summary(fitGX)$coef[-1,2])
    fitGY<-lm(    odat[,4]~  odat[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
    MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
    #mean(  (vdat$true_STE - MRres@Estimate )^2 ) 
    
    
    ###############Decision Tree
    ###DR
    DTRES<-DTfit(odat,vdat,Nb=200,method='DR')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_DR_',scenario_used,'_',gamma_used,'.RData'))
    #ALLRES<-readRDS(file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\S_DR_',scenario_used,'_',gamma_used,'.RData'))
    
    ###R
    DTRES<-DTfit(odat,vdat,Nb=200,method='Residual')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_R_',scenario_used,'_',gamma_used,'.RData'))
    
    ###N
    DTRES<-DTfit(odat,vdat,Nb=200,method='N')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\DT_N_',scenario_used,'_',gamma_used,'.RData'))
    
    
    ################RFQT
    ###DR
    ALLRES<-RFQTfit(odat,vdat,Nb=200,method='DR')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_DR_',scenario_used,'_',gamma_used,'.RData'))
    #ALLRES$MSE_test[length(ALLRES$MSE_test)] 
    
    ###R
    ALLRES<-RFQTfit(odat,vdat,Nb=200,method='Residual')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_R_',scenario_used,'_',gamma_used,'.RData'))
    
    ###N
    ALLRES<-RFQTfit(odat,vdat,Nb=200,method='N')
    saveRDS(ALLRES,file=paste0('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\RFQT_N_',scenario_used,'_',gamma_used,'.RData'))
    
  } #end of all strength of modification for a specific scenario
  
}






###result2: (Figure 7)
#MSE and individual predicted effect curve changing along the number of Q trees

set.seed(1123)
Dat<-getDat(scenario='A', SoM=0.5) #simulated data  #the default setting: scenario='A' and SoM=0.5
odat<-Dat$traning.set  #training set in 全局环境
vdat<-Dat$testing.set

ALLRES<-RFQTfit(odat,vdat,Nb=500,method='DR')#Nb=500

saveRDS(ALLRES,file='D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\sim500.RData')
#ALLRES<-readRDS('D:\\files\\R new\\Precison_Medicine\\ALLRES_rdata\\sim500.RData')
RES<-ALLRES$RES

MSE1<-getMSE( RES , 1 ) ; MSE2<-getMSE( RES , 2 )  #MSE1: OOB error; MSE2: test error
#getMSE()中odat vdat使用全局环境
plot( 1:length(MSE1) , MSE1, type='l',xlab='Number of Q trees', ylab='MSE' )
lines(     1:length(MSE2) , MSE2, col='blue'   )
legend('topright', legend=c('OOB error', 'Test error'),
       col=c("black", "blue"), lty=1, cex=1.0)
#Figure7_left  500*450

predict_matrix2<-getPredict( RES,2     )#2: test set predicts
plot(   1:ncol( predict_matrix2   )   ,  predict_matrix2[1,],  type='l', ylim=c(0.1,1.2),
        xlab='Number of Q trees', ylab='Predicted effect')
lines(     1:ncol( predict_matrix2   )   ,  predict_matrix2[2,], col='blue'   )
lines(     1:ncol( predict_matrix2   )   ,  predict_matrix2[3,], col='red'   )
lines(     1:ncol( predict_matrix2   )   ,  predict_matrix2[4,], col='green'   )

legend('topright', legend=c('Individual One', 'Individual Two', 'Individual Three', 'Individual Four'),
       col=c("black", "blue",'red','green'), lty=1, cex=0.7)
#Figure7_right 500*450


##ggplot version
#left
ggdata<-data.frame( trees= 1:length(MSE1) , MSE=c(MSE1,MSE2),   
                    type=c( rep( 'OOB error' , length(MSE1) ) , rep( 'Test error' , length(MSE2) )   )    )
p<-ggplot(ggdata,aes(trees, MSE,colour=type) ) +
  geom_line( linewidth=1.0 ) +
  xlab('Number of Q trees')+ylab('MSE')+
  geom_vline(xintercept = 200, colour='grey', linewidth=1.0 ,linetype = 2)+
  guides(colour = guide_legend(title = " "))+
  theme(panel.background = element_rect(linewidth=1,colour = "black"),
        legend.position = c(0.85, 0.85),legend.title = element_blank(),
                   legend.background = element_rect( linewidth=0.5, linetype="solid",colour = "black"),
                                                    panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p
ggsave(paste0('Fig7_left.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots'), 
       height = 5, width = 7, units = "in",limitsize=TRUE)

#right
ggdata<-data.frame( trees= 1:ncol( predict_matrix2)  , pe=c(predict_matrix2[1,],predict_matrix2[2,],predict_matrix2[3,],predict_matrix2[4,]),   
                    type=  factor(c( rep( 'Individual One' , ncol(predict_matrix2)  ) ,
                            rep( 'Individual Two' , ncol(predict_matrix2)  ) ,
                            rep( 'Individual Three' , ncol(predict_matrix2)  ) ,
                            rep( 'Individual Four' , ncol(predict_matrix2)  ) )   , 
                            levels=c('Individual One','Individual Two','Individual Three','Individual Four'  )  ) )
p<-ggplot(ggdata,aes(trees, pe,colour=type) ) +
  geom_line( linewidth=1.0 ) +
  xlab('Number of Q trees')+ylab('Predicted effect')+
  geom_vline(xintercept = 200, colour='grey', linewidth=1.0 ,linetype = 2)+
  guides(colour = guide_legend(title = " "))+
  theme(panel.background = element_rect(linewidth=1,colour = "black"),
        legend.position = c(0.85, 0.70),legend.title = element_blank(),
        legend.background = element_rect( linewidth=0.5, linetype="solid",colour = "black"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p
ggsave(paste0('Fig7_right.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots'), 
       height = 5, width = 7, units = "in",limitsize=TRUE)





##Results3: (newFig1)
#scatter plot of true HTE v.s. est HTE
plot(   ALLRES$Predicts_test[,500],vdat$true_STE  )
abline(0,1,col='red')

#ggplot
ggdata<-data.frame(estHTE=ALLRES$Predicts_test[,500], trueHTE=vdat$true_STE )
set.seed(1123)
ggdata<-ggdata[ sample( 1:nrow(ggdata),1000  ) , ]
p<-ggplot(ggdata) +
  geom_abline(intercept = 0, slope = 1, colour='grey', linewidth=1.0 ,linetype = 2)+
  geom_hline(yintercept = 0, colour='grey', linewidth=1.0 ,linetype = 2)+
  geom_vline(xintercept = 0, colour='grey', linewidth=1.0 ,linetype = 2)+
  geom_point(aes(estHTE,trueHTE),size = 1) +
  xlab('Predicted effects')+ylab('True effects')+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p
ggsave(paste0('newFig1.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots'), 
       height = 5, width = 5, units = "in",limitsize=TRUE)


##Results4: (newFig2)
#VI plot
names(odat)
names(odat)[5:24]<-paste0('M',1:20)
#left - VI1
VIindex<-order(ALLRES$VI1)#得出一组position index用以表示最小到最大unit的各自的位置 
ggdata_VI<-data.frame(candidate=factor(names(odat)[VIindex+4],levels =(names(odat)[VIindex+4])[length(VIindex):1] ),  #逆向排序
                      VI=( ALLRES$VI1 )[VIindex]   )  
p<-ggplot(ggdata_VI,aes(x=candidate,y=VI))+geom_col( )+
  xlab('Candidate covariate')+ylab('Variable Importance')+#labs(title=c('Doubly-ranked method', 'Residual method' , 'Naive method'   )[RESindex] )+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  coord_flip()
p
ggsave(paste0('newFig2_left.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots'), 
       height = 6, width = 4, units = "in",limitsize=TRUE)

#right - VI2
VIindex<-order(ALLRES$VI2)#得出一组position index用以表示最小到最大unit的各自的位置 
ggdata_VI<-data.frame(candidate=factor(names(odat)[VIindex+4],levels =(names(odat)[VIindex+4])[length(VIindex):1] ),  #逆向排序
                      VI=( ALLRES$VI2 )[VIindex]   )  
p<-ggplot(ggdata_VI,aes(x=candidate,y=VI))+geom_col( )+
  xlab('Candidate covariate')+ylab('Variable Importance')+#labs(title=c('Doubly-ranked method', 'Residual method' , 'Naive method'   )[RESindex] )+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  coord_flip()
p
ggsave(paste0('newFig2_right.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots'), 
       height = 6, width = 4, units = "in",limitsize=TRUE)



##Result5: (newFig3)
#according to vdat predicted HTE (i.e. v_predict), refit HTE ~ M to visualize the results or for future guidance
#ALLRES就用results2 里面的即可

#ALLRES<-RFQTfit(odat,vdat,Nb=50,method='DR',S=7,endsize=1000)
#RES<-ALLRES$RES
#predict_matrix2<-getPredict( RES,2     )
dim(predict_matrix2)
dim(vdat)

#simple linear regression
fit<-lm(  predict_matrix2[,ncol(predict_matrix2)]~ as.matrix(vdat[,5:24])     )
summary(fit)



#Decision tree
library(rpart)
library(rpart.plot)

dt_data<-cbind(vdat[,5:24],predict_matrix2[,500] ) #decision tree data

dim(dt_data)
names(dt_data)
names(dt_data)<-c(  paste0('M',1:20) , 'estHTE' )


tree <- rpart(estHTE ~. , data = dt_data,control = rpart.control(minsplit  = 2,cp=0.005))#cp可以控制tree的精度
rpart.plot(tree)#共12 leaf/end-node 
#newFig3_upper 650 500
#已验证： 图像里leaf/end node上的值就是各自的mean predicted HTE

###end node specific plot

dim( tree$frame)#31 8 #共12 rows with first col namevariable = <leaf>
which(   tree$frame[,1]=='<leaf>' )#5  6  8  9 12 13 15 16 20 21 23 24 27 28 30 31  #严格按照rpart.plot(tree)从左到右的顺序
length(tree$where) #the row number of frame corresponding to the leaf node that each observation falls into
table( tree$where )
#   5    6    8    9   12   13   15   16   20   21   23   24   27   28   30   31 
#3035 3135 3111 3152 3080 3191 3111 3169 3191 3094 3218 3074 3097 3166 3134 3042 

#rename: number -> ABCD...
ABCnames<-rep(NA,30)
ABCnames[  which(   tree$frame[,1]=='<leaf>' )   ]<-LETTERS[1:length( which(   tree$frame[,1]=='<leaf>' )  )]


ABCleaf<-ABCnames[    tree$where    ]


ggdata<-data.frame( trueHTE=vdat$true_STE , 
                    leaf= factor(   ABCleaf  ,   levels= LETTERS[1:length( which(   tree$frame[,1]=='<leaf>' )  )]       ) )
p<-ggplot(ggdata, aes(leaf, trueHTE))+
  geom_hline(yintercept = 0, colour='grey', linewidth=1.0 ,linetype = 2)+
   geom_boxplot()+
   xlab('Candidate covariate')+ylab('True effects')+
   scale_x_discrete(position = "top") +
   theme_bw()+theme(axis.title.x=element_blank(),
     panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p
ggsave(paste0('newFig3_down.eps' ),   #.eps  #real以后都特指R3： BMI on fev1 for the female
       plot =p ,  #非指定 
       path=paste0('C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots'), 
       height = 5, width = 7, units = "in",limitsize=TRUE)

