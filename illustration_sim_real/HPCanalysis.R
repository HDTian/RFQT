
###post-HPC analysis

#three results: ATE DFfit and RFQTfit with the same seed index: 
#set.seed(10*i*100+iii) : i =1:6 for SOM and iii=1:100 for iterations



#for ATE:
#----------------------------------------------------------------------------------
Gammas<-c( 0, 0.1, 0.2, 0.3 ,0.4 ,0.5    )#共6种SoM情况

ggATE<-c()  #nrow(ggATE) = 3*6*100=1800

for(iii in 1:100){
  for(s in c('A','B','C')){
    scenario_used<-s
    for(i in 1:length(Gammas)){
        set.seed(10*i*100+iii)  #set.seed(10*i)
        gamma_used<-Gammas[i]
        #get data
        Dat<-getDat(scenario=scenario_used, SoM=gamma_used)
        odat<-Dat$traning.set  #training set in 全局环境
        vdat<-Dat$testing.set  #testing set in 全局环境
        
        #ATE
        fitGX<-lm(    odat[,3]~  odat[,2]  )  ; bx<-as.numeric( summary(fitGX)$coef[-1,1]  ); bxse<-as.numeric(  summary(fitGX)$coef[-1,2])
        fitGY<-lm(    odat[,4]~  odat[,2]  )  ; by<-as.numeric( summary(fitGY)$coef[-1,1]  ); byse<-as.numeric(  summary(fitGY)$coef[-1,2])
        MRres<-mr_ivw(mr_input(bx, bxse, by, byse))
        
        ggATE<-rbind(ggATE ,  
                     c( mean(  (vdat$true_STE - MRres@Estimate )^2 )  , gamma_used , scenario_used    )
        )
    }
  }
}

dim(ggATE)  #1800 3
write.csv(ggATE, "C:\\Users\\Haodong\\Documents\\RFQT\\RFQT\\Data\\ggATE.csv"  , row.names=FALSE) #Matrix.b是个矩阵
#ggATE<-read.csv( "C:\\Users\\Haodong\\Documents\\RFQT\\RFQT\\Data\\ggATE.csv" )


ATEmatrix<-matrix( as.numeric(ggATE[,1]), nrow=18, ncol=100  ) #
#ggdata 
ATEggdata<-data.frame( MSE=apply(  ATEmatrix, 1  , median), Gammas=as.numeric(ggATE[1:18,2]), Scenario=ggATE[1:18,3]  )
#ATEggdata<-data.frame( MSE=apply(  ATEmatrix, 1  , mean), Gammas=as.numeric(ggATE[1:18,2]), Scenario=ggATE[1:18,3] )



#for DTfit:
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

#in HPC terminal: $cat vector_newDT_*.txt > newDT.txt  

data0<-read.csv(  'C:\\Users\\Haodong\\Desktop\\newDT.txt' , header = FALSE, sep = "", dec = ".")
okID<-sort(data0[,1])
(1:5400)[!(1:5400)%in%okID]#the killed array task ID  #integer(0)

DTdata<-read.csv(  'C:\\Users\\Haodong\\Documents\\RFQT\\RFQT\\Data\\newDT.txt'   , header = FALSE, sep = "", dec = ".")
DTdata<-DTdata[order( DTdata[,1]  ),]#按照array ID increasing的顺序来

DTdata<-DTdata[,-1]#remove the first column: array task ID

#variable translation
DTvector<-as.vector(  as.matrix(DTdata)  ) #now all elements are string

DTvector[DTvector=='DR']<- 'Doubly-ranked stratification'
DTvector[DTvector=='Residual']<- 'Residual stratification'
DTvector[DTvector=='N']<- 'Naive stratification'
DTvector[DTvector=='DT']<- 'Single tree'

DTdata<-matrix(DTvector , 5400, 5   )
dim(DTdata) #5400 5 


##median or mean calculation

DTmatrix<-matrix( as.numeric(DTdata[,1]), nrow=54, ncol=100  ) #


#ggdata 
DTggdata<-data.frame( MSE=apply(  DTmatrix, 1  , median), Gammas=as.numeric(DTdata[1:54,2]), Scenario=DTdata[1:54,3] , de_collider=DTdata[1:54,4], method_type=DTdata[1:54,5] )
#DTggdata<-data.frame( MSE=apply(  DTmatrix, 1  , mean), Gammas=as.numeric(DTdata[1:54,2]), Scenario=DTdata[1:54,3] , de_collider=DTdata[1:54,4], method_type=DTdata[1:54,5] )


#for RFQTfit:
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

#in HPC terminal: $cat vector_newRFQT_*.txt > newRFQT.txt  

data0<-read.csv(  'C:\\Users\\Haodong\\Desktop\\newRFQT.txt' , header = FALSE, sep = "", dec = ".")
#check which ID is failed in HPC running and then add the running with the modified RFQT.sh with specific array task ID
okID<-sort(data0[,1])
(1:5400)[!(1:5400)%in%okID]#the killed array task ID
#in HPC terminal: $cat vector_newRFQT_*.txt > newRFQT_added.txt  
#then use these array ID to run ansd store the .txt and paste into newRFQT.txt

#then store into the RFQT Git path

RFQTdata<-read.csv(  'C:\\Users\\Haodong\\Documents\\RFQT\\RFQT\\Data\\newRFQT.txt'   , header = FALSE, sep = "", dec = ".")

#ID check
(1:5400)[!(1:5400)%in%RFQTdata[,1]]  #integer(0)
RFQTdata<-RFQTdata[order( RFQTdata[,1]  ),]#按照array ID increasing的顺序来

RFQTdata<-RFQTdata[,-1]#remove the first column: array task ID

#variable translation
RFQTvector<-as.vector(  as.matrix(RFQTdata)  ) #now all elements are string

RFQTvector[RFQTvector=='DR']<- 'Doubly-ranked stratification'
RFQTvector[RFQTvector=='Residual']<- 'Residual stratification'
RFQTvector[RFQTvector=='N']<- 'Naive stratification'
RFQTvector[RFQTvector=='RFQT']<- 'Random forest of Q trees'

RFQTdata<-matrix(RFQTvector , 5400, 5   )
  
dim(RFQTdata) #5400 5 



##median or mean calculation

RFQTmatrix<-matrix( as.numeric(RFQTdata[,1]), nrow=54, ncol=100  ) #
#ggdata 
RFQTggdata<-data.frame( MSE=apply(  RFQTmatrix, 1  , median), Gammas=as.numeric(RFQTdata[1:54,2]), Scenario=RFQTdata[1:54,3] , de_collider=RFQTdata[1:54,4], method_type=RFQTdata[1:54,5] )
#RFQTggdata<-data.frame( MSE=apply(  RFQTmatrix, 1  , mean), Gammas=as.numeric(RFQTdata[1:54,2]), Scenario=RFQTdata[1:54,3] , de_collider=RFQTdata[1:54,4], method_type=RFQTdata[1:54,5] )




###final visualization
ggdata<-rbind(RFQTggdata , DTggdata)
ggATE<-ATEggdata
#log scale for MSE
ggdata$MSE<-log(1+ggdata$MSE)
ggATE$MSE<-log(1+ggATE$MSE)

library(facetscales)
pd <- position_dodge(0.005)
scales_y <- list(
  `A` = scale_y_continuous(limits = c(0, 1.5)),
  `B` = scale_y_continuous(limits = c(0, 8)),
  `C` = scale_y_continuous(limits = c(0, 85))
)


#log-scale version
scales_y <- list(
  `A` = scale_y_continuous(limits = c(0, 1.2)),
  `B` = scale_y_continuous(limits = c(0, 2.5)),
  `C` = scale_y_continuous(limits = c(0, 5))
)

p<- ggplot(data = ggdata, mapping = aes(x = Gammas, y = MSE))+
  geom_line(mapping = aes(color= de_collider , linetype = method_type),size=1.0,alpha=1.0,position = pd)+
  geom_point(mapping = aes(color= de_collider ,  shape=de_collider ,  group = de_collider ) ,size = 2.5,alpha=1.0,position = pd)+
  geom_line(  data=ggATE ,  mapping = aes(x = Gammas, y = MSE), linetype=3 ,linewidth = 1.0 )+
  geom_point(  data=ggATE ,  mapping = aes(x = Gammas, y = MSE),size=2 )+
  labs(x = "Strength of modification", y = "MSE")+
  facet_grid_sc(rows=vars(Scenario), scales = list(y = scales_y))+
  theme_bw()+#theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  guides(  color=guide_legend(title='Stratification Method'),
           shape=guide_legend(title='Stratification Method'),
           group=guide_legend(title='Stratification Method')  , linetype=guide_legend(title='Estimation Method')  )


p<- ggplot(data = ggdata, mapping = aes(x = Gammas, y = MSE))+
  geom_line(mapping = aes(color= de_collider , linetype = method_type),size=1.0,alpha=1.0,position = pd)+
  geom_point(mapping = aes(color= de_collider ,  shape=de_collider ,  group = de_collider ) ,size = 2.5,alpha=1.0,position = pd)+
  geom_line(  data=ggATE ,  mapping = aes(x = Gammas, y = MSE), linetype=3 ,linewidth = 1.0 )+
  geom_point(  data=ggATE ,  mapping = aes(x = Gammas, y = MSE),size=2 )+
  labs(x = "Strength of modification", y = "MSE")+
  facet_grid(rows=vars(Scenario),scales = "free_y")+
  scale_y_continuous(labels = function(y) sprintf("%.2f", exp(y) - 1))+
  theme_bw()+#theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  guides(  color=guide_legend(title='Stratification Method'),
           shape=guide_legend(title='Stratification Method'),
           group=guide_legend(title='Stratification Method')  , linetype=guide_legend(title='Estimation Method')  )  

p

ggsave(paste0('Fig3.eps' ), 
       plot =p ,  #非指定 
       path=paste0('C:\\Users\\Haodong\\Desktop\\Precision_Medicine_new\\newplots'), 
       height = 7, width = 8, units = "in",limitsize=TRUE)









