##DTfit: Decision Tree fit


DTfit<-function(seed=1,
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
                               endsize=1000,##the minimal size of the node of Q-tree allowed to exist
                               const=NA){ 
  
  
}