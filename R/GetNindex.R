###
###GetNindex: get Node index for new input samples
###

#依据inputed rdat给any M 返回Nindex (Nindex: Node index)  (适用于任何M Input matrix)
#input either single individual (M is a vector) or multiple individual (M is a matrix)
GetNindex<-function(M,  #variable information matrix; the M colnumber and order should be same as the training data
                    rdat,#the result from GetTree() #it contains the tree information
                    S=NA #the same depth of the fitted tree
                    ){
  #先根据rdat把S算出来
  if(is.na(S)){ #如果不申明S，那么就会自动判断；如果不放心，就自己设置对的S
    namestring<-names( rdat )[ncol(rdat)]
    S<-as.numeric(  (strsplit(   namestring  , '_' )[[1]])[2]  )
    }else{S<-5  }
  #基本思路是把这个M也跟着tree pass down一下；每次对M进行Nindex的迭代/更新
  theNindex<-rep(0, nrow(M))  #此时是个vector
  for(s in 1:S){  #S: total split times; i.e. the same depth of the fitted tree
    Ls<-as.numeric( levels(factor(      round(  rdat$Nindex,s-1 )    )  ) )#rdat和dat是匹配的！
    for(k in 1:length(Ls)  ){
      rdat_sub<-rdat[ round(  rdat$Nindex,s-1 )  == Ls[k],   ]
      sstring<-strsplit(    levels(factor(     rdat_sub[ ,ncol(rdat  )-S+s ] ) )  , '_' )#splited string
      if(  length(sstring  )==0 ){
        index_added<-0
      }else{
        J_<-as.numeric( (sstring[[1]])[1] )  #此时关注的J_ th M
        cuttingvalue<-(as.numeric( (sstring[[2]])[3] )+as.numeric( (sstring[[1]])[3] ))/2
        M_sub<-M[ round( theNindex,s-1  ) == Ls[k] ,]  #只考虑满足这一类的sub matrix
        #vector赋值(赋vector Nindex)
        index_added<-((M_sub[,J_]>cuttingvalue)+1)*0.1^s #当然是个vector
      }
      theNindex[ round( theNindex,s-1  ) == Ls[k]     ]<-Ls[k]+index_added
    }
  }
  return( theNindex  ) #返回这个vector 结果
}

###examples:
set.seed(60)
res<-getDat() #simulated data
odat<-res$traning.set  #training set
vdat<-res$testing.set  #testing set

rdat<-GetTree(odat)


vdat_Nindex<-GetNindex(vdat[,5:24] ,rdat )#the M colnumber and order should be same as the training data

