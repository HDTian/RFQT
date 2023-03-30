
#get the present values of all the relevant parameters and if not defined the default value 
getPars<-function(empty.argument){
  if( !exists('odat') ){ odat<-NA
  print('No training data correctly defined, please use odat')
  JJ<-'N/A'}else{
    if( is.null(odat$true_STE[1]) ){JJ<-ncol( odat )-4}else{JJ<-ncol( odat )-5}
  }
  if( !exists('vdat') ){ vdat<-NA
  print('No testing data correctly defined, please use vdat')}
  
  ##The 8 parameters
  if( !exists('method') ){  method<-'Not defined (DR)' }
  if( !exists('SoP') ){  SoP<-'Not defined (10)' }
  if( !exists('rate') ){  rate<-'Not defined (0.4)' }
  if( !exists('S') ){  S<-'Not defined (5)' }
  if( !exists('howGX') ){  howGX<-'Not defined (SpecificGX)' }
  if( !exists('endsize') ){  endsize<-'Not defined (5000)' }
  if( !exists('Qthreshold') ){  Qthreshold<-'Not defined (3.0)' }
  
  if( howGX=='const' ){
    if( !exists('const') ){  const<-'Not defined' }
   }else{const<-'N/A'  }
  

  
  results<-c( JJ, nrow(odat), nrow(vdat) , method , SoP , rate, S , howGX, const, endsize,Qthreshold )
  names(results)<-c( 'number.of.candidates.covariate',
                     'training.data.size',
                     'testing.data.size',
                     'stratification.method',
                     'size.of.pre-stratum',
                     'random.proportion',
                     'max.tree.deep' ,
                     'instrument-exposure.style',
                     'GXeffect.value',
                     'min.end.node.size',
                     'Q.value.threshold')
  return(results)
}

