#gmlSEM is under development
#Optionalities are listed below and will be activated upon progression in development

#Rules:
# Runtime: . within the variable names will be eliminated if the variable is not found within the dictionary nor the data set

#Under development:
implementations=list(
  #Multi-Level analysis
  MultiLevel=list(FALSE,test=function(){
    
  }),
  #Inflated distributions
  InflatedDistributions=list(FALSE,test=function(){
    
  }),
  #Truncated distributions
  TruncatedDistributions=list(FALSE,test=function(){
    
  }),
  #Censored distributions
  CensoredDistributions=list(FALSE,test=function(){
    
  }),
  #Categorized distributions
  CategorizedDistributions=list(FALSE,test=function(){
    
  }),
  #Copula Functions
  Copula=list(FALSE,test=function(){
    
  }),
  #Covariance terms
  Covariance=list(FALSE,test=function(){
    
  }),
  #Constraints
  Constraints=list(FALSE,test=function(){
    
  }),
  #Random Effects in regression models
  RandomEffects=list(FALSE,test=function(){
    
  }),
  #functions within regression model i.e., y~center(x)+scale(y)
  ScaleTermRegression=list(FALSE,test=function(){
    
  }),
  #Multi-Group analysis
  MultiGroups=list(FALSE,test=function(){
    
  }),
  #Group can accept finite support random variables
  LatentGroups=list(FALSE,test=function(){
    
  }),
  #Vector-valued random variables , e.g. multinom, dirichlet, tobit of type III, IV, V and etc.
  VectorValuedFamilies=list(FALSE,test=function(){
    
  }),
  #Smooth effects, e.g. y~pspline(x)+pspline(x,y)
  SmoothEffects=list(FALSE,test=function(){
    
  })
  
)
