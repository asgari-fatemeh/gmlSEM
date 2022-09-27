##Arguments
#internal machinery
data.is.provided = FALSE

#Main object returned by the parser
parsedData=list()

#start.idx <- grep("[~=<>:|%]", model.simple)
operators.blocks <- c(
  #G: gmlSEM supported operator
  #L: lavaan supported operator
  
  ## gmlSEM specific ancillary operators 
  # Alias is the first operator to look at
  "<>"    = "alias"     , #G <> substitute for as keyword in aliases
  "<<~"   = "vary"      , #G <<~ substitute for 'vary at level' clause
  "<<"    = "level"     , #G << substitute for within keyword in level: blocks, 
  ">>F>>" = "family"    , #G identify distribution and copula in family: blocks
  ">>G>>" = "group"     , #G Group block
  ">>S>>" ="size"       , #G Sugar syntax for specifying sampe size for simulation prupose
  "^\\w+:" = "block"    ,
  
  ## Other operators
  "=~"   = "measurement", #L#G Latent variable definitions:
  "<~"   = "formative"  , #L#G formative factors as in f5 <~ z1 + z2 + z3 + z4 
  "~~"   = "covariance" , #L#G Covariance
  "~"    = "regression" , #L#G Regression model
  ":="   = "monitor"    , #L#G Defining parameters to monitor
  "<="   = "constraint" , #G   constraint
  ">="   = "constraint" , #G   constraint
  "<"    = "constraint" , #L#G constraint
  ">"    = "constraint" , #L#G constraint
  "=="   = "constraint" , #L#G constraint
  "="    = "constraint"   #G   constraint
  #":"   = "block"         gmlSEM does not expect any other block
  #"\\|" = "nouse"         gmlSEM provides family: blocks as a unified way to define distributions forall variables.
  #                        Hence, we do not use \\| for defining thresholds for categorical variables as in lavaan.
  #                        family: y ordinal(levels=c(1,2,3,4),thresholds=c(0,NA,NA),underlying.family=Gaussian())
  #"~*~" = "nouse"         Scaling factor operator is ommited from gmlSEM syntax
  #"%"   = "nouse"
)

operators <- names(operators.blocks)
lhs.blocks <- c("level","family","heter","group","size")

GROUP_OP <- FALSE
LEVEL_OP <- FALSE
group<-""

add.vars.to.dictionay<-function(name="",vars.list.1=c(),vars.list.2=c(),family=NULL,copula=NULL){
  if(name!="")
    name=get.alias.rhs(name)
  if(length(vars.list)>0){
    vars.list=get.alias.rhs(vars.list)
    N=attr(vars.dictionary,"N")+1
    attr(vars.dictionary,"N")=N
    
  }
    
  ind=which(vars.dictionary)
}

#Parser Output
vars.dictionary <- data.frame(
  name         = character() , 
  var.list     = I(list())   , #The list of variables for vector-valued r.v. or repeated measures
  appeared.as.fa = logical() , #Appeared as a factor in a measrement model
  appeared.as.re = logical() , #Appeared as a random effect in a regression model
  latent       = logical()   , #NA if needs.data.scan
  response     = logical()   ,
  family       = I(list())   ,
  vary         = character() ,
  heter        = logical()   ,  # A heteroskedasticity model?
  reg.model    = character() ,
  heter.params = I(list())   ,
  heter.vars   = I(list())   ,
  reg.params   = I(list())   ,
  reg.vars     = I(list())   ,
  reg.smooth   = logical()   ,
  depends.on   = I(list())
)

attr(vars.dictionary,"N")=0

all.alias       <- data.frame(
  lhs=character(),
  rhs=character(),
  var.name=character(),
  role=character(), #PAR, LV, OV
  all.forms=I(list())
)

covs            <- list()  #lists of lists. one cov matrix per level per group
levels          <- list()
levels.matrix   <- matrix(0,1,1)
rownames(levels.matrix) <- colnames(levels.matrix) <- "1"


#gmlSEM parser environment to resolve parameter values
env=new.env()

group.specs = list(
  list(
    vars.dictionary = vars.dictionary,
    covs            = covs #Group specific covariance structure
  )
)


set.vary<-function(lhs,rhs){
  vars<-split.vars(lhs)
  for(i1 in seq_along(vars)){
    var<-vars[i]
    if(! var %in% vars.dictionary$name)
      vars.dictionary[nrow(vars.dictionary)+1,1]<-var
    vars.dictionary$vary[vars.dictionary$name==var]<-rhs
  }
}

add.levels<-function(lev){
  
  if(length(lev)>1){
    for(i in seq_along(lev))
      add.levels(lev[i])
    return()
  }
  
  if(lev=="")
    return()
  
  lev=trim(lev)
  lev=get.alias.rhs(lev)
  
  # if(lev=="1")
  #   return()  #Do not add the base level to the matrix
  
  rn=get.alias.rhs(rownames(levels.matrix))
  if(lev %in% rn)
    return()
  
  
  levels.matrix<<-rbind(levels.matrix,0)
  levels.matrix<<-cbind(levels.matrix,0)
  colnames(levels.matrix)<<-rownames(levels.matrix)<<-c(rn,lev)
  
}

set.levels.matrix<-function(lhs,rhs){
  
  if(length(lhs)>1){
    for(i in 1:length(lhs))
      for(j in 1:length(rhs))
        set.levels.matrix(lhs[i],rhs[j])
    return()
  }
  
  if(lhs==""||rhs=="")
    return()
  lhs<-get.alias.rhs(lhs)
  if(lhs=="1")
    return()  #ignore the base level
  rhs<-get.alias.rhs(rhs)
  if(rhs=="1")
    stop("gmlSEM error: base level can not be superior to other levels.")
  
  add.levels(rhs)
  levels.matrix[lhs,rhs]<-1
}

get.alias.ind<-function(lbl){
  if(nrow(all.alias)==0)
    return(integer(0))
  inds<-apply(all.alias,1, function(x) lbl%in% x[['all.forms']])
  which(inds)
}

#Reconsider relabeling ability of unlabeled symbols a=a -> a=b
#Also add labels to existing alias collections
add.alias<-function(lhs,rhs=NULL,role=NA){
  
  if(length(lhs)>1){
    if(length(rhs)==1)
      rhs=rep(rhs,length(lhs))
    if(length(role)==1)
      role=rep(role,length(lhs))
    
    for(i in seq_along(lhs))
      add.alias(lhs[i],rhs[i],role[i])
    
    return()
  }
  
  if(is.null(rhs))
    rhs=lhs
  
  if(!is.na(role))
    role=match.arg(role,c("PAR","OV","LV"))
  
  # if(is.null(rhs)){
  #   inds<-get.alias.ind(lhs)
  #   
  #   if(nrow(inds)==0){
  #     ind<-nrow(all.alias)+1
  #     all.alias[ind,]<<-NA
  #     all.alias[ind,1:4]<<-c(lhs,lhs,"",role)
  #     all.alias$all.forms[[ind]]<<-c(lhs)
  #   }
  #   
  #   return()
  # }
  
  if(lhs==""||rhs=="")
    stop("glmSEM error: null alias can not be set for '",lhs,rhs,"'")
  
  i1<-get.alias.ind(lhs)
  i2<-get.alias.ind(rhs)
  
  ind<- unique(c(i1,i2))
  
  if(length(ind)==0){
    ind<-nrow(all.alias)+1
    all.alias[ind,1:4]<<-c(lhs,rhs,"",role)
    all.alias$all.forms[[ind]]<<-unique(c(lhs,rhs))
    
  }else if(length(ind)==1){
    a<-all.alias[ind,1:2]
    if((a[1]==lhs)&(a[2]==rhs)){
      #Do nothing
      return()
    }else if(a[1]==a[2]){#Update the alias
      all.alias[ind,1:4]<<-c(lhs,rhs,"",role)
      all.alias$all.forms[[ind]]<<-unique(c(lhs,rhs,all.alias$all.forms[[ind]]))
    }else{
      if(a[2]==lhs&a[1]==rhs){
        #Ignore the new definition
        return()
      }else{
        #Multiple definions for an alias. Add them to the forms
        all.alias$all.forms[[ind]]<<-unique(c(lhs,rhs,all.alias$all.forms[[ind]]))
        # als=setdiff(all.alias$all.forms[[ind]],all.alias$lhs[ind])
        # warning("\ngmlSEM Parser Error: Multiple definition for alias:\n",
        #      a[1]," is defined as an alias for ",a[2],"\n",
        #      all.alias$lhs[ind]," is alias for '",paste0(als,collapse = ","),"'")
        return()
      }
    }
  }else if(length(ind)>1){
    
    a1<-all.alias[i1,]
    a2<-all.alias[i2,]
    
    if((a1[1]==a1[2])&(a2[1]==a2[2])){
      #remove the two rows and define the alias
      als=unique(c(lhs,rhs,all.alias$all.forms[[i1]],all.alias$all.forms[[i2]]))
      all.alias<<-all.alias[-c(i1,i2),]
      ind<-nrow(all.alias)+1
      all.alias[ind,1:4]<<-c(lhs,rhs,"",role)
      all.alias$all.forms[[ind]]<<-als
    }else{
      b1<-ifelse(a1[1]==lhs,a1[2],a1[1])
      b2<-ifelse(a2[1]==rhs,a2[2],a2[1])
      stop("\ngmlSEM Parser Error: Multiple definition for alias:\n",
           ifelse(b1!=lhs,lhs%+%" is defined as an alias for "%+%b1%+%"\n",""),
           ifelse(b2!=rhs,rhs%+%" is defined as an alias for "%+%b2%+%"\n",""),
           "Thus, ",lhs," and ",rhs," cannot be defined as alias.")
    }
  }
  return()
}

get.alias.lhs<-function(lbl){
  if(length(lbl)>1){
    for(i in seq_along(lbl))
      lbl[i]=get.alias.lhs(lbl[i])
    
    return(lbl)
  }
  
  ind<-get.alias.ind(lbl)
  
  if(length(ind)>0)
    return(all.alias[ind[1],1])
  
  # if(lbl!="1")
  #   add.alias(lbl,lbl)
  
  lbl
}

get.alias.rhs<-function(lbl){
  if(length(lbl)>1){
    for(i in seq_along(lbl))
      lbl[i]=get.alias.rhs(lbl[i])
    
    return(lbl)
  }
  
  ind<-get.alias.ind(lbl)
  
  if(length(ind)>0)
    return(all.alias[ind[1],2])
  
  # if(lbl!="1")
  #   add.alias(lbl,lbl)
  
  lbl
}
