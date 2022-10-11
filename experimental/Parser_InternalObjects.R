##Arguments
#internal machinery
data.is.provided = FALSE

#Main object returned by the parser
parsedData=list()


vars.level      <- data.frame(var = character(), level = character())

covs            <- list()          #lists of lists. one cov matrix per level per group
levels.list     <- list()          #list of levels hierarchy (levels' tree)
levels.matrix   <- matrix(0,1,1)
rownames(levels.matrix) <- colnames(levels.matrix) <- "1"

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

#Accepting models on rhs of the operators
model.operators <- names(operators.blocks)[operators.blocks%in%
                           c("measurement","regression","formative","covariance")]

operators <- names(operators.blocks)
lhs.blocks <- c("level","family","heter","group","size")

GROUP_OP <- FALSE
LEVEL_OP <- FALSE
group<-""

add.vars.to.dictionay<-function(mat, # matrix of variable names. Extracted from 'mat' attribute of the output of expand.ellipsis()
                                     # dim(mat)=n1 or dim(mat)=n1+n2, where n1=family$dim, n2=family$dim.latent
                                family=NULL,copula=NULL){
  if(!is.matrix(mat))
    mat=matrix(mat,ncol=1)
  
  levels=matrix("",nrow = nrow(mat),ncol = ncol(mat))
  
  for(i in 1:nrow(mat))
    for(j in 1:ncol(mat))
      mat[i,j]=get.alias.lhs(mat[i,j])
  
  vs=vars.level[,1]
  for(i in 1:nrow(mat))
    for(j in 1:ncol(mat))
      levels[i,j]=get.level(mat[i,j])
  

  if(!is.null(family)){
    n1=family$dim
    n2=family$dim.latent
    
    dm=dim(mat)
    
    if(dm!=n1 && dm!=(n1+n2)){
      return("Dimension mismatch")
    }
    
    tryCatch({
      for(i in 1:nrow(mat)){
        
        nm=get.alias.lhs(mat[i,])
        
        k=nrow(vars.dictionary)+1
        vars.dictionary               [k,] <<- NA
        vars.dictionary$id            [k]  <<- k
        vars.dictionary$name         [[k]] <<- nm
        vars.dictionary$dim           [k]  <<- n1
        vars.dictionary$dim.latent    [k]  <<- n2
        vars.dictionary$support      [[k]] <<- family$getSupport(family)
        vars.dictionary$is.phantom    [k]  <<- FALSE
        vars.dictionary$family       [[k]] <<- family
        vars.dictionary$hyper.params [[k]] <<- family$runtime.pars.data.frame
        
        if(dm>n1){
          supps=family$getSupport.latent(family)
          fams=family$getFamily.latent(family)
          
          for(j in n1:(n1+n2)){
            k=nrow(vars.dictionary)+1
            nm=get.alias.lhs(mat[i,j])
            
            vars.dictionary               [k,] <<- NA
            vars.dictionary$id            [k]  <<- k
            vars.dictionary$name         [[k]] <<- nm
            vars.dictionary$dim           [k]  <<- 1
            vars.dictionary$dim.latent    [k]  <<- 0
            vars.dictionary$support      [[k]] <<- list(supps[[j-n1+1]])
            vars.dictionary$is.phantom    [k]  <<- TRUE
            vars.dictionary$latent        [k]  <<- TRUE
            vars.dictionary$phantom.for   [k]  <<- k-1
            vars.dictionary$family       [[k]] <<- list(fams[[j-n1+1]])
            vars.dictionary$hyper.params [[k]] <<- fams[[j-n1+1]]$runtime.pars.data.frame
          }
        }
      }
      
    },error=function(e){
      return(paste0("Parser internal error: ",e))
    })
  }
  
  if(!is.null(copula)){
    #Add copula structure to the covs.list
    #covs.list is a named list of covariances. The name is the levels'lhs name, and the
    # value is two elements. The first element
    #covs.list[[lev.name]]=list(copulas,covMat)
    levs=get.level(mat)
    levsc=c(levs)
    if(!all.equal(levsc,rep(levsc[1],length(levsc))))
      return("Parser internal error: Variables belong to different leves. You can assign copula functions to induce conditional dependency for variables vary at the same level.")
    
    lev=levc[1]
    covs.list[[lev]]
  }
    
  return(TRUE)
}

get.level.lhs<-function(lhs){
  lhs=gsub("<~|~|=~","",lhs,perl = TRUE)
  lhs=gsub("^.*:","",lhs,perl = TRUE)
  lhs=trim(lhs)
  lhs=expand.ellipsis(lhs)
  mat=attr(lha,"mat")
  get.level(mat)
}

get.level <- function(mat){
  if(is.matrix(mat)){
    levs=matrix("",nrow(mat),ncol(mat))
    
    for(i in 1:nrow(mat))  
      for(j in 1:ncol(mat))
        levs[i,j]=get.level(mat[i,j])
    return(levs)
  }
  
  if(length(mat)>1){
    levs=rep("",length(mat))
    for(i in seq_along(var))  
      levs[i]=get.level(mat[i])
    return(levs)
  }
  
  if(nrow(vars.level)==0)
    return(NA)
  var=get.alias.lhs(mat)
  ind=which(vars.level[,1]==var)
  if(length(ind)==0)
    return(NA)
  return(vars.level[ind[1],2])
}

add.vars.to.level<-function(mat,lev){
  if(is.matrix(mat)){
    for(i in 1:nrow(mat))  
      for(j in 1:ncol(mat))
        add.vars.to.level(mat[i,j],lev)
    return()
  }
  
  if(length(mat)>1){
    for(i in seq_along(var))  
      add.vars.to.level(mat[i],lev)
    return()
  }
  
  if(nrow(vars.level)==0)
    return(NA)
  
  var=get.alias.lhs(mat)
  lev=get.alias.lhs(lev)
  
  ind=which(vars.level[,1]==var)
  
  if(length(ind)>0){
    if(vars.level[ind,2]==lev)
      return()
    
    stop("gmlSEM error: double level specification for variable '",mat,"'")
  }
  
  n=nrow(vars.level)+1
  vars.level[n,] <<- c(var,lev)
}



#gmlSEM parser environment to resolve parameter values
env=new.env()

#Parser Output
# for each family with a latent generating process, two records are added pointing to each other. One for the observed variable(s), and one for the latent phantom variable(s)
# The phantom record will not be included if it is not named in the family: block
vars.dictionary <- data.frame(
  id           = integer()   ,
  name         = I(list())   , #the name of variable (or vector of variables) 
  dim          = integer()   ,
  dim.latent   = integer()   ,
  support      = I(list())   , # Support of random variable(s)
  is.phantom   = logical()   , #Phantom variables are simply latent variables defines the generating process of a family of distributions
  latent       = logical()   , #NA if needs.data.scan. TRUE for phantom variables.
  phantom.for  = integer()   , #if the variable (or vector of variables) are phantom, to which row they belongs to
  family       = I(list())   ,
  hyper.params = I(list())   , # A named list of hyper parameters e.g. list(argname=varname) which the family of distribution depends on e.g. list("trials"="N") for bin(trials=N)
                               # Hyper parameters must the be resolved either in the data frame or the syntax environment
  appeared.as.fac = logical() , #Appeared as a factor in a measurement model
  appeared.as.ran = logical() , #Appeared as a random effect in a regression model
  appeared.as.res = logical() , #Appeared as a response in a regression model
  appeared.as.res.h = logical() , #Appeared as a response in a het. reg. model
  
  vary         = character() ,  # The level at which the variable(s) vary. NA if needs more speculation
  reg.model    = character() ,
  heter.model  = I(list())   
)


all.alias       <- data.frame(
  lhs=character(),
  rhs=character(),
  var.name=character(),
  role=character(), #PAR, LV, OV
  all.forms=I(list())
)


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
  
  #update level names
  if(nrow(levels.matrix)>0){
    colnames(levels.matrix)<<-rownames(levels.matrix)<<-get.alias.lhs(rownames(levels.matrix))
  }
  
  if(length(lev)>1){
    for(i in seq_along(lev))
      add.levels(lev[i])
    return()
  }
  
  if(lev=="")
    return()
  
  lev=trim(lev)
  
  if(grepl("(",lev,fixed = TRUE)){
    levv=captured.groups(lev," *(?'var'.*) *\\( *(?'varn'.*) *\\)")
    lev=trim(levv$var)
    lev.par=trim(levv$varn)
    add.alias(lev.par,lev,"OV")
  }
  
  lev=get.alias.lhs(lev)
  
  # if(lev=="1")
  #   return()  #Do not add the base level to the matrix
  
  rn=get.alias.lhs(rownames(levels.matrix))
  if(lev %in% rn)
    return()
  
  #Rectify levels.matrix and vars.level based on new aliases
  levels.matrix<<-rbind(levels.matrix,0)
  levels.matrix<<-cbind(levels.matrix,0)
  colnames(levels.matrix)<<-rownames(levels.matrix)<<-c(rn,lev)
  
}

set.levels.matrix<-function(lhs,rhs){
  
  if(length(lhs)>1 || length(rhs)>1){
    for(i in 1:length(lhs))
      for(j in 1:length(rhs))
        set.levels.matrix(lhs[i],rhs[j])
    return()
  }
  
  if(lhs==""||rhs=="")
    return()
  lhs=get.alias.lhs(lhs)
  if(lhs=="1")
    return()  #ignore the base level
  rhs=get.alias.lhs(rhs)
  if(rhs=="1")
    stop("gmlSEM error: base level can not be superior to other levels.")
  
  add.levels(lhs)
  add.levels(rhs)
  levels.matrix[lhs,rhs]<<-1
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
  
  #Removing () and replacing == with single =
  lhs=trim(gsub("(","",gsub(")","",gsub('"',"",gsub("'","",gsub("={2,}","=",lhs,perl = TRUE),fixed = TRUE),fixed = TRUE),fixed = TRUE),fixed = TRUE))
  rhs=trim(gsub("(","",gsub(")","",gsub('"',"",gsub("'","",gsub("={2,}","=",rhs,perl = TRUE),fixed = TRUE),fixed = TRUE),fixed = TRUE),fixed = TRUE))

 if(rhs=="1"){
   #Swap the sides
   rhs=lhs
   lhs="1"
 }
  
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
  
  #rectify vars.level
  if(nrow(vars.level)>0){
    vars.level[,1] <<- get.alias.lhs(vars.level[,1])
    vars.level[,2] <<- get.alias.lhs(vars.level[,2])
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

merge.models<-function(str){
  if(length(str)==1)
    return(str)
  
  str=paste0(str,collapse = "+")
  str=gsub("\\+*\\-+\\+*","-",str,perl = TRUE)
  str=gsub("\\+*\\*+\\+*","*",str,perl = TRUE)
  str=gsub("\\+{2,}","+",str,perl = TRUE)
  str
}

expand.dots.in.syntax<-function(model.syntax){
  pattern  = "(?J)(?>(?>(?>([^:,\\+\\n\\r\\<\\>~\\-\\(\\)|]|(?'cond'\\|))*(?>\\((?>(?!\\.\\.\\.)([^:,\\+\\n\\r\\<\\>~\\-\\(\\)|]|\\k'cond')|(?1))*\\))?){1,2})(?'sep'[,+]))?(?>(?>(?!\\.\\.\\.)([^:,\\+\\n\\r\\<\\>~\\-\\(\\)|]|\\k'cond')*(?:\\((?>[^:,\\+\\n\\r\\<\\>~\\-()|]|(?1))*\\))?){1,2})(?'sep'[,+])\\.\\.\\.\\k'sep'(?>(?>(?!\\.\\.\\.)[^:,\\+\\n\\r\\<\\>~\\-\\(\\)|]*(?>\\((?>(?!\\.\\.\\.)([^:,\\+\\n\\r\\<\\>~\\-\\(\\)|]|\\k'cond')|(?1))*\\))?){1,2})"
  capture.ellipsis=gregexpr(pattern, model.syntax, perl=TRUE,ignore.case =TRUE)
  match.start=capture.ellipsis[[1]]
  match.length=attr(capture.ellipsis[[1]],"match.length")
  for(i in rev(seq_along(match.start))){
    txt.ellipsis=substring(model.syntax,match.start[i],
                            match.start[i]+match.length[i]-1)
    txt.expanded=expand.ellipsis(txt.ellipsis)
    model.syntax=paste0(substring(model.syntax,1,match.start[i]-1),
                         txt.expanded,
                         substring(model.syntax,match.start[i]+match.length[i],nchar(model.syntax)))
  }
  model.syntax
}

extract.modelterms.lhs<-function(lhs){
  terms.match=gregexpr("^(?: *(?'blockname'[\\w\\._]+) *:)? *(?'tterms'[^~<>]+) *(?'op'<~|=~|~) *", lhs,perl = TRUE)[[1]]
  if(terms.match[1]<0)
    stop("\ngmlSEM error: left hand side of model mismatch:\n",lhs)
  
  c.s=attr(terms.match,"capture.start")
  c.l=attr(terms.match,"capture.length")
  
  block.name=""
  if(c.l[1]>0)
    block.name = trim(substring(lhs,c.s[1],c.s[1]+c.l[1]-1))
  yterms = trim(substring(lhs,c.s[2],c.s[2]+c.l[2]-1))
  op     = substring(lhs,c.s[3],c.s[3]+c.l[3]-1)
  
  mat=expand.ellipsis(yterms)
  mat=attr(mat,"mat")
  
  list(yterms=yterms,block.name=block.name,op=op,mat=mat)
}


extract.modelterms.rhs0<-function(tterm){
  re=data.frame(mod    = character(),     #term := mod*varfun(var)@label
                var    = character(),
                varfun = character(), 
                label  = character())
  terms.match=gregexpr("(?J)(?:(?'mod'\\b[\\w\\._]+\\b(?: *\\([\\w\\._=\"']+\\))?) *\\* *)?(?:(?:(?'varfun'[\\w\\._]+) *\\( *(?'var'[\\w\\._]+) *\\))|(?'var'[\\w\\._]+))(?:@(?'lat'[\\w\\._]+))?", tterm,perl = TRUE,ignore.case =TRUE)[[1]]
  c.s=attr(terms.match,"capture.start")
  c.l=attr(terms.match,"capture.length")    
  if(terms.match[1]>=0)
    for(j in 1:nrow(c.s)){
      mod=trim(substr(tterm,c.s[j,"mod"],c.s[j,"mod"]+c.l[j,"mod"]-1))
      varfun=trim(substr(tterm,c.s[j,"varfun"],c.s[j,"varfun"]+c.l[j,"varfun"]-1))
      s.var=ifelse(c.s[j,3]>0,c.s[j,3],c.s[j,4])
      l.var=ifelse(c.s[j,3]>0,c.l[j,3],c.l[j,4])
      var=trim(substr(tterm,s.var,s.var+l.var-1))
      lab=trim(substr(tterm,c.s[j,"lat"],c.s[j,"lat"]+c.l[j,"lat"]-1))
      
      n=nrow(re)+1
      re[n,]<-c(mod,var,varfun,lab,lag)
    }
  
  re
}

extract.random.effects<-function(model.syntax,lhs=""){
  
  lhs=gsub("~","",gsub("=~","",lhs,fixed = TRUE),fixed = TRUE)
  if(grepl("=",lhs,fixed = TRUE)){
    #in the case that lhs is a dummy variable e.g. y=0
    lhs=get.alias.rhs(lhs)
  }
  re=data.frame(mod    = character(),
                var    = character(),    #mod*varfun(var)@label|leve
                varfun = character(), 
                label  = character(),
                level  = character(),
                dummy.label=character())
  
  capture.within <- gregexpr("\\( *(?'term'.*) *\\| *(?'level'[\\w\\._]+) *\\)", model.syntax,perl = TRUE,ignore.case =TRUE)
  match.start=capture.within[[1]]
  match.length=attr(capture.within[[1]],"match.length")
  capture.length=attr(capture.within[[1]],"capture.length")
  capture.start=attr(capture.within[[1]],"capture.start")
  if(match.start[1]>=0){
    for(i in 1:nrow(capture.start)){
      
      tterm=substring(model.syntax,capture.start[i,1],capture.start[i,1]+capture.length[i,1]-1)
      tlev =trim(substring(model.syntax,capture.start[i,2],capture.start[i,2]+capture.length[i,2]-1))
      
      if(!grepl("^[\\w\\._]+$",tlev,perl = TRUE)){
        stop("\ngmlSEM error: levelname mismatch at:\n",
             substring(model.syntax,match.start[i],match.start[i]+match.length[i]-1))
      }
      
      terms=extract.modelterms.rhs(tterm)
      if(nrow(terms)>0)
        for(j in 1:nrow(terms)){
          var=terms$var[j]
          dummylab=""
          if(lhs!="")
            dummylab=paste0(lhs,".",tlev,".",var)
          n=nrow(re)+1
          re[n,]<-c(terms[j,],tlev,dummylab)
        }
    }
  }
  
  re
}


levels.reachUp<-function(levs,lev){
  #go deep to the level lev and return all of its descendant branches
  if(length(levs)==0)
    return(list())
  a=list()
  for(i in 1:length(levs)){
    if(names(levs)[i]==lev){
      return(levs[[i]])
    }else{
      a=reachDown(levs[[i]],lev)
      if(length(a)>0)
        return(a)
    }
  }
  a
}

levels.goDown<-function(levs,lev){
  #Go down in levs until it reach to lev
  #If ound the lev, sen a signal back to the top
  #Otherwise it will echo a zero to the top
  if(length(levs)==0)
    return(0)
  if(lev %in% names(levs))
    return(1)
  
  for(l in levs){
    res=goDown(l,lev)
    if(res!=0)
      return(res)
  }
  return(0)
}

levels.are.consistent<-function(down,up){
  #Return
  #TRUE: if 'down' is a descendent to 'up' in the levels' tree
  #FALSE: otherwise
  
  down=get.alias.lhs(down)
  if(down=="1")
    return(TRUE)
  
  up=get.alias.lhs(up)
  if(up=="1")
    return(FALSE)
  
  if(down==up)
    return(TRUE)
  
  if(! down %in% rownames(levels.matrix))
    stop("\ngmlSEM error: undefined level",down)
  
  if(! up %in% rownames(levels.matrix))
    stop("\ngmlSEM error: undefined level: ",up)
  
  l.up=levels.reachUp(levels,up)
  res=levels.goDown(l.up,down)
  
  if(res==1)
    return(TRUE)
  
  return(FALSE)
}


# Look for any inconsistency in level structure
# Levels are consistent if and only if 
#  any submatrix of levels.matrix contains a zero row
validate.level.hierarchy<-function(){
  n.levels<-nrow(levels.matrix)
  for(i in 1:n.levels){
    cmbn<-combn(1:n.levels,i)
    for(j in 1:ncol(cmbn)){
      .subm <- levels.matrix[cmbn[,j],cmbn[,j]]
      if(i==1){
        if(.subm==1)
          stop(paste0("gmlSEM error in levels: ",rownames(.subm)))
      }else{
        any_zero<-any(sapply(1:nrow(.subm), function(k){sum(.subm[k,])==0}))
        if(!any_zero)
          stop("gmlSEM error: Levels are inconsistent. There is a loop in the level structure")
      }
    }
  }
  
  # Constructing Level structure
  levels.list<<-list()
  .zeros<-sapply(1:nrow(levels.matrix), function(k){sum(levels.matrix[k,])==0})
  .levels<- rownames(levels.matrix)[which(.zeros)]
  for(lev in .levels){
    levels.list[[lev]]<<-getSubLevels(lev)
  }
  
}