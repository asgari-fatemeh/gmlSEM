space.sub= "_______"
tab.sub=   "-------"

"%+%"<-function(a,b){paste0(a,b)}

##The decision depends on the provided data
decide.latent<-function(xname){
  ind<- which(apply(vars.dictionary,1,function(m){xname%in%m$alias}))
  if(length(ind)==0)
    stop("gmlSEM error: variable '",xname,"' not found!")
  m=vars.dictionary[ind[1],]
  if(is.na(m$latent)){ #needs.scan.data
    if(m$appeared.as.fa || m$appeared.as.re)
      return(TRUE)
    
    if(data.is.provided){
      cnames<-colnames(data)
      ind=which(cnames%in%alias)
      if(length(ind)>1)
        stop("glmSEM error: Ambiguty in data. '",
             paste0(alias[ind],collapse = "','"),
             "' are defined as alias while all exists in the data as separte columns!")
      
      length(ind)==0  #Consider it as a latent if it does not exist in the data
    }else{
      NA
    }
  }else{ 
    #Stick to previous decision
    m$latent 
  }
}

resolve<-function(xname=NA,i.row=NA){
  if(!is.na(i.row)){
    m=vars.dictionary[ind[1],]
    alias=m$alias
    latent=m$latent
  }else{
    alias=xname
    ind<- which(apply(vars.dictionary,1,function(m){xname%in%m$alias}))
    if(length(ind)==0)
      stop("gmlSEM error: variable '",xname,"' not found!")
    m=vars.dictionary[ind[1],]
    alias=m$alias
    latent=m$latent
  }
  #if is latent return its support
  if((is.na(latent) || !latent) && !data.is.provided)
    stop("gmlSEM error: to resolve the variable '",paste0(alias,collapse = "~"),' data needs to be provided! ')
  
  if(is.na(latent)){
    #Scan data
    cnames<-colnames(data)
    ind=which(cnames%in%alias)
    if(length(ind)>1)
      stop("glmSEM error: Ambiguty in data. '",
           paste0(alias[ind],collapse = "','"),
           "' are defined as alias while all exists in the data as separte columns!")
    if(length(ind)==0){
      if(is.null(m$family) || is.na(m$family))
        stop("gmlSEM error: No distribution family is assigned to the random variable '",
             m$name,"'")
      return(m$family$support(m$family))
    }
    return(data[,ind])
  }else if(latent){
    if(is.null(m$family) || is.na(m$family))
      stop("gmlSEM error: No distribution family is assigned to the random variable '",
           m$name,"'")
    return(m$family$support(m$family))
  }else{
    if(!data.is.provided)
      stop("gmlSEM error: to resolve the variable '",paste0(alias,collapse = "~"),' data needs to be provided! ')
    cnames<-colnames(data)
    ind=which(cnames%in%alias)
    aa<-setdiff(m$alias,m$name)
    if(length(ind)==0)
      stop("gmlSEM error: Failed to find variable '",m$name,"' ",
           ifelse(length(aa)==0,"",
                  paste0("or either its aliases: '",
                         paste0(aa,collapse="','"),"', ")),
           " in the data")
    return(data[,ind])
  }
}

switch.space<-function(x){
  ifelse(grepl(" ",x,fixed = TRUE),
         gsub(" ",space.sub,gsub("\\t",tab.sub,x,fixed = TRUE),fixed=TRUE),
         gsub(space.sub," ",gsub(tab.sub,"\\t",x,fixed = TRUE),fixed=TRUE)
  )
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

split.vars<-function(x){
  strsplit(txt,"[,\\+]",perl = TRUE)[[1]]
}

