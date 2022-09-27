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

is.valid.varname<-function(x){
  sapply(x, function(xx)grepl("^[\\w\\.][\\w\\._\\-]*",xx,perl = TRUE))
}

captured.groups<-function(txt,pat,...){
  gpr=gregexpr(pat,txt,perl=TRUE,...)[[1]]
  capture.start<-attr(gpr,"capture.start")
  capture.length<-attr(gpr,"capture.length")
  
  groups<-colnames(capture.start)
  gr<-list()
  for(i in seq_along(groups)){
    if(groups[i]=="" | capture.length[1,i]==0)
      next
    gr[[groups[i]]]<-substr(txt,capture.start[1,i],
                            capture.start[1,i]+capture.length[1,i]-1)
  }
  
  gr 
}


captured.groups.list<-function(txt,gpr){
  capture.start<-attr(gpr,"capture.start")
  capture.length<-attr(gpr,"capture.length")
  
  groups<-colnames(capture.start)
  gr<-list()
  for(i in seq_along(groups)){
    if(groups[i]=="" | capture.length[1,i]==0)
      next
    gr[[groups[i]]]<-substr(txt,capture.start[1,i],
                            capture.start[1,i]+capture.length[1,i]-1)
  }
  
  gr 
}

captured.groups.dataframe<-function(txt,gpr){

  capture.start<-attr(gpr,"capture.start")
  capture.length<-attr(gpr,"capture.length")
  
  groups<-colnames(capture.start)
  
  ks=sum(groups!="")
  gr=matrix("",nrow = length(gpr),ncol=ks)
  colnames(gr)=groups[groups!=""]
  for(j in 1:length(gpr)){
  k=0
  vals=c()
  for(i in seq_along(groups)){
    if(groups[i]=="")
      next
     
    k=k+1
    if(capture.length[j,i]>0)
      gr[j,k]=substr(txt,capture.start[j,i],
                              capture.start[j,i]+capture.length[j,i]-1)
   }
  }
  gr 
}


perlsplit<-function(x,pat,drop.captured.sgroups.only=TRUE){
  #expect just one match
  if(!grepl(pat,x,perl=TRUE))
    return(x)
  
  if(!drop.captured.sgroups.only)
    return(strsplit(x,pat,perl= TRUE)[[1]])
  
  text=gregexpr(pat,x,perl= TRUE)[[1]]
  
  cap.start=attr(text,"capture.start")
  cap.len=attr(text,"capture.length")
  
  emptInds=rowSums(cap.start)==0 | rowSums(cap.len)==0
  
  cap.start =  cap.start[!emptInds,]
  cap.len   =    cap.len[!emptInds,]
  
  if(is.null(dim(cap.start))){
    cap.start=matrix(cap.start,nrow=1)
    cap.len=matrix(cap.len,nrow=1)
  }
  
  if(nrow(cap.start)>1)
    return(stopp())
  
  txtt=character()
  st=1
  j=1
  for(i in 1:ncol(cap.start)){
    if(cap.len[i]==0)
      next
    txtt[j]=substr(x,st,cap.start[i]-1)
    st=cap.start[i]+cap.len[i]
    j=j+1
  }
  if(st<=nchar(x))
    txtt[length(txtt)+1]=substr(x,st,nchar(x))
  
  txtt
}


sep.par<-function(x){
  capture.level <- gregexpr("(\\b[\\w\\.]+\\b)(?:\\(([\\w\\.]+)\\))?",x, perl=TRUE)
  capture.start<-attr(capture.level[[1]],"capture.start")
  capture.length<-attr(capture.level[[1]],"capture.length")
  rhs1<-substr(rhs,capture.start[i2,1],capture.start[i2,1]+capture.length[i2,1])
  rhs2<-""
  if(capture.start[i2,2]>0)
    rhs2<-substr(rhs,capture.start[i2,2],capture.start[i2,2]+capture.length[i2,2])
  
  list(var.out=rhs1,var.in=rhs2)
}

any.intersect<-function(list.of.seqs){
  inds=sapply(list.of.seqs, function(x){is.na(x) || is.null(x)})
  inds=which(!inds)
  if(length(inds)<=1)
    return(FALSE)
  seq=list.of.seqs[[inds[1]]]
  for(i in 2:length(inds)){
    if(length(intersect(seq,list.of.seqs[[inds[i]]]))>0)
      return(FALSE)
    seq=union(seq,list.of.seqs[[inds[1]]]  )
  }
  return(TRUE)
}

getSubLevels<-function(level){
  levels<-list()
  .subm<-levels.matrix[,level]
  if(sum(.subm)==0)
    return(levels)
  sub.levels<-rownames(levels.matrix)[which(.subm==1)]
  
  for(lev in sub.levels){
    levels[[lev]]<-getSubLevels(lev)
  }
  levels
}

