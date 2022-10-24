space.sub= "_____----_____"
tab.sub=   "_----____----_"

"%+%"<-function(a,b){paste0(a,b)}

trim <- function (x,omit.quotes=FALSE){
  if(omit.quotes)
    x=gsub("\"|'", "", x)
  
  x=gsub("^\\s+|\\s+$", "", x)
  x
} 

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

switch.space.within.quotes<-function(txt,seal.spaces=NULL,restore.spaces=NULL){
  if(nchar(txt)==0)
    return(txt)
  
  if(is.null(seal.spaces))
    seal.spaces=FALSE
  if(is.null(restore.spaces))
    restore.spaces=FALSE
  
  if(!seal.spaces & !restore.spaces){
    if(grepl(" ",txt,fixed = TRUE)|grepl("\t",txt,fixed = TRUE)){
      seal.spaces=TRUE
    }else if(grepl(space.sub,txt,fixed = TRUE) | grepl(tab.sub,txt,fixed = TRUE)){
      restore.spaces=TRUE
    }else{
      return(txt)
    }
  }
  
  spc=gregexpr("\"[^\"]*\"|'[^']*'",txt,perl = TRUE)[[1]]
  if(spc[1]>=0){
    m.s=spc
    m.l=attr(spc,"match.length")
    for(i in rev(seq_along(spc))){
      txt=paste0(substr(txt,1,m.s[i]),
                   switch.space(substr(txt,m.s[i]+1,m.s[i]+m.l[i]-2),seal.spaces = seal.spaces, restore.spaces = restore.spaces ),
                   substr(txt,m.s[i]+m.l[i]-1,nchar(txt)))
    }
  }
  
  txt
}

switch.space<-function(x,seal.spaces=NULL,restore.spaces=NULL){
  if(nchar(x)==0)
    return(x)
  
  if(is.null(seal.spaces))
    seal.spaces=FALSE
  if(is.null(restore.spaces))
    restore.spaces=FALSE
  
  if(!seal.spaces & !restore.spaces){
    if(grepl(" ",x,fixed = TRUE)|grepl("\t",x,fixed = TRUE)){
      seal.spaces=TRUE
    }else if(grepl(space.sub,x,fixed = TRUE) | grepl(tab.sub,x,fixed = TRUE)){
      restore.spaces=TRUE
    }else{
      return(x)
    }
  }
  
  if(seal.spaces){
    gsub(" ",space.sub,gsub("\\t",tab.sub,x,fixed = TRUE),fixed=TRUE)
  }else if(restore.spaces){
    gsub(space.sub," ",gsub(tab.sub,"\\t",x,fixed = TRUE),fixed=TRUE)
  }
}



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


#get arglist of scaling functions in model syntax i.e. modifs*vardesc@lat where vardesc:=scalefun(xvar(xlag),funargs)
#argslist is the listified of funargs in previous sentence, which is ready to pass into do.call
getArgList<-function(csv,forceNumeric=FALSE,forceCharactersAsNames=FALSE){
  if(nchar(csv)==0)
    return(list())
  
  arglist=list()
  nameList=c()
  
  if(substr(csv,1,1)!=",")
    csv=paste0(",",csv)
  coms=gregexpr(fun.args.format,csv,perl=TRUE)[[1]]
  cm.s=attr(coms,"capture.start")
  cm.l=attr(coms,"capture.length")
  
  for(k in 1:length(coms)){
    funpar=substr(csv,cm.s[k,"funpar"],cm.s[k,"funpar"]+cm.l[k,"funpar"]-1)
    funarg=substr(csv,cm.s[k,"funarg"],cm.s[k,"funarg"]+cm.l[k,"funarg"]-1)
    funarg=str2lang(funarg)
    
    if(is.call(funarg)){
      #Get dependencies
      funargNames=getNameArgs(funarg)
      nameList=unique(c(nameList,funargNames))
    }
    
    if(forceNumeric){
      if(is.call(funarg) && length(funargNames)==0)
        funarg=eval(funarg)
      if(!is.call(funarg) & !is.name(funarg)){
        funargN=suppressWarnings(as.numeric(funarg))
        if(!is.na(funargN)){
          funarg=funargN
        }else if(as.character(funarg)=="NA"){
          funarg=NA
        }
      }
    }
    
    if(forceCharactersAsNames && is.character(funarg)){
      funarg=as.name(funarg)
    }
    
    if(is.name(funarg))
      nameList=unique(c(nameList,funarg))
    
    if(funpar!=""){
      arglist[[funpar]]=funarg
    }else{
      arglist[[length(arglist)+1]]=funarg
    }
  }
  
  list(argList=arglist, nameList=nameList)
}

getNameArgs<-function(strcall){
  anames=c()
  if(length(strcall)<2)
    return(c())
  for(i in 2:length(strcall)){
    if(is.name(strcall[[i]])){
      anames=c(anames,strcall[[i]])
    }else if(is.call(strcall[[i]])){
      anames=c(anames,getNameArgs(strcall[[i]]))
    }
  }
  if(is.null(anames))
    anames=list()
  anames
}


getModifList<-function(modifs){
  #list(start=list(),value=list(), depends.on=list()) #list of values or quotes
  res=list(start=list(),
           value=list(),
           depends.on=c())
  
  if(length(modifs)==0)
    return(res)
  
  res0=list()
  
  for(k in 1:length(modifs)){
    res0[[k]]=list(start=list(),
                   value=list(),
                   depends.on=c())
    
    modif=modifs[k]
    if(substring(modif,nchar(modif))=="*")
      modif=substr(modif,1,nchar(modif)-1)
    #decompose modif against pat=/single.modifier.format/
    gmodif=gregexpr(single.modifier.format,modif,perl=TRUE,ignore.case = TRUE)[[1]]
    cmo.s=attr(gmodif,"capture.start")
    cmo.l=attr(gmodif,"capture.length")
    
    funname=substr(modif,cmo.s[1,"funname"],cmo.s[1,"funname"]+cmo.l[1,"funname"]-1)
    funargs=substr(modif,cmo.s[1,"funargs"],cmo.s[1,"funargs"]+cmo.l[1,"funargs"]-1)
    singlemod=substr(modif,cmo.s[1,"singlemod"],cmo.s[1,"singlemod"]+cmo.l[1,"singlemod"]-1)
    
    if(singlemod!=""){
      arglist=getArgList(singlemod,forceNumeric=TRUE,forceCharactersAsNames=TRUE)
      res0[[k]]$value=arglist$argList
      res0[[k]]$depends.on=unique(c(res0[[k]]$depends.on,arglist$nameList))
    }else{
      arglist=getArgList(funargs,forceNumeric=TRUE,forceCharactersAsNames=TRUE)
      res0[[k]]$depends.on=unique(c(res0[[k]]$depends.on,arglist$nameList))
      
      if(tolower(funname) %in% c("","c","label")){
        #Group definiotn of labels or fixed values
        res0[[k]]$value=arglist$argList
      }else{
        #Starting point
        res0[[k]]$start=arglist$argList
      }
    }
    if(is.null(res0[[k]]$depends.on))
      res0[[k]]$depends.on=c()
  }
  res0
}

combineModifiers<-function(modifs,vardesc=""){
  #modifiersList=list of list(start=list(),value=list(), depends.on=list())
  if(length(modifs)<2)
    return(modifs)
  
  inds=which(sapply(seq_along(modifs), function(mi)is.list(modifs[[i]])))
  
  modifs2=list()
  for(i in inds)
    modifs2[[i]]=modifs[[inds[i]]]
  
  if(length(modifs2)<2)
    return(modifs2)
  
  modifs=modifs2
  
  starts=lapply(modifs, function(m)m$start)
  values=lapply(modifs, function(m)m$value)
  depends.on=unique(unlist(lapply(modifs, function(m)m$depends.on)))
  
  if(is.null(depends.on))
    depends.on=c()
  
  startsp=sapply(starts,function(m)length(m)>0)
  valuesp=sapply(values,function(m)length(m)>0)
  
  err=""
  
  if(sum(startsp)>1)
    err="multiple start modifiers"
  
  if(sum(valuesp)>1)
    err="multiple coefficient modifiers"
  
  if(err!=""){
      gmlSEMerror(paste0(modifs," for ",vardesc))
  }
  
  start=value=list()
  if(sum(startsp)==1)
    start=modifs[[which(startsp)]]$start
  
  if(sum(valuesp)==1)
    value=modifs[[which(valuesp)]]$value
  
  list(start=start,value=value,depends.on=depends.on)
  
}


getPartialModel.for.singleEffect<-function(single.effect,
                                           is.random=FALSE,
                                           level=NA,
                                           level.lag=0){
  
  partial.model.data = new.effects.data.frame()
  
  i=1
  
  partial.model.data[i,] = NA
  partial.model.data$is.random[i] = is.random
  partial.model.data$level[i]     = level
  
  if(!is.na(level))
    partial.model.data$level.lag[i]     = level.lag
  
  if(grepl("(?sm)^\\(*\\-1|0\\)*$",single.effect,perl=TRUE)){
    #Exclude intercept
    partial.model.data$vardesc[i]   = "-1"
    partial.model.data$varName[i]   = "-1"
    partial.model.data$modified[i]  = FALSE
    partial.model.data$scaled[i]    = FALSE
    partial.model.data$var.lag[i]   = 0
  }else{
    m.s=group.terms=gregexpr(singe.effect.elements.format,single.effect,perl = TRUE,ignore.case = TRUE)[[1]]
    m.l=attr(group.terms,"match.length")
    c.s=attr(group.terms,"capture.start")
    c.l=attr(group.terms,"capture.length")
    modifiers=vardesc=varlat=""
    vardesc=substr(single.effect,c.s[1,"elemvar"],c.s[1,"elemvar"]+c.l[1,"elemvar"]-1)
    
    if(grepl("(?sm)^\\(*\\+?1\\)*$",single.effect,perl=TRUE)){
      #Include intercept
      partial.model.data$vardesc[i]     = vardesc
      partial.model.data$varName[i]     = vardesc
      partial.model.data$var.lag[i]     = 0
      partial.model.data$scaled[i]      = FALSE
      
    }else{
     
      varelems=gregexpr(vardesc.format,vardesc,perl=TRUE)[[1]]
      varelems.cap=captured.groups.list(vardesc,varelems)
      
      partial.model.data$vardesc[i]     = trim(vardesc,omit.quotes=TRUE) #for identifying purpose to merge the common terms later
      partial.model.data$varName[i]     = varelems.cap$varname
      partial.model.data$var.lag[i]     = ifelse(is.null(varelems.cap$varlag),0,varelems.cap$varlag)
      partial.model.data$scaled[i]      = !is.null(varelems.cap$varfun)
      
      if(partial.model.data$scaled[i]){
        partial.model.data$scalefun[i]  = varelems.cap$varfun
        partial.model.data$scalefun.args[[i]] = list()
        partial.model.data$sacelfun.depends.on[[i]] = list()
        
        if(!is.null(varelems.cap$varfunargs)){
          arglist=getArgList(varelems.cap$varfunargs,forceNumeric=TRUE,forceCharactersAsNames=TRUE)
          partial.model.data$scalefun.args[[i]]=arglist$argList
          partial.model.data$sacelfun.depends.on[[i]]=arglist$nameList
        }
      }
       
    }
    
    partial.model.data$modified[i]    = FALSE
    
    if(c.s[1,"elemmodifs"]>0)
      modifiers=substr(single.effect,c.s[1,"elemmodifs"],c.s[1,"elemmodifs"]+c.l[1,"elemmodifs"]-1)
    if(c.s[1,"elemlat"]>0)
      varlat=substr(single.effect,c.s[1,"elemlat"],c.s[1,"elemlat"]+c.l[1,"elemlat"]-1)
    
    if(!is.random & varlat!=""){
      gmlSEMerror(paste0("assigning a name to a non-random term is not allowed in ",single.effect))
    }
    if(varlat!=""){
      partial.model.data$latName[i]  = varlat  
    }
    
    if(modifiers!=""){
      #decompose modifiers
      modifs=gregexpr(modifer.split.format,modifiers,perl = TRUE,ignore.case = TRUE)[[1]]
      mo.l=attr(modifs,"match.length")
      
      modif=c()
      for(k in seq_along(modifs)){
        modif[k]=substr(modifiers,modifs[k],modifs[k]+mo.l[k]-1)
      }
      modifs=getModifList(modif)
      modifs=combineModifiers(modifs,vardesc)
      
      partial.model.data$modified[i]=TRUE
      partial.model.data$modifiers[[i]]=modifs
      partial.model.data$modifiers.depends.on[[i]]=if(is.null(modifs$depends.on)){list()}else{modifs$depends.on}
    }
    
  }
  
  partial.model.data$depends.on[[i]]=unique(c(
    partial.model.data$sacelfun.depends.on[[i]],
    partial.model.data$modifiers.depends.on[[i]]))
  
  partial.model.data
}



#The signature for model's effect term data
new.effects.data.frame<-function(){
  data.frame(
    vardesc   = character() , # raw var desc for identifying purpose
    is.random = logical()   ,
    level     = character() ,
    level.lag = integer()   ,
    varName   = character() ,
    latName   = character() ,  #If the term is a random effect, it can be named
    var.lag = integer() , # zero for non-lagged
    scaled    = logical()   ,
    scalefun  = character() , #list(fun.name,fun.args,fun.args.depends.on)
    scalefun.args = I(list()) ,
    sacelfun.depends.on  = I(list()) ,   #symbols which scalefun depends on
    modified  = logical()   ,
    modifiers = I(list())   , #list(label=list(),start=list(),value=list())
    modifiers.depends.on = I(list()) , #symbols which modifiers depends on
    depends.on           = I(list())  #All symbols that the term depends on
  )
}

combine.model.data.rows<-function(effects){
  
  if(!is.data.frame(effects))
    stop("Parsed model data expected to be of type data.frame!")
  
  if(nrow(effects)<2)
    return(effects)
  
  vardesc=effects$vardesc
  level=effects$level
  level.lag=effects$level.lag
  level[is.na(level)]=""
  vardesc=paste0(vardesc,"|",level,"(",level.lag,")")
  vars=unique(vardesc)
  effects2=new.effects.data.frame()
  
  for(i in seq_along(vars)){
    vd=vars[i]
    inds=which(vardesc==vd)
    if(length(inds)==1){
      effects2[i,]=effects[inds,]
    }else{
      effc=effects[inds[1],]
      effc$depends.on[[1]]=unique(unlist(lapply(inds,function(j)effects$depends.on[[j]])))
      modifs=lapply(inds,function(j)effects$modifiers[[j]])
      modifs=combineModifiers(modifs,vd)
      effc$modifiers[[1]]=modifs
      
      latnames=sapply(inds, function(j)effects$latName[j])
      latnames=latnames[!is.na(latnames)]
      latnames=unique(latnames[latnames!=""])
      if(length(latnames)>1)
        gmlSEMerror("multiple names assigned for the same term at level :",effects$level[inds[1]])
      
      effc$latName=latnames[1]
      
      if(!effc$modified[1])
        effc$modifiers[[1]]=NA
      
      effects2[i,]=effc
    }
  }
  effects2
}
