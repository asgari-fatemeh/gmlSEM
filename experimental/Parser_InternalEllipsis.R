# Dose not properly support error handling as in:
# sam8 = '(x11,x12) for (y11,y12),(x21,x22) for (y11,y12),...,(x81,x82) for (y81,y82)' 
# sam9 = '(x11,x12 for y11,y12),(x21,x22 for y21,y22),...,(x81,x82 for y81,y82)' 
# expand.ellipsis(sam8)
# expand.ellipsis(sam9)

"%+%"<-function(a,b){paste0(a,b)}
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

################# expanding arithmetic sequences #################
##### generating simple numeric and char arithmetic sequences ####
# x a sequenc of length 2 or 3, denoting the first (one or two) and the last element of seq.
# type=c("char","numeric")
# ind = 2 or 3, denoting the place of ellipsis in sequence (equals to the length of x)
get_seq<-function(x,type,ind=length(x)){
  #ind is the index of ... in sequence x
  i0<-1
  i1<-ind-1
  i2<-ind
  
  if(type=="char"){
    xs<-sapply(x, function(xi){as.numeric(charToRaw(xi))})
    a<-head(xs,1)
    b<-tail(xs,1)
    
    if(ind==2)
      by<-ifelse(a<b,1,-1)
    if(ind==3)
      by<-ifelse(a<b,xs[2]-xs[1],xs[1]-xs[2])
    
    letters<-rawToChar(as.raw(seq(a,b,by=by)),multiple = TRUE)
    for(i in i0:i1)
      if(letters[i]!=x[i]){
        return(stopp())
      }
        
    
    for(i in i2:length(x))
      if(letters[length(letters)+i-length(x)]!=x[i]){
        return(stopp())
      }
        
    
    
    return(list(list(len=length(letters),seq=letters,amb=FALSE,name.pat=".",alternatives=list(letters))))
  }else if(type=="numeric"){
    res=list()
    ri=1
    
    xs=sapply(x, function(xi){as.numeric(xi)})
    xs2=xs#[-ind]
    if(all(xs2==xs2[1]))
      stopp("The sequence is constant!")
    
    
    ln=nchar(as.character(abs(xs2[1])))
    combs=rev(lapply(0:(ln-1), function(l)combn(ln,l)))
    xs2s=paste0(xs2,collapse=",")
    
    for(cl in combs){
      
      for(j in 1:ncol(cl)){
        fixed.inds=cl[,j]
        name.pat2=strsplit(as.character(abs(xs2[1])),"")[[1]]
        name.pat=rep(".",length(name.pat2))
        name.pat[fixed.inds]=name.pat2[fixed.inds]
        name.pat=paste0(name.pat,collapse = "")
        pat=genPattern(ln,fixed.inds,length(xs2))
        
        if(!grepl(pat,xs2s,perl = TRUE))
          next
        
        m=gregexpr(pat,xs2s,perl = TRUE)[[1]]
        cap.start=attr(m,"capture.start")
        cap.len=attr(m,"capture.length")
        
        cil=(length(cap.start)/length(xs2))
        cij=1
        seqq=list()  
        skip=FALSE
        for(ci in 1:cil){
          if(skip)
            break
          
          
          x1=substr(xs2s,cap.start[ci],cap.start[ci]+cap.len[ci]-1)
          x2=substr(xs2s,cap.start[ci+cil],cap.start[ci+cil]+cap.len[ci+cil]-1)
          x3=substr(xs2s,cap.start[ci+2*cil],cap.start[ci+2*cil]+cap.len[ci+2*cil]-1)
          
          xs=as.numeric(c(x1,x2,x3))
          if(ind==2)
            xs=xs[!is.na(xs)]
          
          if(colnames(cap.start)[ci]==paste0("d",cij)){
            seqq[[ci]]=x1
            cij=cij+1
            next
          }
          
          a=head(xs,1)
          b=tail(xs,1)
          
          #xs must be descending or ascending
          xsd=sort(xs,decreasing = TRUE)
          xsa=sort(xs,decreasing = FALSE)
          
          if(any(xsd!=xs)&&any(xsa!=xs)){
            #Return NULL seq
            break
          }
          
          #xs must be an arithmetic sequence
          if(ind>2){
            a1=xs[1]
            b1=xs[2]-xs[1]
            n1=(tail(xs,1)-a1)/b1
            if(n1!=floor(n1) || n1<2){
              break
            }
          }
          
          
          
          if(ind==2){
            by=ifelse(a<b,1,-1)
          }
          
          if(ind>2){
            by=ifelse(a<b,xs[2]-xs[1],xs[1]-xs[2])
          }
          
          if((abs(a-b)/abs(by))>1000)
            break
          
          skip=FALSE
          letters=seq(a,b,by=by)
          for(i in i0:i1)
            if(letters[i]!=xs[i])
              skip=TRUE
          if(skip)
            break
          for(i in i2:length(xs))
            if(letters[length(letters)+i-length(xs)]!=xs[i])
              skip=TRUE
          if(skip)
            break
          
          seqq[[ci]]=letters
        }
        
        if(skip)
          next
        
        lens=sapply(seqq, function(s)length(s))
        lens=lens[lens!=1]
        if(any(lens!=lens[1]))
          next
        lens=lens[1]
        
        seq1=rep("",lens)
        for(i in 1:length(seqq)){
          seq1=paste0(seq1,seqq[[i]])
        }
        
        res[[ri]]=list(len=lens,seq=seq1,name.pat=name.pat)
        ri=ri+1
      }
    }
    
    
    
  }else{
    return(stopp())
  }
  
  if(length(res)==0){
    return(stopp())
  }
  
  res2=list(res[[1]])
  res2[[1]]$amb=FALSE
  res2[[1]]$alternatives=list(res[[1]]$seq)
  res2[[1]]$name.pat=res[[1]]$name.pat
  jj=1
  
  if(length(res)>1){
    for(i in 2:length(res)){
      skip=FALSE
      
      #Looking for ambiguity
      for(j in 1:length(res2)){
        
        if(res[[i]]$len==res2[[j]]$len){
          skip=TRUE
          dup=FALSE
          for(k in 1:length(res2[[j]]$alternatives)){
            if(all(res2[[j]]$alternatives[[k]]==res[[i]]$seq)){
              dup=TRUE
              break
            }
          }
          if(!dup){
            res2[[j]]$amb=TRUE
            res2[[j]]$alternatives[[length(res2[[j]]$alternatives)+1]]=
              res[[i]]$seq
          }
        }
        if(skip)
          break
      }
      
      
      if(!skip){
        jj=jj+1
        res2[[jj]]=res[[i]]
        res2[[jj]]$amb=FALSE
        res2[[jj]]$alternatives=list(res[[i]]$seq)
      }
    }
  }
  
  res2
  
}


####### generating perl compatible regexp for decomposing number sequences
genPattern<-function(x1len,fixed.inds=c(),seq_leng=2,sep=","){
  pat1=""
  pat2=""
  j=0
  k=0
  for(i in 1:x1len){
    if(i %in% fixed.inds){
      if(k>0){
        pat1=paste0(pat1,"(\\d{",k,"})")
        pat2=paste0(pat2,"(\\d{",k,",})")  
      }
      k=0
      j=j+1
      pat1=paste0(pat1,"(?'d",j,"'\\d{1})")
      pat2=paste0(pat2,"(\\k'd",j,"')")
    }else{
      k=k+1
    }
  }
  
  if(k>0){
    pat1=paste0(pat1,"(\\d{",k,"})")
    pat2=paste0(pat2,"(\\d{",k,",})")  
  }
  
  pat=paste0(pat1,sep,pat2)
  if(seq_leng>2){
    for(i in 3:seq_leng){
      pat=paste0(pat,sep,pat2)
    }
  }
  pat
}

setAttr<-function(txt,sep,prefix="",df=NULL,common.name=""){
  txts=strsplit(txt,sep,fixed=TRUE)[[1]]
  if(prefix!="")
    txt=paste0(prefix,txt)
  attr(txt,"len")=length(txts)
  attr(txt,"n1")=length(txts)
  attr(txt,"n2")=0
  attr(txt,"mat")=matrix(txts,ncol=1)
  if(common.name=="" && !is.null(df))
    common.name=paste0(sapply(1:nrow(df), function(i)ifelse(df$fixed[i],df[i,ncol(df)],".")),collapse = "")
  
  attr(txt,"common.name")=common.name
  txt
}


################### Expanding ellipsis ##########################
#################################################################
expand.ellipsis<-function(txt,multiple=FALSE,ignore.first=FALSE,details=FALSE,sep="",
                          allow.duplicate=FALSE #Allow duplicate elements in vector sequences
                          ){
  
  return.longest.possbile=TRUE  #Change it to FALSE if you want the shortest possible match
  
  latent.keywords="=>|->|for|as|=~|\\^"
  # keyword to introduce latents
  # x1 as (z1,y1),...,x9 as (z9,y9)      #zero-inflated distribution 
  # (x11,x12,x13) as y1,...,(x91,x92,x93) as y9
  # (x11,x12,x13) as (y11,y12),...,(x91,x92,x93) as (y91,y92)
  # (x11,x12,x13) -> y1,...,(x91,x92,x93) -> y9
  # (x11,x12,x13) -> (y11,y12),...,(x91,x92,x93) -> (y91,y92)
  # (x11,x12,x13) ~ y1,...,(x91,x92,x93) -> y9
  # (x11,x12,x13) =~ (y11,y12),...,(x91,x92,x93) -> (y91,y92)
  # (x11,x12,x13 as y11,y12),...,(x91,x92,x93) -> (y91,y92)
  
  #latent variables can be written in the following ways
  # (...) as  ...
  # (...) as (...)
  # (... as ...)
  
  txt.org<-txt
  txt=gsub("\\s+"," ",txt,perl=TRUE)
  
  getKeyWord=function(pat,txt){
    res=gregexpr(pat,txt,perl=TRUE)[[1]]
    c.s=attr(res,"capture.start")
    c.l=attr(res,"capture.length")
    em1=c.s!=0 & c.l!=0
    c.s=c.s[em1]
    c.l[em1]=c.l[em1]
    a=substr(txt,c.s,c.s+c.l-1)
    if(a=="^")
      a="\\^"
  }
  
  stopp<-function(msg=""){
    
    
    if(ignore.first)
      stop()
    
    tryCatch({
      if(!details){
        out=expand.ellipsis(txt.org,ignore.first = TRUE,multiple = multiple)
      }else{
        out=list(list(with.first=TRUE,failed=TRUE))
        out[[2]]=list(with.first=FALSE,failed=FALSE,
                      expand.ellipsis(txt.org,ignore.first = TRUE,multiple = multiple,details = details))
      }
      
      env=parent.env(environment())
      assign("out",out,env)
      assign("txt",out,env)
      # call <- quote(return(out)) 
      # #Return the result all along the caller function
       rlang::eval_bare(quote(return(out)), env = env)
    },error=function(e){
      stop("\ngmlSEM error: cannot parse the sequence\n",txt.org,"\n",msg,ifelse(is.na(msg)||msg=="","","\n")) 
    })
    return(out)
    
  }
  
  varPat="(?:(?:[\\w\\.][\\w\\._=]*\\s*,\\s*)*[\\w\\.][\\w\\._=]*)"     #Variable Pattern
  elem.place.holders=c(                                               #General Place holder patterns
    "\\(%1$s\\)",                                #First pattern represents only a vector, i.e. without latent naming conventions
    "(%4$s%1$s)\\s*%2$s\\s*(%4$s%3$s)",                #Next patterns match to latent naming conventions
    "\\((%4$s%1$s)\\)\\s*%2$s\\s*(%4$s%3$s)", 
    "(%4$s%1$s)\\s*%2$s\\s*\\((%4$s%3$s)\\)", 
    "\\((%4$s%1$s)\\)\\s*%2$s\\s*\\((%4$s%3$s)\\)",
    "\\((%4$s%1$s)\\s*%2$s\\s*(%4$s%3$s)\\)")
  
  
  elem.place.holders.capture=sprintf(elem.place.holders,"%1$s","%2$s","%3$s","")
  elem.place.holders.nocapture=sprintf(elem.place.holders,"%1$s","%2$s","%3$s","?:")
  
  
  latPat0=sprintf(elem.place.holders.nocapture,varPat,"%1$s",varPat)
  latPat=sprintf(elem.place.holders.capture,varPat,"%1$s",varPat)             #Latent definition patterns
  
  
  #Add parantheses around the pattern to capture the exact keyword in the text
  latent.keywords.capture="("%+%latent.keywords%+%")"
  latent.keywords.nocapture="(?:"%+%latent.keywords%+%")"
  
  #gmlSEM accepts 4 types of ellipsis
  #The first three types has the following patterns and exists in family: block statements
  # (x1.1,x1.2,x1.3),...,(x9.1,x9.2,x9.3)
  # (x11,x12,x13) as y1,...,(x91,x92,x93) as y9
  # (x11,x12,x13) as (y11,y12),...,(x91,x92,x93) as (y91,y92)
  
  # The function supress details=TRUE and ignores returning the details, if the template is complicated
  # pat00=sprintf(sprintf("(?:^|\\n)\\s*%1$s(?:\\s*,\\s*%1$s)*(?:\\n|$)",latpats),latent.keyword.capture)
  # pat01=sprintf("(?:^|\\n)\\s*(?:\\([,\\.\\w\\s]+\\)(?!\\s*%1$s)(?:,)\\s*)?\\([,\\.\\w\\s]+\\)(?!\\s*%1$s)\\s*(?:,)\\s*\\.\\.\\.\\s*(?:,)\\s*\\([,\\.\\w\\s]+\\)(?!\\s*%1$s)(?:\\n|$)",latent.keyword.capture)
  # pat02=sprintf("(?:^|\\n)\\s*(?:\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*[\\w\\.]+\\s*(?:,)\\s*)?\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*[\\w\\.]+\\s*(?:,)\\s*\\.\\.\\.\\s*(?:,)\\s*\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*[\\w\\.]+(?:\\n|$)",latent.keyword.capture)
  # pat03=sprintf("(?:^|\\n)\\s*(?:\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*\\([,\\.\\w\\s]+\\)\\s*(?:,)\\s*)?\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*\\([,\\.\\w\\s]+\\)\\s*(?:,)\\s*\\.\\.\\.\\s*(?:,)\\s*\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*\\([,\\.\\w\\s]+\\)(?:\\n|$)",latent.keyword.capture)
  
  elem.pat0=paste0("\\s*",sprintf(latPat0,latent.keywords.nocapture),"\\s*")
  elem.pat0.flatten=paste0("(?:",paste0("(?:",elem.pat0,")"),")")
  seq.pat0=sprintf("(?:%1$s\\s*(?:\\s*(,)\\s*%1$s)*\\s*(,)\\s*\\.\\.\\.\\s*(,)\\s*%1$s\\s*(?:\\s*(,)\\s*%1$s\\s*)*)|(?:%1$s\\s*(?:\\s*(,)\\s*%1$s)*\\s*)|\\s*(?:%1$s)\\s*",elem.pat0.flatten) 
  pats0=paste0("(?:^|\\n)(?:\\s*",seq.pat0,"\\s*)(?:\\n|$)")
  
  
  pm=sapply(pats0, function(p)grepl(p,txt,perl=TRUE))
  
  
  if(any(pm)){
    pmi=which(pm)[1]
    
    
    if(!grepl("...",txt,fixed=TRUE)){
      pats0=paste0(elem.pat0.flatten[pmi],"\\s*(,)")
    }else{
      pats0=pats0[pmi]
    }
    
    text=perlsplit(txt,pats0)
    
    
    if(pm[1]){
      n1=length(strsplit(text[1],",")[[1]])
      n2=0
    }else{
      latPatKeywordIdent=sprintf(latPat0[pmi],latent.keywords.capture)
      latent.keyword=getKeyWord(latPatKeywordIdent,text[1])
      latent.keyword.nocapture="(?:"%+%latent.keyword%+%")"
      latent.keyword.capture="("%+%latent.keyword%+%")"
      
      txt0=perlsplit(text[1],latPatKeywordIdent)
      n1=length(strsplit(txt0[1],",")[[1]])
      n2=length(strsplit(txt0[2],",")[[1]])  
    }
    
    mat=matrix(".",nrow = length(text),ncol=n1+n2 )
    
    
    if(!grepl("...",txt,fixed=TRUE)){
      for(i in 1:length(text)){
        txt0=perlsplit(text[i],latPatKeywordIdent)
        txt0=trim(strsplit(gsub("(","",gsub(")","",paste0(txt0,collapse = ","),fixed = TRUE),fixed=TRUE),",")[[1]])
        mat[i,]=txt0
        for(j in 1:n1){
          if(grepl("^\\.+$",mat[i,j],perl=TRUE)){
            stopp("Missing variables are only accepted for underlying latent variables in the family: block.")
          }
        }
        for(j in (n1+1):(n1+n2)){
          if(grepl("^\\.*$",mat[i,j],perl=TRUE)){
            mat[i,j]=gsub("^\\.*$",".",mat[i,j],perl = TRUE)
          }
        }
      }
      
      attr(txt,"len")=length(text)
      attr(txt,"n1")=n1
      attr(txt,"n2")=n2
      attr(txt,"mat")=mat
      attr(txt,"common.name")=""
      return(txt)
    }
    
    txts=matrix(character(),nrow=length(text),ncol=n1+n2)
    for(i in 1:length(text)){
      txt1=trim(text[i])
      if(txt1=="..."){
        txts[i,]="..."
        next
      }
      
      if(pmi>1){
        txt1=sapply(txt1, function(tx){xx=perlsplit(tx,latPatKeywordIdent);paste0(xx,collapse = ",")})        
      }

      txt1=trim(gsub("\\s","",gsub(")","",gsub("(","",txt1,fixed = TRUE),fixed = TRUE)))
      txts[i,]=strsplit(txt1,",")[[1]]
      
      if(any(txt1==".")){
        if(i<=n1){
          stopp("Missing variables are only accepted for underlying latent variables in the family: block.")
        }
      }
    }
    res=list()
    for(i in 1:ncol(txts)){
      txt0=paste0(txts[,i],collapse=",")
      if(txts[1,i]=="."){
        res[[i]]=list(list(failed=FALSE,lens=NA,missing=TRUE),
                      list(failed=FALSE,lens=NA,missing=TRUE))
      }else{
        res[[i]]=expand.ellipsis(txt0,ind,details=TRUE)  
      }
      
    }
    
    lens1=sort(unique(unlist(lapply(res, function(r)if(r[[1]]$failed){NA}else{r[[1]]$lens}))))
    lens2=sort(unique(unlist(lapply(res, function(r)if(r[[2]]$failed){NA}else{r[[2]]$lens}))))
    
    if(is.null(lens2)){
      lens2=integer(0)
    }
    
    for(i in seq_along(res)){
      if(exists('missing',res[[i]][[1]]))
        next
      
      lens1=intersect(lens1,res[[i]][[1]]$lens)
      if(length(lens2)>0)
        lens2=intersect(lens2,res[[i]][[2]]$lens)
    }
    
    if(allow.duplicate){
      if(return.longest.possbile){
        #Return the longest possible sequence
        lens1=sort(lens1,decreasing = TRUE)[1]
        lens2=ifelse(length(lens2)>0,sort(lens2,decreasing = TRUE)[1],NA)  
      }else{
        #Return the shortes possible sequence
        lens1=sort(lens1,decreasing = FALSE)[1]
        lens2=ifelse(length(lens2)>0,sort(lens2,decreasing = FALSE)[1],NA)   
      }
    }else{
      if(return.longest.possbile){
        #Return the longest possible sequence
        lens11=sort(lens1,decreasing = TRUE)
        lens22=if(length(lens2)>0){sort(lens2,decreasing = TRUE)}else{NA}
      }else{
        #Return the shortest possible sequence
        lens11=sort(lens1,decreasing = FALSE)
        lens22=if(length(lens2)>0){sort(lens2,decreasing = FALSE)}else{NA}
      }
        lens1=lens2=NA
        ki=1
        for(jm in rev(seq_along(lens11))){
          lens=lens11[jm]
          seqs=lapply(res, function(r)if(r[[ki]]$failed){NA}else{r[[ki]]$seqs[[which(sapply(r[[ki]]$seqs,function(s)s$len==lens))[1]]]})
          if(!any.intersect(seqs)){
            lens1=lens
            break
          }
        }
        ki=2
        if(length(lens2)>0)
        for(jm in rev(seq_along(lens22))){
          lens=lens22[jm]
          seqs=lapply(res, function(r)if(r[[ki]]$failed){NA}else{r[[ki]]$seqs[[which(sapply(r[[ki]]$seqs,function(s)s$len==lens))[1]]]})
          if(!any.intersect(seqs)){
            lens2=lens
            break
          }
        }
      
    }
    
    
    
    wf1=sapply(res, function(r)r[[1]]$failed)
    wf2=if(!is.na(lens2)){sapply(res, function(r)r[[2]]$failed)}else{NA}
    
    if(all(!is.na(wf1))&&!is.null(lens1)){
      ki=1
      lens=lens1
    }else if(all(!is.na(wf2))&&!is.null(lens2)){
      ki=2
      lens=lens2
    }else{
      return(stopp())
    }
    
    amb=FALSE
    xss=list()
    xss.amb=list()
    for(i in 1:(n1+n2)){
      if(exists('missing',res[[i]][[ki]])){
        xss[[i]]=rep(".",lens)
        xss.amb[[i]]=rep(".",lens)
      }else{
        ind=which(sapply(res[[i]][[ki]]$seqs,function(s)s$len==lens))
        amb=amb||res[[i]][[ki]]$seqs[[ind]]$amb
        xss[[i]]=res[[i]][[ki]]$seqs[[ind]]$seq
        xss.amb[[i]]=res[[i]][[ki]]$seqs[[ind]]$seqamb  
      }
    }
    
    tx=rep("",length(xss))
    mat=matrix(".",nrow = length(xss[[1]]),ncol=n1+n2)
    .tx=rep("",length(xss))
    .tx.amb=rep("",length(xss))
    
    for(j in 1:length(xss[[1]])){
      .tx[j]=.tx.amb[j]=tx[j]="("
      for(i in 1:n1){
        mat[j,i]=xss[[i]][j]
        sep=ifelse(i==n1,"",",")
        tx[j]=tx[j]%+%xss[[i]][j]%+%sep
        .tx[j]=.tx[j]%+%xss[[i]][j]%+%sep
        .tx.amb[j]=.tx.amb[j]%+%xss.amb[[i]][j]%+%sep
      }
      tx[j]=tx[j]%+%")"  
      .tx[j]=.tx[j]%+%")"  
      .tx.amb[j]=.tx.amb[j]%+%"" 
      
      if(n2==0)
        next
      
      tx[j]=tx[j]%+%"^("  
      .tx[j]=.tx[j]%+%" "%+%latent.keyword%+%" ("  
      .tx.amb[j]=.tx.amb[j]%+%" "%+%latent.keyword%+%" ("  
      for(i in (n1+1):(n1+n2)){
        mat[j,i]=xss[[i]][j]
        sep=ifelse(i==(n1+n2),"",",")
        tx[j]=tx[j]%+%xss[[i]][j]%+%sep
        .tx[j]=.tx[j]%+%xss[[i]][j]%+%sep
        .tx.amb[j]=.tx.amb[j]%+%xss.amb[[i]][j]%+%sep
      }
      tx[j]=tx[j]%+%")"  
      .tx[j]=.tx[j]%+%")"  
      .tx.amb[j]=.tx.amb[j]%+%")"  
    }
    
    if(amb){
      stopp(paste0("The sequence is ambiguist. It has at least two non-distinguishable alternatives:\n",
                   paste0(.tx,collapse = ","),"\n",
                   paste0(.tx.amb,collapse = ",")))
    }
    lentx=length(tx)
    tx=paste0(tx,collapse = ",")
    
    attr(tx,"len")=lentx
    attr(tx,"mat")=mat
    attr(tx,"n1")=n1
    attr(tx,"n2")=n2
    return(tx)
    
  }
  
  
  
  if(sep==""){
    grp<-gregexpr("(?J)(?>(?'sep'[,\\+\\*/])?\\s*\\.\\.\\.\\s*\\k'sep'|[^\\s,+]+(?'sep'[,\\+\\*/])(?: *\\k'sep' *[^\\s,+]+)*)",txt,perl=TRUE)
    capture.start<-attr(grp[[1]],"capture.start")
    ind=which(capture.start>0)[1]
    sep<-substr(txt,capture.start[ind],capture.start[ind])  
    if(is.na(sep)){
      sep="," #Does not harm the sequencing. The input simply is a single element
    }
      
  }
  
  
  
  #The sequence might be started with the sep char, in the case that it has splited from a longer sequence
  prefix=""
  if(substring(txt,1,1)==sep){
    prefix=sep
    txt=substring(txt,2)
  }
  
  # if(!grepl("...",txt,fixed=TRUE)){
  #   txt=setAttr(txt,sep,prefix)
  #   return(txt)
  # }
  
  txt<-gsub("\\s","",txt,perl=TRUE)
  text<-strsplit(txt,sep,fixed = TRUE)[[1]]
  elems<-list()
  
  #At most two elements before ... will be used
  #and we only proceed for one elem after ...
  foundEllip<-FALSE
  AddedNextElem<-FALSE
  ki<-which(text=="...")[1]
  if(!is.na(ki) && (ki==1|ki==length(text)))
    return(stopp())
  ki=max(ki-2,1)
  ki2=ifelse(is.na(ki),1,ki)
  prefix<-paste0(prefix,
                 ifelse(is.na(ki) || ki==1,"",paste0(paste0(text[1:(ki-1)],collapse=sep),sep)))
  
  
  
  suffix<-""
  i<-0
  for(ik in ki2:length(text)){
    i<-i+1
    txt<-text[ik]
    if(txt=="..."){
      foundEllip<-TRUE
      next
    }
    if(foundEllip&AddedNextElem){
      suffix<-paste0(sep,paste0(text[ik:length(text)],collapse=sep))
      break
    }
    if(foundEllip)
      AddedNextElem<-TRUE
    
    
    elems[[i]]<-list()
    elems[[i]]$txt<-txt
    grp<-gregexpr("[a-zA-Z]|_|\\.|\\d+|\\*|\\(|\\)|\\[|\\]|\\\"|\\'|-|=|\\+|@|\\^|%|\\!|\\$|\\&|\\\\|/|:|\\;",txt,perl=TRUE)
    match.length<-attr(grp[[1]],"match.length")
    match.start<-grp[[1]]
    
    elems[[i]]$parts<-sapply(seq_along(match.start),
                             function(j){substr(txt,match.start[j],match.start[j]+match.length[j]-1)})
  }
  
  
  count<-sapply(seq_along(elems), 
                function(i)ifelse(is.null(elems[[i]]),NA,length(elems[[i]]$parts)))
  count<-count[!is.na(count)]
  
  if(ignore.first)
    count<-count[-1]
  
  all.equals<-function(x){
    all(x==rep(x[1],length(x)))
  }
  
  if(!all.equals(count))
    return(stopp())
  
  ellipsis.idx<-which(sapply(seq_along(elems),function(i)is.null(elems[[i]])))
  
  #If failed, retrying with dropping first element
  
  getTypes<-function(xs){
    sapply(xs, function(x){
      if(!is.na(suppressWarnings(as.numeric(x))))
        return("numeric")
      if(grepl("[a-zA-Z]",x,perl = TRUE))
        return("char")
      return(x)
    })
  }
  
  i0<-ifelse(ignore.first,2,1)
  df<-data.frame(id=seq_along(elems[[i0]]$parts),
                 val=elems[[i0]]$parts,
                 type=getTypes(elems[[i0]]$parts),
                 fixed=TRUE)
  
  k<-0
  kk<-ncol(df)
  kkk<-ifelse(ignore.first,1,0)
  seq_elems= seq_along(elems)
  if(!is.na(ki))
    seq_elems=seq_elems[-ellipsis.idx]
  for(i in seq_elems){
    k<-k+1
    el<-elems[[i]]
    df[,k+kk]<-NA
    colnames(df)[k+kk]<-i
    if(i<i0)
      next
    for(j in 1:nrow(df))
      df[j,k+kk]<-ifelse(df$type[j]=="numeric",as.numeric(el$parts[j]),el$parts[j])
  }
  
  for(i in 1:nrow(df)){
    if(!all.equals(unname(unlist(df[i,-(1:(kk+kkk))])))){
      df$fixed[i]<-FALSE
    }
  }
  
  get_seq_df<-function(df,i){
    #tryCatch({
      x<-unname(unlist(df[i,-(1:kk)]))
      if(!ignore.first){
        seqq<-get_seq(x,df$type[i],ellipsis.idx)    
      }else{
        seqq<-get_seq(x[-1],df$type[i],ellipsis.idx-1)    
      }
      
      return(seqq)
    #},error=function(e)return(list(list(len=0,seq=numeric(),amb=TRUE,alternatives=list(numeric())))))
  }
  
  # if(is.na(len)){
  #   id1=which(!df$fixed)[1]  
  #   len=length(get_seq_df(df,id1))
  # }
  
  # seqq=rep("",len)
  
  if(is.na(ki)){  #No ellipsis to expand
    txt=setAttr(txt,sep,prefix,df=df)
    return(txt)
  }
  
  seqqs=list()
  lens=list()
  common.name=""
  for(i in 1:nrow(df)){
    if(df$fixed[i]){
      #seqq<-paste0(seqq,df[i,ncol(df)])
      seqqs[[i]]=list(fixed=TRUE,lens=NA,val=df[i,ncol(df)],res=list())
      common.name=paste0(common.name,df[i,ncol(df)])
      lens[[i]]=NA
    }else{
      tryCatch({
      seq2<-get_seq_df(df,i)
      },error=function(e){
        return(stopp())
        })
      # if(length(seq2)!=len)
      #   return(stopp())
      #seqq<-paste0(seqq,seq2)
      common.name=paste0(common.name,seq2$common.name)
      
      lens[[i]]=sapply(seq2, function(s)s$len)
      seqqs[[i]]=list(fixed=FALSE,lens=lens[[i]],val=NA,res=seq2)
    }
  }
  
  expand.all.seqs<-function(xs){
    d=unlist(lapply(xs, function(x)x$lens))
    d=d[!is.na(d)]
    if(length(d)==0)
      return(stopp())
    
    for(i in 1:length(xs)){
      if(xs[[i]]$fixed)
        next
      
      d=intersect(d,xs[[i]]$lens)    
    }
    
    d=sort(d)
    
    inds=which(sapply(xs, function(x) !x$fixed))
    xs2=lapply(xs, function(x){
      if(x$fixed){
        ls=list()
        for(di in seq_along(d)){
          ls[[di]]=list(amb=FALSE,name.pat=x$val,seq=rep(x$val,d[di]),seqamb=rep(x$val,d[di]))
        }
        return(ls)
      }
      
      
      ls=list()
      for(di in seq_along(d)){
        ind=which(sapply(x$res, function(xr)xr$len==d[di]))
        ind2=1
        if(length(x$res[[ind]]$alternatives)>1)
          ind2=2
        ls[[di]]=list(amb=x$res[[ind]]$amb,name.pat=x$res[[ind]]$name.pat,seq=x$res[[ind]]$seq,seqamb=x$res[[ind]]$alternatives[[ind2]])
      }
      ls
    })
    
    res=list()
    
    for(di in seq_along(d)){
      amb=FALSE
      seq=""
      seqamb=""
      name.pat=""
      for(i in 1:length(xs2)){
        amb=amb||xs2[[i]][[di]]$amb
        seq=paste0(seq,xs2[[i]][[di]]$seq)
        seqamb=paste0(seqamb,xs2[[i]][[di]]$seqamb)
        name.pat=paste0(name.pat,xs2[[i]][[di]]$name.pat)
      }
      if(ignore.first){
        seq=c(elems[[1]]$txt,seq)
        seqamb=c(elems[[1]]$txt,seqamb)
      }
      
      res[[di]]=list(amb=amb,len=length(seq),seq=seq,seqamb=seqamb,common.name=name.pat)
    }
    res
  }
  
  res=expand.all.seqs(seqqs)
  
  if(!details){
    if(return.longest.possbile){
      #Return the longest possible sequence
      ind=which.max(sapply(res, function(r)r$len))   
    }else{
      #Return the shortes possible sequence
      ind=which.min(sapply(res, function(r)r$len))   
    }
    
    
    seqq=res[[ind]]$seq
    common.name=res[[ind]]$common.name
  }else{
    if(ignore.first)
      return(res)
    lens1=sapply(res,function(r)r$len)
    if(length(elems)>3){
      res2=expand.ellipsis(txt.org,ignore.first = TRUE,details = TRUE)
      lens2=sapply(res2,function(r)r$len)  
      res0=list(list(with.first=TRUE,failed=FALSE,seqs=res,lens=lens1),
                list(with.first=FALSE,failed=FALSE,seqs=res2,lens=lens2))
    }else{
      res0=list(list(with.first=TRUE,failed=FALSE,seqs=res,lens=lens1),
                list(with.first=FALSE,failed=TRUE,seqs=c(),lens=list()))
    }
    
    
    return(res0)
  }
  
  if(multiple)
    return(seqq)
  
  txt<-paste0(prefix,
              paste0(seqq,collapse = sep),
              suffix)
  
  #txt<-expand.ellipsis(txt,multiple = FALSE)
  
  #if(multiple)
  #  return(strsplit(txt,sep))
  #attr(txt,"len")=length(seqq)
  
  if(details==TRUE &&grepl("...",txt,fixed = TRUE)){
    stopp("The parser cannot expand ellipsis. Try to use a more simple syntax or avoid ambiguity!")
  }
  k=0
  while(grepl("...",txt,fixed = TRUE)){
    k=k+1
    txt=expand.ellipsis(txt)
    if(k>10)
      stopp("The parser cannot expand ellipsis. Try to use a more simple syntax or avoid ambiguity!")
  }
  
  
  if(!details)
    txt=setAttr(txt,sep,df=df,common.name=common.name)
  return(txt)
  
}

