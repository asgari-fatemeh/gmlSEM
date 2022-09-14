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
      if(letters[i]!=x[i])
        stopp()
    
    for(i in i2:length(x))
      if(letters[length(letters)+i-length(x)]!=x[i])
        stopp()
    
    
    return(list(list(len=length(letters),seq=letters,amb=FALSE,alternatives=list(letters))))
  }else if(type=="numeric"){
    res=list()
    ri=1
    
    xs=sapply(x, function(xi){as.numeric(xi)})
    xs2=xs#[-ind]
    if(all(xs2==xs2[1]))
      stopp("The sequence is constant!")
    
    
    ln=nchar(as.character(abs(xs2[1])))
    combs=lapply(0:(ln-1), function(l)combn(ln,l))
    xs2s=paste0(xs2,collapse=",")
    
    for(cl in combs){
      
      for(j in 1:ncol(cl)){
        fixed.inds=cl[,j]
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
        
        res[[ri]]=list(len=lens,seq=seq1)
        ri=ri+1
      }
    }
    
    
    
  }else{
    stopp()
  }
  
  res2=list(res[[1]])
  res2[[1]]$amb=FALSE
  res2[[1]]$alternatives=list(res[[1]]$seq)
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





################### Expanding ellipsis ##########################
#################################################################
expand.ellipsis<-function(txt,multiple=FALSE,ignore.first=FALSE,details=FALSE){
  
  return.longest.possbile=TRUE  #Change it to FALSE if you want the shortest possible match
  
  latent.keyword="=>|->|as|=~|=|~"
  # keyword to introduce latents
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
  
  varPat="(?:(?:[\\w\\.][\\w\\._]*\\s*,\\s*)*[\\w\\.][\\w\\._]*)"
  latpats=sprintf(c("(%1$s)\\s*%2$s\\s*(%1$s)",
                    "\\((%1$s)\\)\\s*%2$s\\s*(%1$s)",
                    "\\((%1$s)\\)\\s*%2$s\\s*\\((%1$s)\\)",
                    "\\((%1$s)\\)\\s*%2$s\\s*\\((%1$s)\\)",
                    "\\((%1$s)\\s*%2$s\\s*(%1$s)\\)"),varPat,"%1$s")
  
  
  #Add parantheses around the pattern to capture the exact keyword in the text
  latent.keyword.capture="("%+%latent.keyword%+%")"
  
  
  txt.org<-txt
  if(!grepl("...",txt,fixed=TRUE)){
    return(txt)
  }
  txt<-gsub("\\s","",txt,perl=TRUE)
  sep=""
  grp<-gregexpr("(?'sep'[,\\+\\*/])?\\s*\\.\\.\\.\\s*\\k'sep'",txt,perl=TRUE)
  capture.start<-attr(grp[[1]],"capture.start")
  sep<-substr(txt,capture.start[1],capture.start[1])
  
  #The sequence might be started with the sep char, in the case that it has splited from a longer sequence
  prefix=""
  if(substring(txt,1,1)==sep){
    prefix=sep
    txt=substring(txt,2)
  }
  
  len<-NA
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
      assign("out",out,parent.frame())
      call <- rlang::expr(return(out)) 
      #Return the result all along the caller function
      rlang::eval_bare(call, env = parent.frame())
    },error=function(e){
      stop("\ngmlSEM error: cannot parse the sequence\n",txt.org,"\n",msg,ifelse(is.na(msg)||msg=="","","\n")) 
    })
    
  }
  if(!grepl("(?'sep'[,\\+\\*/])\\s*\\.\\.\\.\\s*\\k'sep'",txt,perl=TRUE))
    stopp()
  
  
  txt=gsub("\\s+"," ",txt,perl=TRUE)
  
  #gmlSEM accepts 4 types of ellipsis
  #The first three types has the following patterns and exists in family: block statements
  # (x1.1,x1.2,x1.3),...,(x9.1,x9.2,x9.3)
  # (x11,x12,x13) as y1,...,(x91,x92,x93) as y9
  # (x11,x12,x13) as (y11,y12),...,(x91,x92,x93) as (y91,y92)
  
  # The function supress details=TRUE and ignores returning the details, if the template is complicated
  pat00=sprintf(sprintf("(?:^|\\n)\\s*%1$s(?:\\s*,\\s*%1$s)*(?:\\n|$)",latpats),latent.keyword.capture)
  pat01=sprintf("(?:^|\\n)\\s*(?:\\([,\\.\\w\\s]+\\)(?!\\s*%1$s)(?:,)\\s*)?\\([,\\.\\w\\s]+\\)(?!\\s*%1$s)\\s*(?:,)\\s*\\.\\.\\.\\s*(?:,)\\s*\\([,\\.\\w\\s]+\\)(?!\\s*%1$s)(?:\\n|$)",latent.keyword.capture)
  pat02=sprintf("(?:^|\\n)\\s*(?:\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*[\\w\\.]+\\s*(?:,)\\s*)?\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*[\\w\\.]+\\s*(?:,)\\s*\\.\\.\\.\\s*(?:,)\\s*\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*[\\w\\.]+(?:\\n|$)",latent.keyword.capture)
  pat03=sprintf("(?:^|\\n)\\s*(?:\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*\\([,\\.\\w\\s]+\\)\\s*(?:,)\\s*)?\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*\\([,\\.\\w\\s]+\\)\\s*(?:,)\\s*\\.\\.\\.\\s*(?:,)\\s*\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*\\([,\\.\\w\\s]+\\)(?:\\n|$)",latent.keyword.capture)
  
  
  perlsplit<-function(x,pat,drop.captured.sgroups.only=TRUE){
    #expect just one match
    if(!grepl(pat,x,perl=TRUE))
      return(x)
    
    if(!drop.captured.sgroups.only)
      return(strsplit(x,pat,perl= TRUE)[[1]])
    
    text=gregexpr(pat,x,perl= TRUE)[[1]]
    
    cap.start=attr(text,"capture.start")
    cap.len=attr(text,"capture.length")
    if(nrow(cap.start)>1)
      stopp()
    
    txtt=character()
    st=1
    j=1
    for(i in 1:ncol(cap.start)){
      if(cap.len[i]==0)
        next
      txtt[j]=substr(x,st,cap.start[i]-1)
      st=cap.start[i]+1
      j=j+1
    }
    if(st<=nchar(x))
      txtt[length(txtt)+1]=substr(x,st,nchar(x))
    
    txtt
  }
  
  p1m=grepl(pat01,txt,perl = TRUE)
  p2m=grepl(pat02,txt,perl = TRUE)
  p3m=grepl(pat03,txt,perl = TRUE)
  
  
  if(p1m||p2m||p3m){
    
    getKeyWord=function(pat){
      res=gregexpr(pat,txt,perl=TRUE)[[1]]
      c.s=attr(res,"capture.start")[2]
      c.l=attr(res,"capture.length")[2]
      substr(txt,c.s,c.s+c.l-1)
    }
    
    text=if(p1m){
      perlsplit(txt,pat1)
    }else if(p2m){
      latent.keyword=getKeyWord(pat02)
      latent.keyword.nocapture="(?:"%+%latent.keyword%+%")"
      latpats=sprintf(latpats,latent.keyword.nocapture)
      pat2=sprintf("(?:^|\\n)\\s*(?:\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*[\\w\\.]+\\s*(,)\\s*)?\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*[\\w\\.]+\\s*(,)\\s*\\.\\.\\.\\s*(,)\\s*\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*[\\w\\.]+(?:\\n|$)",latent.keyword.nocapture)
      perlsplit(txt,pat2)
    }else{
      latent.keyword=getKeyWord(pat03)
      latent.keyword.nocapture="(?:"%+%latent.keyword%+%")"
      pat3=sprintf("(?:^|\\n)\\s*(?:\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*\\([,\\.\\w\\s]+\\)\\s*(,)\\s*)?\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*\\([,\\.\\w\\s]+\\)\\s*(,)\\s*\\.\\.\\.\\s*(,)\\s*\\([,\\.\\w\\s]+\\)\\s*%1$s\\s*\\([,\\.\\w\\s]+\\)(?:\\n|$)",latent.keyword.nocapture)
      perlsplit(txt,pat3)
    }
    
    
    if(p1m){
      n1=length(strsplit(text,",")[[1]])
      n2=0
    }else{
      txt0=strsplit(text[1],"(?:\\s+|\\))as(?:\\s+|\\()",perl=TRUE)[[1]]
      n1=length(strsplit(txt0[1],",")[[1]])
      n2=length(strsplit(txt0[2],",")[[1]])  
    }
    
    txts=matrix(character(),nrow=length(text),ncol=n1+n2)
    for(i in 1:length(text)){
      txt1=trim(text[i])
      if(txt1=="..."){
        txts[i,]="..."
        next
      }
      txt1=gsub("(?:\\s+|\\))as(?:\\s+|\\()",",",txt1,perl = TRUE)
      txt1=trim(gsub("\\s","",gsub(")","",gsub("(","",txt1,fixed = TRUE),fixed = TRUE)))
      txts[i,]=strsplit(txt1,",")[[1]]
    }
    res=list()
    for(i in 1:ncol(txts)){
      txt0=paste0(txts[,i],collapse=",")
      res[[i]]=expand.ellipsis(txt0,ind,details=TRUE)
    }
    
    lens1=sort(unique(unlist(lapply(res, function(r)if(r[[1]]$failed){NA}else{r[[1]]$lens}))))
    lens2=sort(unique(unlist(lapply(res, function(r)if(r[[1]]$failed){NA}else{r[[2]]$lens}))))
    
    if(is.null(lens2)){
      lens2=c()
    }
    
    for(i in seq_along(res)){
      lens1=intersect(lens1,res[[i]][[1]]$lens)
      if(!is.null(lens2))
        lens2=intersect(lens2,res[[i]][[2]]$lens)
    }
    
    if(return.longest.possbile){
      #Return the longest possible sequence
      lens1=sort(lens1,decreasing = TRUE)[1]
      lens2=ifelse(length(lens2)>0,sort(lens2,decreasing = TRUE)[1],NA)
    }else{
      #Return the shortes possible sequence
      lens1=sort(lens1,decreasing = FALSE)[1]
      lens2=ifelse(length(lens2)>0,sort(lens2,decreasing = FALSE)[1],NA)   
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
      stopp()
    }
    
    amb=FALSE
    xss=list()
    xss.amb=list()
    for(i in 1:(n1+n2)){
      ind=which(sapply(res[[i]][[ki]]$seqs,function(s)s$len==lens))
      amb=amb||res[[i]][[ki]]$seqs[[ind]]$amb
      xss[[i]]=res[[i]][[ki]]$seqs[[ind]]$seq
      xss.amb[[i]]=res[[i]][[ki]]$seqs[[ind]]$seqamb
    }
    
    tx=rep("",length(xss))
    .tx=rep("",length(xss))
    .tx.amb=rep("",length(xss))
    
    for(j in 1:length(xss[[1]])){
      .tx[j]=.tx.amb[j]=tx[j]="("
      for(i in 1:n1){
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
    tx=paste0(tx,collapse = ",")
    
    return(tx)
    
  }
  
  text<-strsplit(txt,sep,fixed = TRUE)[[1]]
  elems<-list()
  
  #At most two elements before ... will be used
  #and we only proceed for one elem after ...
  foundEllip<-FALSE
  AddedNextElem<-FALSE
  ki<-which(text=="...")[1]
  if(ki==1|ki==length(text))
    stopp()
  ki<-max(ki-2,1)
  prefix<-paste0(prefix,
                 ifelse(ki==1,"",paste0(paste0(text[1:(ki-1)],collapse=sep),sep)))
  
  
  
  suffix<-""
  i<-0
  for(ik in ki:length(text)){
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
    stopp()
  
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
  for(i in seq_along(elems)[-ellipsis.idx]){
    k<-k+1
    el<-elems[[i]]
    df[,k+kk]<-NA
    colnames(df)[k+kk]<-i
    for(j in 1:nrow(df))
      df[j,k+kk]<-ifelse(df$type[j]=="numeric",as.numeric(el$parts[j]),el$parts[j])
  }
  
  for(i in 1:nrow(df)){
    if(!all.equals(unname(unlist(df[i,-(1:(kk+kkk))])))){
      df$fixed[i]<-FALSE
    }
  }
  
  get_seq_df<-function(df,i){
    tryCatch({
      x<-unname(unlist(df[i,-(1:kk)]))
      if(!ignore.first){
        seqq<-get_seq(x,df$type[i],ellipsis.idx)    
      }else{
        seqq<-get_seq(x[-1],df$type[i],ellipsis.idx-1)    
      }
      
      return(seqq)
    },error=function(e)stopp())
  }
  
  # if(is.na(len)){
  #   id1=which(!df$fixed)[1]  
  #   len=length(get_seq_df(df,id1))
  # }
  
  # seqq=rep("",len)
  seqqs=list()
  lens=list()
  for(i in 1:nrow(df)){
    if(df$fixed[i]){
      #seqq<-paste0(seqq,df[i,ncol(df)])
      seqqs[[i]]=list(fixed=TRUE,lens=NA,val=df[i,ncol(df)],res=list())
      lens[[i]]=NA
    }else{
      seq2<-get_seq_df(df,i)
      # if(length(seq2)!=len)
      #   stopp()
      #seqq<-paste0(seqq,seq2)
      lens[[i]]=sapply(seq2, function(s)s$len)
      seqqs[[i]]=list(fixed=FALSE,lens=lens[[i]],val=NA,res=seq2)
    }
  }
  
  expand.all.seqs<-function(xs){
    d=unlist(lapply(xs, function(x)x$lens))
    d=d[!is.na(d)]
    if(length(d)==0)
      stopp()
    
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
          ls[[di]]=list(amb=FALSE,seq=rep(x$val,d[di]),seqamb=rep(x$val,d[di]))
        }
        return(ls)
      }
      
      
      ls=list()
      for(di in seq_along(d)){
        ind=which(sapply(x$res, function(xr)xr$len==d[di]))
        ind2=1
        if(length(x$res[[ind]]$alternatives)>1)
          ind2=2
        ls[[di]]=list(amb=x$res[[ind]]$amb,seq=x$res[[ind]]$seq,seqamb=x$res[[ind]]$alternatives[[ind2]])
      }
      ls
    })
    
    res=list()
    for(di in seq_along(d)){
      amb=FALSE
      seq=""
      seqamb=""
      for(i in 1:length(xs2)){
        amb=amb||xs2[[i]][[di]]$amb
        seq=paste0(seq,xs2[[i]][[di]]$seq)
        seqamb=paste0(seqamb,xs2[[i]][[di]]$seqamb)
      }
      if(ignore.first){
        seq=c(elems[[1]]$txt,seq)
        seqamb=c(elems[[1]]$txt,seqamb)
      }
      
      res[[di]]=list(amb=amb,len=length(seq),seq=seq,seqamb=seqamb)
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
  
  return(txt)
  
}

