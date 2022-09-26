# TODO: Family does not fully support vector-valued distributions

source("experimental/Parser_InternalMisc.R")
source("experimental/Parser_InternalEllipsis.R")
source("experimental/Parser_InternalFamilyClass.R")
source("experimental/Parser_InternalFamilyCollection.R")
source("experimental/Parser_InternalObjects.R")

sample_syntax_file       = 'experimental/sample_gmlSEMScript.txt'
fileName = 'experimental/sample_gmlSEMScript.txt'
model.syntax      = readChar(sample_syntax_file, file.info(sample_syntax_file)$size)

#sample_syntax_extra_file = 'experimental/sample_gmlSEMScript_extra.txt'
#simulation.syntax = readChar(sample_syntax_extra_file, file.info(sample_syntax_extra_file)$size)
#model.syntax<-paste0(model.syntax,"\n",simulation.syntax)


#############################################################
##################### Parser Code ###########################


# check for empty syntax
if(length(model.syntax) == 0) {
  stop("gmlSEM ERROR: empty model syntax")
}

# remove comments prior to split:
# match from comment character to newline, but don't eliminate newline
model.syntax <- gsub("[#!].*(?=\n|$)", "", model.syntax, perl = TRUE)

model.syntax <- gsub(";", "\n", model.syntax, fixed = TRUE)
model.syntax <- gsub("\r\n", "\n", model.syntax, fixed = TRUE)
model.syntax <- gsub("\n(\\s*\n)?", "\n", model.syntax, perl = TRUE)
model.syntax <- gsub(pattern = "\u02dc", replacement = "~", model.syntax)

model.org=model.syntax

# Change block names to simpler  identifiable ones
getBlockPat<-function(block){#Create a pattern matches to block name and its abbreviations
 b=strsplit(block,"")[[1]] 
 pat="(?<=^|\n) *(?=.{1,})"
 pat=paste0(pat,b[1],"?")
 if((length(b)==2)){
   pat=paste0(pat,sprintf("(?:(?<=%s)%s|\\.)?",b[1],b[2]))
 }else if(length(b)>2){
   for(i in 2:(length(b)-1))
     pat=paste0(pat,sprintf("(?:(?<=%s)%s)?",substr(block,1,i-1),b[i]))
   i=length(b)
   pat=paste0(pat,sprintf("(?:(?<=%s)%s|\\.)?",substr(block,1,i-1),b[i]))
 }
pat=paste0(pat," *:")   
 
}

#model.syntax <- gsub("(?<=^|\n) *group:", ">>G>>", model.syntax, perl = TRUE,ignore.case =TRUE)
model.syntax <- gsub(getBlockPat("group"), ">>G>>:", model.syntax, perl = TRUE,ignore.case =TRUE)
model.syntax <- gsub(getBlockPat("level"), ">>L>>:", model.syntax, perl = TRUE,ignore.case =TRUE)
model.syntax <- gsub(getBlockPat("repeated"), ">>R>>:", model.syntax, perl = TRUE,ignore.case =TRUE)
model.syntax <- gsub(getBlockPat("heteroskedasticity"), ">>H>>:", model.syntax, perl = TRUE,ignore.case =TRUE)
model.syntax <- gsub(getBlockPat("family"), ">>F>>:", model.syntax, perl = TRUE,ignore.case =TRUE)
model.syntax <- gsub(getBlockPat("size"), ">>S>>:", model.syntax, perl = TRUE,ignore.case =TRUE)


# Sugar syntax for simulation
# Sample size can be set within size: blocks
pattern1="(?J)(?:(?<=^|\\n) *>>S>>:)(?'phrase'\\s*(?:(?'N'[^\\s:>]+)(?=\\b) *(?'level'[\\w\\.]+)(?=\\b)) *(?:(?:(?'keyword'per|in each|each) *(?'plevel'[\\w\\.]+)|(?'keyword'|in) *(?'plevel'[\\w\\.]+) *(?'plevelid'\\d+)))?)\\g'phrase'*"
pattern2="(?J)(?'phrase'\\s*(?'N'[^\\s:>]*) *(?'level'[\\w\\.]+)(?=\\b) *(?:((?'keyword'per|in each|each) *(?'plevel'[\\w\\.]+)|(?'keyword'in) *(?'plevel'[\\w\\.]+) *(?'plevelid'\\d+)))?)[,;]?"
res<-gregexpr(pattern1,model.syntax,ignore.case = TRUE, perl = TRUE)[[1]]
match.length<-attr(res,"match.length")
for(i in rev(seq_along(res))){
  txt<-substr(model.syntax,res[i],res[i]+match.length[i]-1)
  txt<-gsub(">>S>>:","",txt,ignore.case = TRUE)
  txt<-gsub("\n","",txt,fixed = TRUE)
  if(grepl(",\\s*\\.\\.\\.\\s*,",txt,perl=TRUE)){
    txt<-gsub(",\\s*\\.\\.\\.\\s*,",",...,",txt,perl=TRUE)
    txt<-trim(txt)
    txt<-switch.space(txt)
    txt<-expand.ellipsis(txt)
  }
  
  res2=gregexpr(pattern2,txt,ignore.case = TRUE, perl = TRUE)[[1]]
  vals=captured.groups.dataframe(txt,res2)
  vals[,1]=trim(vals[,1])
  parsedData[["size"]]=rbind(vals,parsedData[["size"]])
  txt<-sapply(txt, function(x)paste0(">>S>>:",vals[,1]))
  txt<-switch.space(txt)
  txt<-paste(txt,collapse = "\n")
  model.syntax<-paste0(substring(model.syntax,1,res[i]-1),
                       txt,
                       substring(model.syntax,res[i]+match.length[i]))
}



## Replacing keywords and phrases with operators and commands in gmlSEM syntax
#[as .] suffix name assign operator
capture.within <- gregexpr("\\[\\s*?as\\s+(.*?)\\s*\\]", model.syntax,perl = TRUE,ignore.case =TRUE)
match.start=capture.within[[1]]
match.length=attr(capture.within[[1]],"match.length")
capture.length=attr(capture.within[[1]],"capture.length")
capture.start<-attr(capture.within[[1]],"capture.start")
if(length(capture.start)>0){
  for(i in seq(length(capture.start),1,by=-1))
    model.syntax<-paste0(substring(model.syntax,1,match.start[i]-1),"@",
                         substring(model.syntax,capture.start[i],capture.start[i]+capture.length[i]-1),
                         substring(model.syntax,match.start[i]+match.length[i]))
}

#alias: 'as' keyword
capture.within <- gregexpr("(?:^|\n) *?((?>\\b|\\.)[\\w\\.\\+\\, ]+?\\b) *( as ) *(?! *,)(.*)", model.syntax, perl = TRUE,ignore.case =TRUE)
match.start=capture.within[[1]]
match.length=attr(capture.within[[1]],"match.length")
capture.length=attr(capture.within[[1]],"capture.length")
capture.length=as.matrix(capture.length,ncol=2)
capture.start<-attr(capture.within[[1]],"capture.start")
capture.start=as.matrix(capture.start,ncol=2)
if(nrow(capture.start)>0){
  for(i in seq(nrow(capture.start),1,by=-1)){
    lhs=trim(substring(model.syntax,capture.start[i,1],capture.start[i,1]+capture.length[i,1]))
    rhs=trim(substring(model.syntax,capture.start[i,3],capture.start[i,3]+capture.length[i,3]-1))
    
    lhs=gsub("\\s{1,}"," ",lhs,perl=TRUE)
    rhs=gsub("'","",gsub('"',"",gsub("\\s{1,}"," ",rhs,perl=TRUE))) #Removing quotation marks
    
    lhs=expand.ellipsis(lhs,sep=",")   #In alias statement it is expected that the elements are separated by a comma
    rhs=expand.ellipsis(rhs,sep=",")
    
    lhsm=attr(lhs,"mat")
    rhsm=attr(rhs,"mat")
    
    if(length(lhsm)!=length(rhsm))
      stop("\ngmlSEM error: Error in parsing alias statement:\n",
           substring(model.syntax,capture.start[i,1],capture.start[i,3]+capture.length[i,3]),
           "\nThe length of left hand side and right hand side are not equal")
    
    parsedData[["alias"]]=rbind(data.frame(lhs=lhsm,rhs=rhsm),
                                parsedData[["alias"]])
    
    model.syntax<-paste0(substring(model.syntax,1,capture.start[i,1]-1),">>A>>:",
                         lhs,"<>",switch.space(rhs),
                         substring(model.syntax,match.start[i]+match.length[i]))
    
    add.alias(lhsm,rhsm)
    
  }
}


#replacing '(vary|varies|varying) at level' with <<~
capture.within <- gregexpr("(?<=^|\n) *([\\w\\.\\,\\+ ]+) *(?:(?:vary|varies|varying) at level) *(.*) *(?=$|\n)", model.syntax, perl=TRUE,ignore.case =TRUE)
match.start<-capture.within[[1]]
if(capture.within[[1]][1]>0){
  capture.start=attr(capture.within[[1]],"capture.start")
  capture.length=attr(capture.within[[1]],"capture.length")
  
  if(length(capture.start)>0){
    for(i in seq(nrow(capture.start),1,by=-1)){
      if(any(capture.length[i,]==0))
        stop("\ngmlSEM error: invalid syntax in alias statement\n",
             substring(model.syntax,match.start[i],capture.start[i,2]+capture.length[i,2]-1))
      
     txt=trim(substring(model.syntax,capture.start[i,1],capture.start[i,1]+capture.length[i,1]-1))
     txt=expand.ellipsis(txt,sep=",") #It is expected that elements to be seperated by comma
     lev=trim(substring(model.syntax,capture.start[i,2],capture.start[i,2]+capture.length[i,2]-1))
       
      txtm=trim(c(attr(txt,"mat")))
      txtm.par=txtm
      lev.par=lev
      if(grepl("(",lev,fixed = TRUE)){
        levv=captured.groups(lev," *(?'var'.*) *\\( *(?'varn'.*) *\\)")
        lev=levv$var
        lev.par=levv$varn
      }
      
      for(j in 1:length(txtm)){
        if(grepl("(",txtm[j],fixed=TRUE)){
          levv=captured.groups(txtm[j]," *(?'var'.*) *\\( *(?'varn'.*) *\\)")
          txtm[j]=levv$var
          txtm.par[j]=levv$varn
        }
      }
    
      
      txtv1=!is.valid.varname(txtm)
      txtv2=!is.valid.varname(txtm.par)
      lev1=!is.valid.varname(lev)
      lev2=!is.valid.varname(lev.par)
      
      if(any(txtv1)||any(txtv2)||lev1||lev2){
        tt=c(txtm,txtm.par,lev,lev.par)
        ttt=c(txtv1,txtv2,lev1,lev2)
        tx=tt[which(ttt)[1]]
        stop("\ngmlSEM error: Invalid variable name '",tx,"' in alias statement:\n",
             substring(model.syntax,match.start[i],capture.start[i,2]+capture.length[i,2]-1))
      }
      
      model.syntax<-paste0(substring(model.syntax,1,match.start[i]-1),">>V>>:",
                           txt,"<<~",lev,
                           substring(model.syntax,capture.start[i,2]+capture.length[i,2]))
      
      lev=get.alias.lhs(lev)
      add.alias(lev,lev.par,"OV")
      add.levels(lev)
      
      parsedData[["vars.level"]]=rbind(data.frame(var=txtm,var.inpar=txtm.par,level=lev),
                                       parsedData[["vars.level"]])
    }
    
  }  
}


#'within' keyword in level: command
while(TRUE){
  capture.within <- gregexpr("(?<=>>L>>:)(\\s*[\\w\\.]+ *(?:\\( *[\\w\\.]+ *\\))?)+(?=$|\n)", model.syntax,perl = TRUE,ignore.case =TRUE)
  if(capture.within[[1]][1]==-1)
    break
  
  capture.start<-attr(capture.within[[1]],"capture.start")
  capture.length<-attr(capture.within[[1]],"capture.length")
  match.length<-attr(capture.within[[1]],"match.length")
  match.start<-capture.within[[1]]
  if(length(match.start)>0){
    for(i in seq(length(match.start),1,by=-1)){
      level_text<-substring(model.syntax,capture.start[i],capture.start[i]+capture.length[i])
      level_text<-gsub(" within ","<<",level_text,fixed=TRUE)
      model.syntax<-paste0(substring(model.syntax,1,capture.start[i]-1),"<<",
                           level_text,
                           substring(model.syntax,match.start[i]+match.length[i]))
    }
    
  }  
}

#Rectify lev names with respect to new aliases
if(!is.null(nrow(parsedData[["vars.level"]]))){
  parsedData[["vars.level"]]$level=get.alias.lhs(parsedData[["vars.level"]]$level)
}
if(nrow(levels.matrix)>0){
  colnames(levels.matrix)=rownames(levels.matrix)=get.alias.lhs(rownames(levels.matrix))
}


#marking distributions with tag >>F>>:
#template allow breaklines after the list of variables
pattern='(?<=\\n|^) *>>F>>:\\s*([\\w\\._\\,\\+\\s()\\*]*?)\\s*([\\w\\.][\\w\\._]*?(?:\\((?:[^()]|(?2))*\\)) *?)\\s*([\\w\\.][\\w\\._]*?(?:\\((?:[^()]|(?3))*\\)) *?)?(?=\\n|$)'
capture.within <- gregexpr(pattern, model.syntax, perl=TRUE,ignore.case =TRUE)
match.start<-capture.within[[1]]
match.length<-attr(capture.within[[1]],"match.length")
capture.start<-attr(capture.within[[1]],"capture.start")
capture.length<-attr(capture.within[[1]],"capture.length")
capture.start=as.matrix(capture.start,ncol=3)
capture.length=as.matrix(capture.length,ncol=3)
if(nrow(capture.start)>0){
  for(i in seq(nrow(capture.start),1,by=-1)){
    txt2<-""
    if(capture.start[i,3]>0)
      txt2<-paste0(">>F>>",substring(model.syntax,capture.start[i,3],capture.start[i,3]+capture.length[i,3]-1))
    
      txt1<-paste0(">>F>>",substring(model.syntax,capture.start[i,2],capture.start[i,2]+capture.length[i,2]-1))
      
      seqq=trim(substring(model.syntax,capture.start[i,1],capture.start[i,1]+capture.length[i,1]-1))
      seqq=strsplit(seqq,"\n",fixed = TRUE)[[1]]
      for(j in seq_along(seqq)){
        seqq[j]=gsub("(^,|,$)","",trim(seqq[j]),perl = TRUE)
        seqq[j]=expand.ellipsis(seqq[j])
      }
      seqq=paste0(seqq,collapse = ",")
      txt<-paste0(seqq,txt1,txt2)
      txt<-gsub("\n","",txt,perl  = TRUE)
      
    model.syntax<-paste0(substring(model.syntax,1,match.start[i]-1),">>F>>:",txt,"\n",
                         substring(model.syntax,match.start[i]+match.length[i]))    
  }

}

# remove all whitespace prior to split
model.syntax <- gsub("[ \t]+", "", model.syntax, perl = TRUE)
# remove any occurrence of >= 2 consecutive newlines to eliminate
# blank statements; this retains a blank newline at the beginning,
# if such exists, but parser will not choke because of start.idx
model.syntax <- gsub("\n{2,}", "\n", model.syntax, perl = TRUE)



#Expanding three-dots
#pattern = "[\\w\\.][^:,\\+\\n\\b\\b\\<\\>~\\-]*(?'sep'[,\\+])(?:[^:,\\+\\n\\b\\b\\<\\>~\\|\\-]+\\k'sep')?\\.\\.\\.\\k'sep'[^:,\\+\\n\\b\\b\\<\\>~\\|\\-]*[\\w\\.]"
pattern  = "(?J)(?>(?>(?>([^:,\\+\\n\\r\\<\\>~\\-\\(\\)|]|(?'cond'\\|))*(?>\\((?>(?!\\.\\.\\.)([^:,\\+\\n\\r\\<\\>~\\-\\(\\)|]|\\k'cond')|(?1))*\\))?){1,2})(?'sep'[,+]))?(?>(?>(?!\\.\\.\\.)([^:,\\+\\n\\r\\<\\>~\\-\\(\\)|]|\\k'cond')*(?:\\((?>[^:,\\+\\n\\r\\<\\>~\\-()|]|(?1))*\\))?){1,2})(?'sep'[,+])\\.\\.\\.\\k'sep'(?>(?>(?!\\.\\.\\.)[^:,\\+\\n\\r\\<\\>~\\-\\(\\)|]*(?>\\((?>(?!\\.\\.\\.)([^:,\\+\\n\\r\\<\\>~\\-\\(\\)|]|\\k'cond')|(?1))*\\))?){1,2})"
capture.ellipsis <- gregexpr(pattern, model.syntax, perl=TRUE,ignore.case =TRUE)
match.start<-capture.ellipsis[[1]]
match.length<-attr(capture.ellipsis[[1]],"match.length")
for(i in rev(seq_along(match.start))){
  txt.ellipsis<-substring(model.syntax,match.start[i],
                          match.start[i]+match.length[i]-1)
  txt.expanded<-expand.ellipsis(txt.ellipsis)
  model.syntax<-paste0(substring(model.syntax,1,match.start[i]-1),
                       txt.expanded,
                       substring(model.syntax,match.start[i]+match.length[i],nchar(model.syntax)))
}


if(grepl("...",model.syntax,fixed = TRUE)){
  #There is still ellipsis which does not match with the current pattern
  
  match.start<-gregexpr("...",model.syntax,fixed = TRUE)[[1]]
  txt<-substr(model.syntax,match.start[1]-20,match.start[1]+23)
  stop("gmlSEM error: Parser can not expand ellipsis in: \n",
       txt)
}

# break up in lines
model <- unlist( strsplit(model.syntax, "\n") )

# check for multi-line formulas: they contain no operator symbol
# but before we do that, we remove all strings between double quotes
# to avoid confusion with for example equal("f1=~x1") statements
# model.simple <- gsub("\\(.*\\)\\*", "MODIFIER*", model)
model.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", model)



start.idx <- grep(paste(operators, collapse = "|"), model.simple)

# check for empty start.idx: no operator found (new in 0.6-1)
if(length(start.idx) == 0L) {
  stop("ERROR: model does not contain gmlSEM syntax (no operators found)")
}


# check for non-empty string, without an operator in the first lines
# (new in 0.6-1)
if(start.idx[1] > 1L) {
  # two possibilities:
  # - we have an empty line (ok)
  # - the element contains no operator (warn!)
  for(el in 1:(start.idx[1] - 1L)) {
    # not empty?
    if(nchar(model.simple[el]) > 0L) {
      warning("gmlSEM WARNING: no operator found in this syntax line: ", model.simple[el], "\n", "                  This syntax line will be ignored!")
    }
  }
}

end.idx <- c( start.idx[-1]-1, length(model) )
model.orig    <- model
model <- character( length(start.idx) )
for(i in 1:length(start.idx)) {
  model[i] <- paste(model.orig[start.idx[i]:end.idx[i]], collapse="")
}

# ok, in all remaining lines, we should have an operator outside the ""
model.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", model)
idx.wrong <- which(!grepl(paste(operators, collapse = "|"),
                          model.simple))
#idx.wrong <- which(!grepl("[~=<>:|%]", model.simple))
if(length(idx.wrong) > 0) {
  cat("gmlSEM: missing operator in formula(s):\n")
  print(model[idx.wrong])
  stop("ERROR: syntax error in gmlSEM model syntax")
}

# but perhaps we have a '+' as the first character?
idx.wrong <- which(grepl("^\\+", model))
if(length(idx.wrong) > 0) {
  cat("gmlSEM: some formula(s) start with a plus (+) sign:\n")
  print(model[idx.wrong])
  stop("ERROR: syntax error in gmlSEM model syntax")
}


#Sort model by operator rank
ind<-sapply(model.simple, function(m)which(sapply(operators,function(o)grepl(o,m,fixed=TRUE)))[1])
ind<-sort(ind,index.return=T)$ix
model<-model[ind]
model.simple<-model.simple[ind]



for(i in 1:length(model)) {
  x <- model[i]
  if(debug) {
    cat("formula to parse:\n"); print(x); cat("\n")
  }
  
  # 1. which operator is used?
  line.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", x)
  
  op=""
  for(opp in operators)
    if(grepl(opp, line.simple, fixed=TRUE)){
      op<-opp
      break
    }
  if(op==""){
    stop("unknown operator in ", model[i])
  }
  bl=operators.blocks[op][[1]]
  x<-gsub(paste0("^(",bl,":)"),"",x,perl=TRUE)    

  # 2. split by operator (only the *first* occurence!)
  # check first if equal/label modifier has been used on the LEFT!
  if(substr(x,1,6) == "label(")
    stop("label modifier can not be used on the left-hand side of the operator")
  
  op.idx <- regexpr(op, x)
  
  lhs <- substr(x, 1L, op.idx-1L)
  
  # right-hand side string
  rhs <- substr(x, op.idx+attr(op.idx, "match.length"), nchar(x))
  
  # check if first character of rhs is '+'; if so, remove silently
  # (for those who copied multiline R input from a website/pdf)
  if(substr(rhs, 1, 1) == "+") {
    rhs <- substr(rhs, 2, nchar(rhs))
  }
  
  
  group=quote({
    #Look for group
    #Group is the first operator parser looks for,
    #gmlSEM expect zero or one group: block
    if(GROUP_OP)
      stop("gmlSEM Parser Eror: At most one group: statement is expected; Either remove extra blocks, or combine them into one factor.")
    
    if(grepl("(",rhs,fixed = TRUE)){
      
      i1<-gregexpr("(",rhs,fixed = TRUE)[[1]][1]
      g<-substring(rhs,1,i1-1)
      group    <- g
      
      txt<-substring(rhs,i1)
      if(!grepl("\\((?>[\\w\\.]+=[\\w\\.]+(?:,[\\w\\.]+=[\\w\\.]+)*|[\\w\\.]+(?:,[\\w\\.]+=[\\w\\.]+)*)\\)",rhs,perl=TRUE)){
        stop("gmlSEM Parser Error: The group: block is not of correct format.\n",
             "Template is\n group: groupLabel(label1=level1,label2=level2,...) or\n",
             "group: groupLabel(varname,label1=level1,label2=level2,...) or\n",
             "group: groupLabel(varname)")
      }else{
        res<-gregexpr("(?J)\\((?>(?'labels'[\\w\\.]+=[\\w\\.]+(?:,[\\w\\.]+=[\\w\\.]+)*)|(?'alias'[\\w\\.]+)(?:,(?'labels'[\\w\\.]+=[\\w\\.]+(?:,[\\w\\.]+=[\\w\\.]+)*))?)\\)",txt,perl=TRUE)[[1]]
        res<-captured.groups.list(txt,res)
        if(res$alias!="")
          add.alias(res$alias,g)
        txt<-strsplit(res$labels,",")[[1]]
        lbls<-sapply(txt, function(x){strsplit(x,"=",fixed=TRUE)[[1]]})
        lbls<-as.data.frame(t(lbls))
        rownames(lbls)<-NULL
        colnames(lbls)<-c("label","level")
        attr(group,"labels")<-lbls
        
        
        
      }
      
    }else{
      group    <- rhs  
      attr(group,"labels")<-NULL
    }
    
    
    GROUP_OP <- TRUE
  })
  
  alias=quote({
    add.alias(lhs,rhs)
  })
  
  
  
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
  
  level=quote({
    
    LEVEL_OP <- TRUE
    
    rhss <- strsplit(rhs,"<<")
    for(i1 in seq_len(rhss)){
      rhs<-rhss[[i1]]
      seps<-sep.par(rhs)
      rhs<-seps$var.out
      rhs2<-seps$var.in
    
      add.alias(rhs2,rhs1)
      add.levels(rhs)
      set.levels.matrix(lhs,rhs)<-1 
      lhs<-rhs
    }
    
  })
  
  vary=quote({

    
    LEVEL_OP <- TRUE
    
    seps<-sep.par(rhs)
    rhs<-seps$var.out
    rhs2<-seps$var.in
    
    add.alias(rhs2,rhs)
    add.levels(rhs)
    
    set.vary(lhs,rhs)
    
    
    
  })
  
  constraint=quote({
    
  })
  
  
  
  size=quote({
    
    rhs<-switch.space(rhs)
    pattern<-"(?J)(?'phrase'\\s*(?'N'[^\\s]*) *(?'level'[\\w\\.]+)(?=\\b) *(?:((?'keyword'per|in each|each) *(?'plevel'[\\w\\.]+)|(?'keyword'in) *(?'plevel'[\\w\\.]+) *(?'plevelid'\\d+)))?)"
    res<-gregexpr(pattern,rhs,ignore.case = TRUE, perl = TRUE)[[1]]
    d<-captured.groups.list(rhs,res)
    
    
  })
  
  family=quote({})
  regression=quote({})
  measurement=quote({})
  formative=quote({})
  covariance=quote({})
  monitor=quote({})
  unknown=quote({
    #Other, maybe a problem in syntax
  })
  
  task<-list(
    constraint=constraint,
    alias=alias,
    vary=vary,
    level=level,
    family=family,
    regression=regression,
    measurement=measurement,
    formative=formative,
    covariance=covariance,
    monitor=monitor,
    group=group,
    size=size,
    unknown=unknown)

  eval(task[[bl]])
  
  
}

if(LEVEL_OP){
  #Drop first row and column which are NA
  levels.matrix <- levels.matrix[-1,-1]
}else{
  levels.matrix <- matrix(nrow = 0,ncol = 0)
}

# Look for any inconsistency in level structure
# Levels are consistent if and only if 
#  any submatrix of levels.matrix contains a zero row
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

# Constructing Level structure
  .zeros<-sapply(1:nrow(levels.matrix), function(k){sum(levels.matrix[k,])==0})
  .levels<- rownames(levels.matrix)[which(.zeros)]
  for(lev in .levels){
    levels[[lev]]<-getSubLevels(lev)
  }




