# TODO: Family does not fully support vector-valued distributions

source("experimental/Parser_InternalMisc.R")
source("experimental/Parser_InternalEllipsis.R")
source("experimental/Parser_InternalFamilyClass.R")
source("experimental/Parser_InternalFamilyCollection.R")

sample_syntax_file       = 'experimental/sample_gmlSEMScript.txt'
fileName = 'experimental/sample_gmlSEMScript.txt'
model.syntax      = readChar(sample_syntax_file, file.info(sample_syntax_file)$size)

#sample_syntax_extra_file = 'experimental/sample_gmlSEMScript_extra.txt'
#simulation.syntax = readChar(sample_syntax_extra_file, file.info(sample_syntax_extra_file)$size)
#model.syntax<-paste0(model.syntax,"\n",simulation.syntax)

##Arguments
#internal machinery
data.is.provided = FALSE

#############################################################
##################### Parser Code ###########################

# check for empty syntax
if(length(model.syntax) == 0) {
  stop("gmlSEM ERROR: empty model syntax")
}

# remove comments prior to split:
# match from comment character to newline, but don't eliminate newline
model.syntax <- gsub("[#!].*(?=\n)", "", model.syntax, perl = TRUE)

# replace semicolons with newlines prior to split
model.syntax <- gsub(";", "\n", model.syntax, fixed = TRUE)
model.syntax <- gsub("\r\n", "\n", model.syntax, fixed = TRUE)
model.syntax <- gsub(pattern = "\u02dc", replacement = "~", model.syntax)

model.syntax <- gsub("(?<=^|\r\n|\n) *group:", ">>G>>", model.syntax, perl = TRUE,ignore.case =TRUE)


# Sugar syntax for simulation
# Sample size can be set within size: blocks
pattern1<-"(?J)(?:(?<=^|\r\n|\n) *size:)(?'phrase'\\s*(?'N'[^\\s]*) *(?'level'[\\w\\.]+)(?=\\b) *(?:((?'keyword'per|in each|each) *(?'plevel'[\\w\\.]+)|(?'keyword'in) *(?'plevel'[\\w\\.]+) *(?'plevelid'\\d+)))?)(?:[,;]\\g'phrase')*"
pattern2<-"(?J)(?'phrase'\\s*(?'N'[^\\s]*) *(?'level'[\\w\\.]+)(?=\\b) *(?:((?'keyword'per|in each|each) *(?'plevel'[\\w\\.]+)|(?'keyword'in) *(?'plevel'[\\w\\.]+) *(?'plevelid'\\d+)))?)[,;]?"
res<-gregexpr(pattern1,model.syntax,ignore.case = TRUE, perl = TRUE)[[1]]
match.length<-attr(res,"match.length")
for(i in rev(seq_along(res))){
  txt<-substr(model.syntax,res[i],res[i]+match.length[i]-1)
  txt<-gsub("size:","",txt,ignore.case = TRUE)
  txt<-gsub("\n","",txt,fixed = TRUE)
  if(grepl(",\\s*\\.\\.\\.\\s*,",txt,perl=TRUE)){
    ###### 2Do: 
    ###### Postpone expanding ellipsis until variables are resolved
    txt<-gsub(",\\s*\\.\\.\\.\\s*,",",...,",txt,perl=TRUE)
    txt<-trim(txt)
    txt<-switch.space(txt)
    txt<-expand.ellipsis(txt)
  }else{
    txt<-switch.space(txt)  
  }
  txt<-strsplit(txt,",",fixed=TRUE)
  txt<-sapply(txt, function(x)paste0(">>S>>",x))
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
capture.within <- gregexpr("(?:^|\r\n|\n) *?((?>\\b|\\.)[\\w\\.\\+\\, ]+?\\b) *( as ) *(?! *,)(.*)", model.syntax, perl = TRUE,ignore.case =TRUE)
match.start=capture.within[[1]]
match.length=attr(capture.within[[1]],"match.length")
capture.length=attr(capture.within[[1]],"capture.length")
capture.length=as.matrix(capture.length,ncol=2)
capture.start<-attr(capture.within[[1]],"capture.start")
capture.start=as.matrix(capture.start,ncol=2)
if(nrow(capture.start)>0){
  for(i in seq(nrow(capture.start),1,by=-1))
    model.syntax<-paste0(substring(model.syntax,1,capture.start[i,1]-1),"alias:",
                         expand.ellipsis(substring(model.syntax,capture.start[i,1],capture.start[i,1]+capture.length[i,1])),"<>",
                         expand.ellipsis(switch.space(trim(substring(model.syntax,capture.start[i,3],capture.start[i,3]+capture.length[i,3]-1)))),
                         substring(model.syntax,match.start[i]+match.length[i]))
}


#'within' keyword in level: command
while(TRUE){
  capture.within <- gregexpr("(?<=level:)([^<]+?)(?:$|\r\n|\n)", model.syntax,perl = TRUE,ignore.case =TRUE)
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


#replacing '(vary|varies|varying) at level' with <<~
capture.within <- gregexpr("(?:^|\r\n|\n)\\s*?[\\w\\.\\,\\+ ]+\\s*? ((vary|varies|varying) at level) ", model.syntax, perl=TRUE,ignore.case =TRUE)
if(capture.within[[1]][1]>0){
  capture.start=attr(capture.within[[1]],"capture.start")
  capture.length=attr(capture.within[[1]],"capture.length")
  
  capture.start=capture.start[,1]
  capture.length=capture.length[,1]
  
  if(length(capture.start)>0){
    for(i in seq(length(capture.start),1,by=-1))
      model.syntax<-paste0(substring(model.syntax,1,capture.start[i]-1),"<<~",
                           substring(model.syntax,capture.start[i]+capture.length[i]))
  }  
}


#marking distributions with tag >>F>>
#template allow breaklines after the list of variables
pattern='(?<=\\n|^) *family:\\s*([\\w\\._\\,\\+\\s()\\*]*?)\\s*([\\w\\.][\\w\\._]*?(?:\\((?:[^()]|(?2))*\\)) *?)\\s*([\\w\\.][\\w\\._]*?(?:\\((?:[^()]|(?3))*\\)) *?)?(?=\\n|$)'
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
      
    model.syntax<-paste0(substring(model.syntax,1,match.start[i]-1),"family:",txt,"\n",
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


GROUP_OP <- FALSE
LEVEL_OP <- FALSE
group<-""

#Parser Output
vars.dictionary <- data.frame(
                       name         = character() ,
                       alias        = I(list())   ,
                       nameInData   = character() ,
                       nameToPrint  = character() ,
                       appeared.as.fa = logical() , #Appeared as a factor in a measrement model
                       appeared.as.re = logical() , #Appeared as a random effect in a regression model
                       latent       = logical()   , #NA if needs.data.scan
                       response     = logical()   ,
                       family       = I(list())   ,
                       vary         = character() ,
                       heter        = logical()   ,
                       heter.params = I(list())   ,
                       heter.vars   = I(list())   ,
                       reg.params   = I(list())   ,
                       reg.vars     = I(list())   ,
                       reg.smooth   = logical()   ,
                       depends.on   = I(list())
                       )
all.alias       <- data.frame(
  lhs=character(),
  rhs=character(),
  var.name=character(),
  role=character(), #PAR, LV, OV
  all.forms=I(list())
)

covs            <- list()
levels          <- list()
levels.matrix   <- matrix()
rownames(levels.matrix) <- colnames(levels.matrix) <- NA


#gmlSEM parser environment to resolve parameter values
env=new.env()

group.specs     <- list(
  "" = list(
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
  if(lev=="")
    return()
  
  rhs2<-get.alias.rhs(lev)
  if(lev=="1"|rhs2=="1")
    return()  #Do not add the base level to the matrix
  
  if(! rhs2 %in%  rownames(levels.matrix)){
    levels.matrix<-rbind(levels.matrix,0)
    levels.matrix<-cbind(levels.matrix,0)
    rownames(levels.matrix)<-c(rownames(levels.matrix),rhs2)
    colnames(levels.matrix)<-c(colnames(levels.matrix),rhs2)
  }
  
}
set.levels.matrix<-function(lhs,rhs){
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

get.alias.ind<-function(x){
    inds<-apply(all.alias,1, function(x) lhs%in% x[['all.forms']])
  which(inds)
}
#Reconsider relabling ability of unlabled symbols a=a -> a=b
add.alias<-function(lhs,rhs=NULL,role=NA){
  
  if(is.null(rhs)){
    inds<-get.alias.ind(lhs)
    
    if(nrow(inds)==0){
      ind<<-nrow(all.alias)+1
      all.alias[ind,]<<-NA
      all.alias[ind,1:4]<<-c(lhs,lhs,"",role)
      all.alias$all.forms[[ind]]<<-c(lhs)
    }
    
    return()
  }
  
  if(lhs==""||rhs=="")
    stop("glmSEM error: null alias can not be set for",lhs,rhs)
  
  ind<- unique(c(get.alias.ind(lhs),get.alias.ind(rhs)))

  if(length(ind)==0){
    ind<-nrow(all.alias)+1
    all.alias[ind,1:4]<<-c(lhs,rhs,"",role)
    all.alias$all.forms[[ind]]<<-c(lhs,rhs)
    
  }else if(length(ind)>1){
    i1<-get.alias.ind(lhs)
    i2<-get.alias.ind(rhs)
    
    a1<-all.alias[i1,]
    a2<-all.alias[i2,]
    
    if((a1[1]==a1[2])&(a2[1]==a2[2])){
      #remove the two rows and define the alias
      all.alias<<-all.alias[-c(i1,i2),]
      ind<-nrow(all.alias)+1
      all.alias[ind,1:4]<<-c(lhs,rhs,"",role)
      all.alias$all.forms[[ind]]<<-unique(c(lhs,rhs,a1$all.forms[[1]],a2$all.forms[[1]]))
    }else{
      b1<-ifelse(a1[1]==lhs,a1[2],a1[1])
      b2<-ifelse(a2[1]==rhs,a2[2],a2[1])
      stop("\ngmlSEM Parser Error: Multiple definition for alias:\n",
           ifelse(b1!=lhs,lhs%+%" is defined as an alias for "%+%b1%+%"\n",""),
           ifelse(b2!=rhs,rhs%+%" is defined as an alias for "%+%b2%+%"\n",""),
           "Thus, ",lhs," and ",rhs," cannot be defined as alias.")
    }
  }else{
    a<-all.alias[ind,]
    if((a[1]=="lhs")&(a[2]==rhs)){
      #Do nothing
      return()
    }else if(a[1]==a[2]){#Update the alias
      all.alias[ind,1:4]<<-c(lhs,rhs,"",role)
      all.alias$all.forms[ind]<-unique(c(lhs,rhs,all.alias$all.forms[ind]))
    }else{
      if(a[2]==lhs&a[1]==rhs){
        warning("\ngmlSEM Parser Error: Multiple definition for alias:\n",
                a[1]," is defined as alias for ",a[2], 
                " while ",a[2], "is defined on the right hand side.",
                '\nNew definition for "',lhs,"' as '",rhs,"' with '",rhs,"' on the right hand side will be ignored.")
      }else{
        stop("\ngmlSEM Parser Error: Multiple definition for alias:\n",
             a[1]," is defined as an alias for ",a[2],"\n",
             "Thus, ",lhs," and ",rhs," cannot be defined as new aliases.")
      }
      
    }
  }
  
}

add.alias.new.forms<-function(al,newform){
  #Add alias if not exists
  al<-get.alias.lhs(al)
  ind<-get.alias.ind(al)
  all.alias$all.forms[[ind]]<<-unique(c(all.alias$all.forms[[ind]],newform))
}

get.alias.lhs<-function(lbl){
  ind<-get.alias.ind(lbl)
  
  if(length(ind)>0)
    return(all.alias[ind,1])
  
  if(lbl!="1")
    add.alias(lbl,lbl)
  
  lbl
}

get.alias.rhs<-function(lbl){
  ind<-get.alias.ind(lbl)
  
  if(length(ind)>0)
    return(all.alias[ind,2])
  
  if(lbl!="1")
    add.alias(lbl,lbl)
  
  lbl
}


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
        res<-captured.groups(txt,res)
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
  
  
  captured.groups<-function(txt,gpr){
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
  
  size=quote({
    
    rhs<-switch.space(rhs)
    pattern<-"(?J)(?'phrase'\\s*(?'N'[^\\s]*) *(?'level'[\\w\\.]+)(?=\\b) *(?:((?'keyword'per|in each|each) *(?'plevel'[\\w\\.]+)|(?'keyword'in) *(?'plevel'[\\w\\.]+) *(?'plevelid'\\d+)))?)"
    res<-gregexpr(pattern,rhs,ignore.case = TRUE, perl = TRUE)[[1]]
    d<-captured.groups(rhs,res)
    
    
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




