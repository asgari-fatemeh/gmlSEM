# TODO: Family does not fully support vector-valued distributions

source("experimental/Parser_InternalMisc.R")
source("experimental/Parser_InternalEllipsis.R")
source("experimental/Parser_InternalFamilyClass.R")
source("experimental/Parser_InternalFamilyCollection.R")
debugSource("experimental/Parser_InternalObjects.R")

sample_syntax_file       = 'experimental/sample_gmlSEMScript.txt'
fileName = 'experimental/sample_gmlSEMScript.txt'
model.syntax      = readChar(sample_syntax_file, file.info(sample_syntax_file)$size)

#sample_syntax_extra_file = 'experimental/sample_gmlSEMScript_extra.txt'
#simulation.syntax = readChar(sample_syntax_extra_file, file.info(sample_syntax_extra_file)$size)
#model.syntax<-paste0(model.syntax,"\n",simulation.syntax)


#############################################################
##################### Parser Code ###########################
add.levels('1') 

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

#remove spaces from the begining of the lines
model.syntax <- gsub("(?<=\n|^)\\s*", "", model.syntax, perl = TRUE)

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



#Searching for random effects and adding them to the vars.level data.frame
#Before that we first replace the name assign operator i.e. [as .] with a simple syntax then we will expand three dots
capture.within <- gregexpr("\\[ *as +(?'alias'.*?) *\\]", model.syntax,perl = TRUE,ignore.case =TRUE)
match.start=capture.within[[1]]
match.length=attr(capture.within[[1]],"match.length")
capture.length=attr(capture.within[[1]],"capture.length")
capture.start<-attr(capture.within[[1]],"capture.start")
if(nrow(capture.start)>0){
  for(i in seq(nrow(capture.start),1,by=-1)){
    model.syntax<-paste0(substring(model.syntax,1,match.start[i]-1),"@",
                         substring(model.syntax,capture.start[i],capture.start[i]+capture.length[i]-1),
                         substring(model.syntax,match.start[i]+match.length[i]))
  }
}

# finding regression and expanding three dots
# This pattern supports multi-lines models
ran.list=list()
capture.within <- gregexpr("(?<=^|\n)[^\n][^~]*(?:<~|~|=~)(?'model'(?:[^~:\n<>]*[+\\-*/]\\s*)*[^~:\n<>]*(?=$|\n))", model.syntax,perl = TRUE,ignore.case =TRUE)
match.start=capture.within[[1]]
match.length=attr(capture.within[[1]],"match.length")
capture.length=attr(capture.within[[1]],"capture.length")
capture.start<-attr(capture.within[[1]],"capture.start")
if(nrow(capture.start)>0){
  for(i in seq(nrow(capture.start),1,by=-1)){
    
    lhs=trim(substring(model.syntax,match.start[i],capture.start[i]-1))
    mdl=substring(model.syntax,capture.start[i],capture.start[i]+capture.length[i]-1)
    mdl=trim(strsplit(mdl,"\n")[[1]])
    mdl=expand.dots.in.syntax(merge.models(mdl))
    
    model.syntax<-paste0(substring(model.syntax,1,match.start[i]-1),lhs,mdl,
                         substring(model.syntax,match.start[i]+match.length[i]))
    
    #Add random effects if any to ran.list to later adding them to vars.level 
    res=extract.modelterms.rhs(mdl,lhs)
    ran.list[[length(ran.list)+1]]=list(mdl=mdl,lhs=lhs)
    
    # res=extract.random.effects(mdl,lhs)
    # if(nrow(res)>0){
    #   ran.list[[length(ran.list)+1]]=list(mdl=mdl,lhs=lhs)
    # }
    
  }
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
    LEVEL_OP <- TRUE
    
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
        lev=trim(levv$var)
        lev.par=trim(levv$varn)
      }
      
      for(j in 1:length(txtm)){
        if(grepl("(",txtm[j],fixed=TRUE)){
          levv=captured.groups(txtm[j]," *(?'var'.*) *\\( *(?'varn'.*) *\\)")
          txtm[j]=trim(levv$var)
          txtm.par[j]=trim(levv$varn)
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
      
      add.alias(lev.par,lev,"OV")  #Level is an observed variable
      add.levels(lev)
      
      add.alias(txtm.par,txtm)
      add.vars.to.level(txtm,lev)
    }
    
  }  
}


#'within' keyword in level: block
#'template allow multi line definition of level structure
while(TRUE){
  capture.within <- gregexpr("(?<=>>L>>:)((?:\\s*,? *[\\w\\.][\\w\\._]* *(?:\\( *[\\w\\.][\\w\\._]* *\\))?)+)(?=$|\n)", model.syntax,perl = TRUE,ignore.case =TRUE)
  if(capture.within[[1]][1]==-1)
    break
  
  capture.start<-attr(capture.within[[1]],"capture.start")
  capture.length<-attr(capture.within[[1]],"capture.length")
  match.length<-attr(capture.within[[1]],"match.length")
  match.start<-capture.within[[1]]
  if(length(match.start)>0){
    
    LEVEL_OP <- TRUE
    
    for(i in seq(length(match.start),1,by=-1)){
      level_text_ml<-trim(substring(model.syntax,capture.start[i],capture.start[i]+capture.length[i]))
      level_text_ml=strsplit(level_text_ml,"\n")[[1]]  #Multi line statement
      for(m in seq_along(level_text_ml)){
        level_text = trim(level_text_ml[m])
        level_text = trim(gsub(" within ","<<",level_text,fixed=TRUE))
        
        rhss = strsplit(level_text,"<<")[[1]]
        
        for(i1 in seq_along(rhss)){
          levs.par=levs=levs.in=trim(strsplit(rhss[i1],",")[[1]])
          for(j1 in seq_along(levs.in)){
            rhs=trim(levs.in[j1])
            if(grepl("(",rhs,fixed=TRUE)){
              levv=captured.groups(rhs," *(?'var'.*) *\\( *(?'varn'.*) *\\)")
              levs[j1]=trim(levv$var)
              levs.par[j1]=trim(levv$varn)
            }
            add.alias(levs.par,levs,"OV")  #Level is an observed variable
            add.levels(levs)  
          }
          
          if(i1>1)
            set.levels.matrix(levs.old,levs)
          levs.old=levs
        }
      }
      model.syntax<-paste0(substring(model.syntax,1,capture.start[i]-1),"<<",
                           level_text,
                           substring(model.syntax,match.start[i]+match.length[i]))
      
    }
    
  }  
}

#########################################################
################# Processing Levels #####################

  # Look for any inconsistency in level structure
  # Levels are consistent if and only if 
  #  any submatrix of levels.matrix contains a zero row
  validate.level.hierarchy()

  #Now that the alias list is completely loaded we can add random effects to vars.level
  #We postponed this procedure for using aliases in naming conventions
  if(length(ran.list)>0)
    for(i in 1:length(ran.list)){
      res=extract.random.effects(ran.list[[i]]$mdl,ran.list[[i]]$lhs)
      if(nrow(res)>0){
        for(j in 1:nrow(res)){
          lev.lhs=get.level.lhs(lhs)
          if(!is.na(lev.lhs) && any(sapply(unique(lev.lhs),function(lev)!levels.are.consistent(lev,res$level[j]))))
            stop("\ngmlSEM error: levels are inconsistent at the regression model '",lhs,"'\n'",
                 lhs,"' varies at the level '",paste0(unique(lev.lhs),collapse=","),"' while the random effect '",
                 res$var,"' varies at the level '",res$level[j],"'.\n")
          if(res$label[j]!=""){
            add.alias(res$dummy.label[j],res$label[j])
            add.vars.to.level(res$label[j],res$level[j])
          }else{
            add.vars.to.level(res$dummy.label[j],res$level[j])
          }
        }
      }
    }
  
  ######## Every level has its own conditional covariance structure 
  covs.list=list()
  for(nm in names(levels)){
    covs.list[[nm]]=matrix(NA,0,0)
  }

###########################

#At this step we which variables vary at which level. 
#Now, we can specify copula (if any) to the variables vary at the same level

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
    cop=""
    fam=""
    if(capture.start[i,3]>0){
      cop=trim(substring(model.syntax,capture.start[i,3],capture.start[i,3]+capture.length[i,3]-1))
      txt2<-paste0(">>F>>",cop)
      cop=srt2gmlSEMfamily(cop)
      if(cop$name!="copula"){
        stop("\ngmlSEM error: Copula was expected at line:\n",
             substring(model.syntax,1,match.start[i]-1,match.start[i]+match.length[i]))
      }
    }
      
      fmn=substring(model.syntax,capture.start[i,2],capture.start[i,2]+capture.length[i,2]-1)
      txt1<-paste0(">>F>>",fmn)
      fmn=srt2gmlSEMfamily(fmn)
      if(fmn$name=="copula"){
        if(txt2!=""){
          #Two copula statement
          stop("\ngmlSEM error: Two copula statement is read at line:\n",
               substring(model.syntax,1,match.start[i]-1,match.start[i]+match.length[i]))
        }
        cop=fmn
        fmn=NULL
      }
      
      if(!is.null(fmn)){
       n1=fmn$dim  
       n2=fmn$dim.latent
      }
      
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




parsedData[["vars.level"]] = vars.level

