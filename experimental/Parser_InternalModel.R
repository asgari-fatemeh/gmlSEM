
parseModel<-function(model){
  
}

GeneralPat="(?xms)
#Defining blocks
(?\'var\'[a-zA-Z\\._][\\w\\._]*){0}
(?\'reservedmodfunnames\'label|equal|start|c|(?=\\()){0}
(?\'balp\'\\((?:[^()\\\"\\\']*|\\\"[^\\\"]*\\\"|\\\'[^\\\']*\\\'|\\g\'balp\')*\\)){0}   #Balanced parantehsis 
(?\'eval\'\\b(?!\\g\'reservedmodfunnames\'\\b)\\g\'var\'\\g\'balp\'){0}   #functions needed evaluations 
(?\'modformat\'[\\-\\+]?\\d+|[a-zA-z\\._][\\w\\._]*|\\\"[^\\\"]*\\\"|\\\'[^\\\']*\\\'|\\g\'eval\'){0}
(?\'modfunarg\'(?:(?\'par\'\\g\'var\')=)?(?\'val\'\\g\'modformat\'|)){0}
(?\'modfunargs\'\\g\'modfunarg\'(,\\g\'modfunarg\')*){0}
(?\'modfun\'\\g\'reservedmodfunnames\'\\((?:\\g\'modfunargs\'|c?\\(\\g\'modfunargs\'\\))\\)){0}
(?\'modif\'\\g\'modformat\'|\\g\'modfun\'){0}
(?\'modifs\'\\g\'modif\'(?:\\*\\g\'modif\')*){0}
(?\'lag\'[\\+\\-]?\\d+){0}
(?\'varlag\'\\((?:\\g\'lag\'|\\g\'varlag\')*\\)){0}
(?\'varfunwithlag\'\\((?:\\g\'var\'\\g\'varlag\'?|\\g\'varfunwithlag\')*\\)){0}
(?\'variable\'1|\\g\'var\'\\g\'varfunwithlag\'|\\g\'var\'\\g\'varlag\'?){0}
(?\'single\'(?:\\-1|0|(?\'modifiers\'\\g\'modifs\'\\*)*\\g\'variable\'(?\'varlat\'@\\g\'var\')?)(?!\\*)){0}
(?\'singlep\'\\g\'single\'|\\(\\g\'singlep\'\\)){0} #extra paranthesis around terms
(?\'sep\'\\+\\n?|(?=\\-)){0}
(?\'multiple\'\\g\'singlep\'(?:\\g\'sep\'\\g\'singlep\')*){0}
(?\'multiplep\'\\g\'multiple\'|\\(\\g\'multiplep\'\\)){0}
(?\'multiplep2\'\\g\'multiplep\'(?:\\g\'sep\'\\g\'multiplep\')*|\\(\\g\'multiplep2\'\\)){0}
(?\'randomp\'\\((?:\\g\'multiplep2\'\\|(?\'level\'\\g\'var\')|\\(\\g\'randomp\'\\))\\)){0}#extra paranthesis around terms
(?\'effect\'\\g\'multiplep\'|\\g\'randomp\'){0}
(?\'multipleps\'\\g\'effect\'(?:\\g\'sep\'\\g\'effect\')*|\\(\\g\'multipleps\'\\)){0}
#(?\'multipleEffects\'\\g\'effect\'(?:\\g\'sep\'\\g\'effect\')*){0}
(?\'multipleEffects\'\\g\'multipleps\'(?:\\g\'sep\'\\g\'multipleps\')*|\\(\\g\'multipleEffects\'\\)){0}
(?Jx)\\g\'multipleEffects\'"

effect.groups.pattern="(?xms)
#Defining blocks
(?\'var\'[a-zA-Z\\._][\\w\\._]*){0}
(?\'reservedmodfunnames\'label|equal|start|c|(?=\\()){0}
(?\'balp\'\\((?:[^()\\\"\\\']*|\\\"[^\\\"]*\\\"|\\\'[^\\\']*\\\'|\\g\'balp\')*\\)){0}   #Balanced parantehsis 
(?\'eval\'\\b(?!\\g\'reservedmodfunnames\'\\b)\\g\'var\'\\g\'balp\'){0}   #functions needed evaluations 
(?\'modformat\'[\\-\\+]?\\d+|[a-zA-z\\._][\\w\\._]*|\\\"[^\\\"]*\\\"|\\\'[^\\\']*\\\'|\\g\'eval\'){0}
(?\'modfunarg\'(?:(?\'par\'\\g\'var\')=)?(?\'val\'\\g\'modformat\'|)){0}
(?\'modfunargs\'\\g\'modfunarg\'(,\\g\'modfunarg\')*){0}
(?\'modfun\'\\g\'reservedmodfunnames\'\\((?:\\g\'modfunargs\'|c?\\(\\g\'modfunargs\'\\))\\)){0}
(?\'modif\'\\g\'modformat\'|\\g\'modfun\'){0}
(?\'modifs\'\\g\'modif\'(?:\\*\\g\'modif\')*){0}
(?\'lag\'[\\+\\-]?\\d+){0}
(?\'varlag\'\\((?:\\g\'lag\'|\\g\'varlag\')*\\)){0}
(?\'varfunwithlag\'\\((?:\\g\'var\'\\g\'varlag\'?|\\g\'varfunwithlag\')*\\)){0}
(?\'variable\'1|\\g\'var\'\\g\'varfunwithlag\'|\\g\'var\'\\g\'varlag\'?){0}
(?\'single\'(?:\\-1|0|(?\'modifiers\'\\g\'modifs\'\\*)*\\g\'variable\'(?\'varlat\'@\\g\'var\')?)(?!\\*)){0}
(?\'singlep\'\\g\'single\'|\\(\\g\'singlep\'\\)){0} #extra paranthesis around terms
(?\'sep\'\\+\\n?|(?=\\-)){0}
(?\'multiple\'\\g\'singlep\'(?:\\g\'sep\'\\g\'singlep\')*){0}
(?\'multiplep\'\\g\'multiple\'|\\(\\g\'multiplep\'\\)){0}
(?\'multiplep2\'\\g\'multiplep\'(?:\\g\'sep\'\\g\'multiplep\')*|\\(\\g\'multiplep2\'\\)){0}
(?\'randomp\'\\((?:\\g\'multiplep2\'\\|(?\'level\'\\g\'var\')|\\(\\g\'randomp\'\\))\\)){0}#extra paranthesis around terms
(?\'effect\'\\g\'multiplep\'|\\g\'randomp\')"

#replace multiple spaces with one space character
model=gsub(" {1,}"," ",model,perl=TRUE)

#Seal the spaces between " " and ' '
spc=gregexpr("\"[^\"]*\"|'[^']*'",model,perl = TRUE)
if(spc[1]>=0){
  m.s=spc
  m.l=attr(spc,"match.length")
  for(i in rev(seq_along(spc))){
    model=paste0(substr(model,1,m.s[i]-1),
                 switch.space(substr(model,m.s[i],m.s[i]+m.l[i]-1),omit.spaces = TRUE),
                 substr(model,m.s[i]+m.l[i],nchar(model)))
  }
}

#Remove all other spaces
#The patterns do not match to the text with spaces
model=gsub(" ","",model,fixed=TRUE)

seps=paste0(model.operators,collapse = "|")
if(grepl(seps,model,perl = TRUE)){
  op=gregexpr(seps,model,perl=TRUE)[[1]]
  if(length(op)>1){
    stop("gmlSEM error: multiple model operators in one statement\n",model)
  }
  c.s=op
  c.l=attr(op,"match.length")
  op=substr(model,c.s,c.s+c.l-1)
  lhs=trim(substr(model,1,c.s-1))
  rhs=trim(substr(model,c.s+c.l,nchar(model)))
}else{
  op=""
  lhs=""
  rhs=model
}
lhs.match=rhs.match=TRUE
if(lhs!="")
  lhs.match=grepl(GeneralPat,lhs,perl = TRUE,ignore.case = TRUE)
rhs.match=grepl(GeneralPat,rhs,perl = TRUE,ignore.case = TRUE)
if(!lhs.match){
  stop("gmlSEM error: left hand side mismatch with expected pattern\n",model)
}
if(!rhs.match){
  stop("gmlSEM error: right hand side mismatch with expected pattern\n",model)
}

#match the pattern
model.data = data.frame(
  is.random = logical()   ,
  level     = character() ,
  varName   = character() ,
  latName   = character() ,  #If the term is a random effect, it can be named
  is.lagged = logical()   ,
  lag.integer = integer() ,
  scaled    = logical()   ,
  modified  = logical()   ,
  modifirs  = I(list())   #list(label=list(),start=list(),value=list()) #list of values or quotes
)

#Devide effects into bunch of fixed and random effect groups
effs=effect.groups=gregexpr(effect.groups.pattern,rhs,perl = TRUE,ignore.case = TRUE)[[1]]
for(i in 1:)







