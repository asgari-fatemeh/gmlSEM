

patterns.definition.blocks="(?xmsni)
#Support simple expressions
#support indexes in modifiers
# beta[s(,t)?]|beta[group[t]]
#Defining blocks
(?'var'[a-zA-Z\\._][\\w\\._]*){0} 
(?'reservedmodfunnames'(?>label|equal|start|c|(?=\\())){0}
(?'balp'\\((?>\\g'balp'|\\\"[^\\\"]*\\\"|\\'[^\\']*\\'|[^()\\\"\\']*)*\\)){0}   #Balanced parantehsis 
(?'eval'\\b(?!\\g'reservedmodfunnames'\\b)\\g'var'\\g'balp'){0}   #functions needed evaluations 
(?'modformat'(?>\\\"[^\\\"]*\\\"|\\'[^\\']*\\'|[\\-\\+]?\\d+|\\g'varwithlag'|[a-zA-z\\._][\\w\\._]*)|\\g'eval'){0}
(?'modfunarg'(?:(?'par'\\g'var')=)?(?'val'\\g'modformat')){0}
(?'modfunargs'\\g'modfunarg'(,\\g'modfunarg')*){0}
(?'modfun'\\g'reservedmodfunnames'\\((?:\\g'modfunargs'|c?\\(\\g'modfunargs'\\))\\)){0}
(?'modif'\\g'modformat'|\\g'modfun'){0}
(?'lag'(?:\\g'simpleExppVar'|\\g'varwithlag')){0}
(?'varwithlag'\\g'var'(?:\\[(?:\\g'lag'(?:,\\g'lag')*)?\\])?){0}
(?'varfunwithlag'(?>(?'varfun'\\g'var')\\(\\g'varwithlag'(?'varfunargs',\\g'modfunarg')*\\)|\\g'varwithlag')(?!\\*)){0}
(?'variable'(?>1|\\g'varfunwithlag'|\\(\\g'simpleExpression'\\))){0}
(?'simpleExp'(?>\\g'varwithlag'|[\\+\\-]?\\d+)([\\+\\-\\*\\/\\^](?>\\g'varwithlag'|\\d+))*){0}
(?'simpleExpVar'(?>\\g'varwithlag'|[\\+\\-]?\\d+)([\\+\\-\\*\\/\\^](?>\\g'varwithlag'|\\d+))*){0}
(?'simpleExpp'(?>\\(\\g'simpleExpp'\\)|\\g'simpleExp'|\\(\\g'simpleExpression'\\))){0}
(?'simpleExppVar'(?>\\(\\g'simpleExppVar'\\)|\\g'simpleExpVar')){0}
(?'simpleExpression'\\g'simpleExpp'([\\+\\-\\*\\/\\^]\\g'simpleExpp')*){0}
(?'singleEffect'(?:-1|0|(?:(?'modifiers'\\g'modif'\\*)*(?'vardesc'\\g'variable')(?:@(?'varlat'\\g'varwithlag'))?))(?!\\*)){0}
(?'singlep'(?:\\g'singleEffect'|\\(\\g'singlep'\\))){0} #extra paranthesis around terms
(?'sep'\\+\\n?|(?=\\-)){0}
(?'multiple'\\g'singlep'(?:\\g'sep'\\g'singlep')*){0}
(?'multiplep'(?>\\g'multiple'|\\(\\g'multiplep'\\))){0}
(?'multiplep2'\\g'multiplep'(?:\\g'sep'\\g'multiplep')*|\\(\\g'multiplep2'\\)){0}
(?'randomp'\\((?:\\g'multiplep2'\\|(?'level'\\g'varwithlag')|\\(\\g'randomp'\\))\\)){0}#extra paranthesis around terms
(?'effectGroup'\\g'multiplep'|\\g'randomp'){0}
(?'multipleEGs'\\g'effectGroup'(?:\\g'sep'\\g'effectGroup')*|\\(\\g'multipleEGs'\\)){0}
(?'generalPattern'\\g'multipleEGs'(?:\\g'sep'\\g'multipleEGs')*|\\(\\g'generalPattern'\\)){0}
\\K
\\g'generalPattern'"

# General pattern to test the model
general.pattern                = patterns.definition.blocks %+% "^\\g'generalPattern'$"

# Separating group of effects (random terms and fixed-terms)
effect.groups.pattern          = patterns.definition.blocks %+% "\\g'effectGroup'"

# Test if a group of terms are random effects
level.format                   = patterns.definition.blocks %+% "(?xms)(?'ranparbal'\\((?>(?'groupEffects'[^\\|]*)\\|\\b(?'grouplevel'\\g'var')\\b(?:\\((?'grouplevellag'\\g'lag')\\))?\\))|\\g'ranparbal'\\))$"

single.effect.format           = patterns.definition.blocks %+% "\\g'singleEffect'"
singe.effect.elements.format   = patterns.definition.blocks %+% "(?>(?'intercept'(?>-1|0))|(?'elemmodifs'\\g'modifs'\\*)*(?'elemvar'\\g'variable')(?:@(?'elemlat'\\g'var'))?)(?!\\*)"
fun.args.format                = patterns.definition.blocks %+% "(?:,(?:(?'funpar'\\g'var')=)?(?'funarg'\\g'modformat'))"
modifer.split.format           = patterns.definition.blocks %+% "(?:(\\g'modif')\\*)"
single.modifier.format         = patterns.definition.blocks %+% "(?>(?'funname'\\g'reservedmodfunnames')\\((?'copt'c?\\()?(?'funargs'\\g'modfunargs')(?('copt')\\))\\)|(?'singlemod'\\g'modformat'))"
vardesc.format                 = patterns.definition.blocks %+% "(?'variabeldef'(?>1|\\g'varfunwithlag'))"

## Parsing the model
# 0) validate the model against the pat=/general.pattern/
# 1) splitting group of effects into random and fixed terms using pat=/effect.groups.pattern/
# 2) For the random terms (includes a vertical bar | ) -> identify the level using pat=/level.format/ and follow 3+ for effect group
# 3) group of Effects will split using pat=/single.effect.format/
# 4) signle terms will decomposed into into elemmodifs*elemvar@elemlat using pat=/singe.effect.elements.format/
# 5) vardesc will be decomposed into fun(x(lag),args) using pat=/vardesc.format/
# 6) modifiers will split using pat=/modifer.split.format/
# 7) single modifier will decomposed int reservedfun(args) using pat=/single.modifier.format/


parseModel<-function(model){
  #replace multiple spaces with one space character
  model = gsub(" {1,}"," ",model,perl=TRUE)
  
  #Seal the spaces between " " and ' '
  model = switch.space.within.quotes(model,seal.spaces = TRUE)
  
  #Remove all other spaces
  #The patterns do not match a text with spaces
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
    lhs.match=grepl(general.pattern,lhs,perl = TRUE,ignore.case = TRUE)
  rhs.match=grepl(general.pattern,rhs,perl = TRUE,ignore.case = TRUE)
  if(!lhs.match){
    stop("gmlSEM error: left hand side mismatch with expected pattern\n",model)
  }
  if(!rhs.match){
    stop("gmlSEM error: right hand side mismatch with expected pattern\n",model)
  }
  
  model=gsub("\n","",model,fixed = TRUE)
  #removing duplicate plus signs
  model=gsub("\\+{1,}","\\+",model,perl=TRUE)
  
  #match the pattern
  model.data = new.effects.data.frame()
  
  
  #Divide effects into bunch of fixed and random effect groups
  effect.groups=gregexpr(effect.groups.pattern,rhs,perl = TRUE,ignore.case = TRUE)[[1]]
  m.s=effect.groups
  m.l=attr(effect.groups,"match.length")
  
  if(effect.groups[1]<0)
    stop("gmlSEM error: no terms on the right hand side of model\n",model)
  for(i in seq_along(m.s)){
    effect.group=substr(rhs,m.s[i],m.s[i]+m.l[i]-1)
    level=NA
    level.lag=NA
    
    if(grepl("|",effect.group,fixed = TRUE)){
      #It is a random effect
      lev=gregexpr(level.format,effect.group,perl = TRUE)[[1]]
      lev=captured.groups.list(effect.group,lev)
      level=lev$grouplevel
      level.lag=lev$grouplevellag
      level.lag=ifelse(is.null(level.lag),0,level.lag)
      if(level=="")
        gmlSEMerror("level is not specified at \n",effect.group)
      
      effect.group=lev$groupEffects
    }
    
    is.random=ifelse(is.na(level),FALSE,TRUE)
    
    single.effects = gregexpr(single.effect.format,effect.group,perl=TRUE,ignore.case = TRUE)[[1]]
    m.s.s=single.effects
    m.s.l=attr(single.effects,"match.length")
    
    for(j in seq_along(single.effects)){
      single.effect=substr(effect.group,m.s.s[j],m.s.s[j]+m.s.l[j]-1)
      #extract data for single.effect
      single.effect.data=getPartialModel.for.singleEffect(single.effect,is.random = is.random,level = level,level.lag = level.lag)
      model.data[nrow(model.data)+1,]=single.effect.data
    }
  }
  
  model.data=combine.model.data.rows(model.data)
 
  list(lhs=strsplit(lhs,"+",fixed = TRUE)[[1]],
       op=op,
       rhs.raw=rhs,
       model.data=model.data)
}




