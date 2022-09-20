
## As a fundamental rule: Dist specific arguments maybe quoted for lazy loadings,
# The name must be aware of that, and to evaluate the params run time

setClass("support",slots=c(
  type    = "character"  ,
  is.integer = "numeric" ,
  support = "numeric"    ,   
  include.lhs ="numeric" ,
  include.rhs ="numeric" ,
  value   = "numeric"
))

setValidity("support", function(object){
  type=object@type
  support=object@support
  value=object@value
  is.integer=object@is.integer
  include.lhs=object@include.lhs
  include.rhs=object@include.rhs
  
  if(length(type)!=1 | length(support)==0 | length(value)>1 |
     length(is.integer)!=1|length(include.rhs)!=1|length(include.lhs)!=1)
    return("Error in parameter")
  if(type=="interval" & length(support)!=2)
    return("Parameter support for interval-valued parameter must be of length 2")
  if(length(value)>0){
    if(!is.na(value)){
    if((type=="discrete" & !value %in% support))
      return("Parameter value not in the parameter space")
    if(type=="interval" && (
       ((include.rhs && value>support[2]) || (!include.lhs & value>=support[2])) || 
       ((include.lhs && value<support[1]) || (!include.lhs && value<=support[1]))
          ))
      return("Parameter value not in the parameter space")
    if(type=="interval" & is.integer==1 & (floor(value)!=value))
      return("Parameter value not in the parameter space")
  }}
  return(TRUE)
})

setMethod("initialize", "support", 
          function(.Object, ...) {
            .Object <- callNextMethod()
            
            if(length(.Object@type) != 1 )
              stop("Type of parameter must be specified")
            
            if(length(.Object@support) ==0 )
              stop("Support of parameter must be specified")
            
            if(length(.Object@is.integer) != 1 )
              stop("Support of parameter must be specified")
            
            if(.Object@type=="interval" & length(.Object@support)!=2)
              stop("Parameter support for interval-valued parameter must be of length 2")
            
            
            if(length(.Object@include.lhs)==0)
              .Object@include.lhs=0
            
            if(length(.Object@include.rhs)==0)
              .Object@include.rhs=0
            
            if(.Object@type=="discrete" && all(sapply(.Object@support,function(x){floor(x)==x})))
              .Object@is.integer=1
            
            if(.Object@type=="discrete"){
              .Object@include.lhs=.Object@include.rhs=1
            }
              
            
            if(.Object@is.integer==1 && .Object@type=="interval"){
              if(.Object@support[1]!=-Inf)
                .Object@include.lhs=1
              if(.Object@support[2]!=Inf)
                .Object@include.rhs=1
            }
            
            if(length(.Object@value)==0)
              .Object@value=NA_real_
            
            .Object
          })


newSupport<-function(type=c("interval","discrete"),
               support=NULL,
               include.lhs=FALSE,include.rhs=FALSE,
               is.integer=FALSE){
  include.lhs=1*include.lhs
  include.rhs=1*include.rhs
  is.integer=1*is.integer
  type=match.arg(type)
  new("support",type=type,is.integer=is.integer,
      support=support,
      include.lhs=include.lhs,
      include.rhs=include.rhs)
}

validate.value.support<-function(dom,value){
  if(dom@is.integer==1){
    if(dom@support[1]!=-Inf)
      dom@include.lhs==1
    if(dom@support[length(dom@support)]!=Inf)
      dom@include.rhs==1
  }
    if((dom@type=="discrete" & !value %in% dom@support))
      return(FALSE)
    if(dom@type=="interval" && (
      ((dom@include.rhs && value>dom@support[2]) || (!dom@include.lhs & value>=dom@support[2])) || 
      ((dom@include.lhs && value<dom@support[1]) || (!dom@include.lhs && value<=dom@support[1]))
    ))
      return(FALSE)
    if(dom@type=="interval" & dom@is.integer==1 & (floor(value)!=value))
      return(FALSE)
  TRUE
}

ArgSymbol<-function(f){  #Extracting symbols usd in f
  f2=rlang::call_args(f)
  args=unlist(lapply(f2, function(x){if(is.call(x)){ArgSymbol(x)}else{x}}))
  args2=sapply(args, function(e){(is.name(e))})
  args[args2]  
}

#Support and params can be fixed objects of 'support' class,
# or functions with no argument returns an object of 'support' depends on other parameters 
# on.args.change is a function that will react to change in the name specific arguments
# This functions if of signature function(...){} which can read and write internal parameters and name specific parameters. The arguments which has changed will be passed to the argument
newFamily<-function(fname,
                    dim=1,          #support's dimension: number of elements in vector-valued random variables
                    support=NULL, 
                    params=NULL,
                    link.mean=NULL,
                    link.scale=NULL,
                    on.args.change=NULL,
                    args.mandatory=NULL, #name of arguments with no default values that must be supplemented in a gmlSEM syntax
                    ...){
  gmlSEMfamily.classname<-"gmlSEMfamily"
  
  if(is.null(on.args.change)){
    on.args.change=list()
  }else if(is.function(on.args.change)){
    on.args.change=list(on.args.change)
  }
  
  #Looking for latent process family declaerations, and put a constraint on them
  dim.latent=0   #dimension of latent generator process
  lst=list(...)
  inds=sapply(lst, function(l) gmlSEMfamily.classname %in% class(l))
  k=length(on.args.change)
  if(any(inds)){
    inds=which(inds)
    dim.latent=sum(sapply(inds, function(i)lst[[i]]$dim))
    for(i in seq_along(inds)){
      ind=inds[i]
      nm=names(lst)[ind]
      fn='function(%1$s){
            if(is.null(%1$s))
              return()
            if(%1$s$dim.latent>0)
              stop("gmlSEM error: the underlying generating process of the family %2$s cannot no be defined as a latent generating process itself.")
            }'
      fn=eval(str2lang(sprintf(fn,nm,paste0("'",fname[1],"'"))))
      
      k=k+1
      on.args.change[[k]]=fn
    }
  }
  
 
  
  nls.onargc=ifelse(is.null(on.args.change),list(),
                     ifelse(is.function(on.args.change),
                        list(names(formals(on.args.change))),
                        lapply(on.args.change, function(o)names(formals(o)))))
  
  if(is.null(on.args.change) || (length(nls.onargc)==1 && is.null(nls.onargc[[1]]))){
    nls.onargc=list()
  }
  if(is.function(on.args.change)){
    on.args.change=list(on.args.change)
  }
  if(length(nls.onargc)>0 && any(sapply(nls.onargc, function(o)length(o)==0)))
    stop("gmlSEM error: (developer warning!) on.args.change must have been defined with arguments intends to capture")
  
  
  registeredArg<-function(gmlSEMfamily.obj,args){
    if(is.list(args))
      args=names(args)
    
    nls.onargc=gmlSEMfamily.obj$nls.onargc
    
    if(length(nls.onargc)==0)
      return(numeric())
    
    res=sapply(nls.onargc, function(o){
      any(args%in%o)
    })
    
    which(res)
  }
  
  #Parameters which needs to be evaluated at run time.
  #These are parameters that fails to be evaluated at parsing time
  runtime.pars=c()
  
  #Family specific arguments are exposed
  args.exposed=names(list(...)) 
    #arguments truncated,censored,categorized, and inflated are always exposed
  args.exposed<-unique(c(args.exposed,"truncated","censored",
                         "categorized","inflated",args.exposed))
  
  if(is.null(link.mean)||length(link.mean)==0){
    args.exposed<-setdiff(args.exposed,c("link.mean","link"))
  }else{
    args.exposed<-unique(c("link.mean","link",args.exposed))
  }
  
  if(is.null(link.scale)||length(link.scale)==0){
    args.exposed<-setdiff(args.exposed,"link.scale")
  }else{
    args.exposed<-unique(c("link.scale",args.exposed))
  }
  
 
  interval.intersect<-function(x,y){
    z<-c(max(x[1],y[1]),min(x[2],y[2]))
    if(z[2]<=z[1])
      stop("two interval with no intersect!")
    z
  }
  
  alias=fname
  fname=fname[1]
  
  
  #extending the name
  extend2<-function(gmlSEMfamily.obj,...){
    
    lst<-list(...)
    
    if(length(lst)==0){
      gmlSEMfamily.obj[['extend']]=NULL
      return(gmlSEMfamily.obj)
    }
      
    
    lst.names=names(lst)
    
    if(!is.null(gmlSEMfamily.obj$args.mandatory) &&
       length(gmlSEMfamily.obj$args.mandatory)>0 &&
       !lst.names%in%gmlSEMfamily.obj$args.mandatory){
      args.absent=setdiff(gmlSEMfamily.obj$args.mandatory,
                          lst.names)
      stop("gmlSEM error: The following argument(s) must be pass to the family '",
           gmlSEMfamily.obj$name,"':\n",
           args.absent)
    }
    
    lst.names2=lst.names[!lst.names%in%gmlSEMfamily.obj$args.exposed]
    if(length(lst.names2)>0)
      stop("gmlSEM error: invalid argument(s) '",
           paste0(lst.names2,collapse=", ",sep="'"),"' passed to the family '",
           gmlSEMfamily.obj$name,"'")
    
    #link and link.mean are equivalent
    if("link"%in%lst.names && "link.mean"%in%lst.names)
      stop("gmlSEM error: link and link.mean arguments are equivalent. Both argument are provided for the family: '",
           gmlSEMfamily.obj$name,"'")
    
    if("link"%in%lst.names){
      ind=which(lst.names=="link")
      names(lst)[ind]<-"link.mean"
    }
    
    
    if("link.mean"%in%lst.names){
      if(length(gmlSEMfamily.obj$link.mean)==0)
        stop("gmlSEM error: you can not provide link function for the family '",
             gmlSEMfamily.obj$name,"'")
      
      link.mean<-lst[["link.mean"]]
      if(length(link.mean)>1)
        stop("gmlSEM error: multiple link function provided for family: '",
             gmlSEMfamily.obj$name,"'")
      if(!link.mean%in%gmlSEMfamily.obj$link.mean)
        stop("gmlSEM error: invalid link function provided for family: '",
             gmlSEMfamily.obj$name,"'. Valid link functions are:\n",
             paste0("'",paste0(gmlSEMfamily.obj$link.mean,sep="'",collapse=", "),"'"))
    }
    
    if("link.scale"%in%lst.names){
      if(length(gmlSEMfamily.obj$link.scale)==0)
        stop("gmlSEM error: you can not provide link function for the scale prameter in the family '",
             gmlSEMfamily.obj$name,"'")
      
      link.mean<-lst[["link.scale"]]
      if(length(link.mean)>1)
        stop("gmlSEM error: multiple link.scale function provided for family: '",
             gmlSEMfamily.obj$name,"'")
      if(!link.mean%in%gmlSEMfamily.obj$link.scale)
        stop("gmlSEM error: invalid link.scale function provided for family: '",
             gmlSEMfamily.obj$name,"'. Valid link functions for parameter 'scale' are:\n",
             paste0("'",paste0(gmlSEMfamily.obj$link.scale,sep="'",collapse=", "),"'"))
    }
    
    ###############################
    # Looking for quote arguments and register them in env
    # both 'name' and 'call' arguments can be accepted
    
    # if(any(sapply(lst, function(x)is.call(x)))){
    #   lst1=which(lapply(lst, function(x)is.call(x)))
    #   lst1=paste0(names(lst1),collapse = ",")
    #   stop("gmlSEM error: you can not pass call to the the distribution name '",gmlSEMfamily.obj$name,"' arg(s):",lst1)
    # }
    
    inds=sapply(lst, function(x) is.call(x) || is.name(x))
    
      inds=which(inds)
      runtime.pars=names(lst)[inds]
      fix.parts=setdiff(names(lst),runtime.pars)
      
      runtime.pars=intersect(gmlSEMfamily.obj$args.exposed,runtime.pars)
      fix.parts=intersect(gmlSEMfamily.obj$args.exposed,fix.parts)
      
      gmlSEMfamily.obj$runtime.pars=runtime.pars
        #setdiff(
        #unique(c(gmlSEMfamily.obj$runtime.pars,runtime.pars)),
        #fix.parts)
      
      for(i in  seq_along(runtime.pars)){
        s=runtime.pars[i]
        gmlSEMfamily.obj$env[[s]]=list(quoteValue=lst[[s]],computedValue=NA)
        if(is.name(lst[[s]])){
          gmlSEMfamily.obj$names.env[[as.character(lst[[s]])]]=NA
          attr(gmlSEMfamily.obj$names.env[[as.character(lst[[s]])]],"arg")=
            unique(c(s,
                     attr(gmlSEMfamily.obj$names.env[[as.character(lst[[s]])]],"arg")))
        }
          
        if(is.call(lst[[s]])){
          name.args=ArgSymbol(lst[[s]])
          for(s in name.args)
            gmlSEMfamily.obj$names.env[[as.character(s)]]=NA
          attr(gmlSEMfamily.obj$names.env[[as.character(s)]],"arg")=
            unique(c(s,
                     attr(gmlSEMfamily.obj$names.env[[as.character(s)]],"arg")))
          
        }
          
      }
        
      
        
      
      gmlSEMfamily.obj=updateArg(gmlSEMfamily.obj,...)
    
    
      gmlSEMfamily.obj[['extend']]=NULL #Do not allow deeper extensions
      gmlSEMfamily.obj
  }

  setRunTimeParams<-function(gmlSEMfamily.obj,...){
    lst=list(...)
    fx2=intersect(names(gmlSEMfamily.obj$names.env),names(lst))
    if(length(fx2)==0)
      return(gmlSEMfamily.obj)
    
    lst2=list()
    for(i in seq_along(lst)){
      if(names(lst)[i]%in%fx2)
        lst2[[names(lst)[i]]]=lst[[i]]
    }
    lst=lst2
    lst2=sapply(lst, function(x)is.name(x)||is.call(x)||is.null(x)||is.na(x))
    lst2=which(lst2)
    if(length(lst2)>0)  #No call or symbol is exoected as argument
      stop("gmlSEM error: invalid value for runtime parameter(s): ",paste0(names(lst2),collapse = ","))
    
    lst2=list()
    for(i in seq_along(lst)){
      s=names(lst)[[i]]  # runtime parameter name
      args=attr(gmlSEMfamily.obj$names.env[[s]],"arg")
      gmlSEMfamily.obj$names.env[[s]]=lst[[i]]
      attr(gmlSEMfamily.obj$names.env[[s]],"arg")=args #reset the arg attr
      
      for(arg in args){
        gmlSEMfamily.obj$env[[arg]]$computedValue=
          eval(gmlSEMfamily.obj$env[[arg]]$quoteValue,envir = gmlSEMfamily.obj$names.env)
        lst2[[arg]]=gmlSEMfamily.obj$env[[arg]]$computedValue
      }
    }
      
    lst2$gmlSEMfamily.obj=gmlSEMfamily.obj
    gmlSEMfamily.obj=with(gmlSEMfamily.obj,do.call(gmlSEMfamily.obj$updateArg,lst2))
    
    gmlSEMfamily.obj
  }
  
  updateArg<-function(gmlSEMfamily.obj,...){
    
    #Modify object with new arguments
    lst<-list(...)
    lst.names=names(lst)
    if(length(lst.names)==0)
      return(gmlSEMfamily.obj)
    
    lst.names2=lst.names[!lst.names%in%gmlSEMfamily.obj$args.exposed]
    
    if(length(lst.names2)>0)
      stop("gmlSEM error: invalid argument(s) '",
           paste0(lst.names2,collapse=", "),"' passed to the family: ",
           gmlSEMfamily.obj$name)
    

        
      on.argc.trig.inds=gmlSEMfamily.obj$registeredArg(gmlSEMfamily.obj,lst)
      on.argc.trig=length(on.argc.trig.inds)>0
      fx<-intersect(gmlSEMfamily.obj$args.exposed,names(lst))
      if(length(fx)==0)
        return(gmlSEMfamily.obj)
      
      lst2=lst
      lst=list()
      
      env1=env2=env3=NULL
      
      if(on.argc.trig){
        env1=lapply(on.argc.trig.inds, function(ii)environment(gmlSEMfamily.obj$on.args.change[[ii]]))
      }
        
      if(!is.null(gmlSEMfamily.obj$params) && is.function(gmlSEMfamily.obj$params)){
        env2=environment(gmlSEMfamily.obj$params)      
      }
        
      if(!is.null(gmlSEMfamily.obj$support) && is.function(gmlSEMfamily.obj$support)){
        env3=environment(gmlSEMfamily.obj$support)
      }
        
      
      for(s in fx){
        lst[[s]]=lst2[[s]]
        if(!is.null(env1)){
          for(ii in seq_along(env1))
            assign(s,lst[[s]],env1[[ii]])  
        }
          
        if(!is.null(env2))
          assign(s,lst[[s]],env2)  
        if(!is.null(env3))
          assign(s,lst[[s]],env3)  
      }
      
      lastChangedPars=lst
      
      if(is.null(env1)&&
         !is.function(support)&&
         !is.function(params))
        return(gmlSEMfamily.obj)
      
      if(on.argc.trig){
        for(ii in sort(on.argc.trig.inds)){
          lst1= gmlSEMfamily.obj$call(gmlSEMfamily.obj,on.args.change[[ii]],lastChangedPars,saveScope = TRUE)
          if(is.environment(lst1)){
            if("params"%in%names(lst1))
              gmlSEMfamily.obj$computed$params=lst1[["params"]]
            if("support"%in%names(lst1))
              gmlSEMfamily.obj$computed$support=lst1[["support"]]      
            lst1[['params']]=lst1[['support']]=NULL
            gmlSEMfamily.obj=modifyList(gmlSEMfamily.obj,as.list(lst1))
            
            other.new.params=names(lst1)
            if(length(other.new.params)==0)
              next
            lst1=as.list(lst1)
            lst1=modifyList(lst1,lst1) #remove NULL elements
            lastChangedPars=modifyList(lastChangedPars,lst1)
            
            other.new.params.inds=gmlSEMfamily.obj$registeredArg(gmlSEMfamily.obj,other.new.params)
            other.new.params.inds=setdiff(other.new.params.inds,ii) #avoid recursive runs
            if(length(other.new.params.inds)==0)
              next
            
            #trigger the events once more for new params
            for(ij in sort(other.new.params.inds)){
              lst1=gmlSEMfamily.obj$call(gmlSEMfamily.obj,on.args.change[[ij]],lastChangedPars,saveScope = TRUE)
              if(is.environment(lst1)){
                if("params"%in%names(lst1))
                  gmlSEMfamily.obj$computed$params=lst1[["params"]]
                if("support"%in%names(lst1))
                  gmlSEMfamily.obj$computed$support=lst1[["support"]]      
                lst1[['params']]=lst1[['support']]=NULL
              
                lst1=as.list(lst1)
                lst1=modifyList(lst1,lst1) #remove NULL elements  
                gmlSEMfamily.obj=modifyList(gmlSEMfamily.obj,as.list(lst1))
                lastChangedPars=modifyList(lastChangedPars,lst1)
              }
          }
        }
        }  
      }
        
      gmlSEMfamily.obj$lastChangedPars=lastChangedPars
      
      if(!is.null(env2)){
        lst1 = gmlSEMfamily.obj$call(gmlSEMfamily.obj,params,lastChangedPars)
        if(is.list(lst1)||isS4(lst1))
          gmlSEMfamily.obj$computed$params=lst1
      }  
        
      if(!is.null(env3)){
        lst1=gmlSEMfamily.obj$call(gmlSEMfamily.obj,support,lastChangedPars)
        if(is.list(lst1)||isS4(lst1))
          gmlSEMfamily.obj$computed$support=lst1
      }
        
    gmlSEMfamily.obj
  }
  
  
  #Support is dependent on other arguments
  #Thus it has to be retrieved dynamically 
  support2<-function(gmlSEMfamily.obj){
    
    if(!is.null(gmlSEMfamily.obj$computed$support)){
      gmlSEMfamily.obj$computed$support
    }else if(isS4(gmlSEMfamily.obj$support)&&class(gmlSEMfamily.obj$support)=="support"){
      gmlSEMfamily.obj$support
    }else{
      stop("Bad support specification in family: ",gmlSEMfamily.obj$name)
    }
  }
  
  params2<-function(gmlSEMfamily.obj,par=NA){
    
    if(!is.null(gmlSEMfamily.obj$computed$params)){
      if(is.na(par))
        gmlSEMfamily.obj$computed$params
      else
        gmlSEMfamily.obj$params[[par]]
    }else if(isS4(gmlSEMfamily.obj$params[[1]])&&class(gmlSEMfamily.obj$params[[1]])=="support"){
      if(is.na(par))
        gmlSEMfamily.obj$params
      else
        gmlSEMfamily.obj$params[[par]]
    }else{
      stop("Bad params specification in family: ",gmlSEMfamily.obj$name)
    }
    
}
  
  call<-function(gmlSEMfamily.obj,fun,arglist,saveScope=FALSE){
    if(saveScope){
      nl=length(body(fun))
      body(fun)[[nl+1]]=quote({environment()})
    }
    if(is.null(arglist))
      arglist=list()
    
    
    forms=formals(fun)
    fargs=names(forms)
    
    missing.args=setdiff(fargs,names(arglist))
    
    lst.miss=list()
    
    if(length(missing.args)>0){
      for(m in missing.args){
        if(m=="...")
          next
        if(forms[[m]]!=""){
          stop("gmlSEM error: (developer warning!) Can not use default value for the arguments")
        }else{
          lastVal=gmlSEMfamily.obj[[m]]
          if(is.call(lastVal) || is.name(lastVal) || is.na(lastVal))
            return(NA)  #Parameters has not been resolved yet
          
          lst.miss[[m]]=gmlSEMfamily.obj[[m]]
        }
      }
    }
    
    if(length(lst.miss)>0)
      arglist=modifyList(arglist,lst.miss)
    
    
    nam=names(arglist)
    if(!"..."%in%fargs)
      nam=intersect(fargs,names(arglist))
    
    lst=list()
    for(s in nam)
      lst[[s]]=arglist[[s]]
    
    if(any(sapply(lst, function(m)is.call(m)||is.name(m)||is.na(m))))
      return(NA) #Parameters has not been resolved yet
      
    
    en1=as.environment(gmlSEMfamily.obj)
    parent.env(en1)=parent.frame()
    environment(fun)=en1
    
    if(is.null(fargs[1])){
      do.call(fun,list())
    }else{
      do.call(fun,lst)
    }
    
  }
  
  resolve.symbols<-function(gmlSEMfamily.obj,env=parent.env()){
    if(length(gmlSEMfamily.obj$env)==0)
      return(gmlSEMfamily.obj)
    
    for(i in seq_along(gmlSEMfamily.obj$env)){
      s=names(gmlSEMfamily.obj$env)[i]
      qt=gmlSEMfamily.obj$env[[i]]$quoteValue
        tryCatch({
          newv=eval(s,envir=env)  
          gmlSEMfamily.obj$env[[i]]$computedValue=newv
        },error=function(e){
          #Stick to the last resolved value for the parameter
          if(is.na(gmlSEMfamily.obj$env[[i]]$computedValue))
            stop("gmlSEM error: can not resolve runtime parameter '",
                 s,"' passed to the family: ",
                 gmlSEMfamily.obj$name)
        })
    }
    gmlSEMfamily.obj
  }
 
  
  getSupport<-function(fam){
    
    supps=list(a=NA,b=NA,vals=NA,discrete=FALSE,integer=FALSE,inflated=FALSE)
    attrs=list(a="",b="")
    
    #If distribution is vector-valued then we return a list of supports whose elements are support of each dimension
    if("support"%in%class(fam)){
      dom1=list(fam)
      n.f=1
    }else if(is.list(fam) && "support" %in% class(fam[[1]])){
      dom1=fam
      n.f=length(fam)
    }else{
      dom1=fam$getSupportObject(fam)
      n.f=ifelse("support" %in% class(dom1),1,length(dom1))
      if(!is.list(dom1))
        dom1=list(dom1)
    }
    
    supps.list=list()
    for(n.i in 1:n.f){
      dom=dom1[[n.i]]
      
      if(dom@type=="discrete"){
        supps$vals=supp=dom@support
        supps$discrete=TRUE
      }
      
      
      if(dom@type=="interval"){
        a<-dom@support[1]
        b<-dom@support[2]
        supp=c(a,b)
        supps$a=a
        supps$b=b
        attrs$a=ifelse(dom@include.lhs==1,"[","(")
        attrs$b=ifelse(dom@include.rhs==1,"]",")")
        supps$discrete=FALSE
        
        supp2=supp
        if(!is.logical(fam$truncated[1]))
          supp2=interval.intersect(supp,fam$truncated)
        if(!is.logical(fam$censored[1]))
          supp2=interval.intersect(supp,fam$censored)
        
        if(!is.logical(fam$truncated[1]) || !is.logical(fam$censored[1])){
          supps$a=supp2[1]
          supps$b=supp2[2]
          attrs$a=ifelse(supp[1]==supp2[1],attrs$a,"[")
          attrs$b=ifelse(supp[2]==supp2[2],attrs$b,"]")
          supps$discrete=FALSE
        }
        
        if(dom@is.integer==1){
          supps$discrete=TRUE
          supps$integer=TRUE
          supps$a=ceiling(supps$a)
          supps$a=floor(supps$b)
          attrs$a=ifelse(!is.infinite(supp[1]),"[","(")
          attrs$b=ifelse(!is.infinite(supp[2]),"]",")")
        }
        
        
      }else{
        if(!is.logical(fam$truncated[1]))
          supp=supp[supp>=fam$truncated[1] & supp<=fam$truncated[2]]
        if(!is.logical(fam$censored[1]))
          supp=supp[supp>=fam$censored[1] & supp<=fam$censored[2]]
        
        if(!is.logical(fam$truncated[1]) || !is.logical(fam$censored[1])){
          supps$vals=supp
          supps$discrete=TRUE
        }
        
        if(all(floor(supp)==supp))
          supps$integer=TRUE
      }
      
      if(!is.logical(fam$categorized[1])){
        
        supp=fam$categorized
        supps$vals=supp
        supps$discrete=TRUE
        
        if(all(floor(supp)==supp))
          supps$integer=TRUE
        
      }
      
      
      if(!is.logical(fam$inflated[1])){
        supps$inflated=TRUE
        if(supps$discrete){
          supp=unique(c(supp,fam$inflated))  
          supp=sort(supp)
          supps$vals=supp
        }else{
          supp=sort(supp)
          supps$vals=attr(supp,"inflated")=fam$inflated
        }
      }
      
      supps$atag=attrs$a
      supps$btag=attrs$b
      
      supps.list[[n.i]]=supps
    }
    
    if(n.f==1)
      supps.list=supps.list[[1]]
    
    supps.list
  }
  
  fn<-list(#internal slots
           #Set once, while defining new fmaily
           name          = fname       ,
           alias         = alias       ,
           support       = support     ,
           dim           = dim         ,  #Dimension of the support
           dim.latent    = dim.latent  ,  #sum of the dimension of underlying families' supports (for tobit distribution and psudo families like logit, probit, etc (processes which can be described using an underlying latent continuous process)
           params        = params      ,
           on.args.change=on.args.change,
           
           #set once, internally
           getSupportObject   = support2   ,
           getParamsObject    = params2    ,
           getSupport         = getSupport ,
           
           ########## Important functions and slots ###########
           updateArg          = updateArg        ,  #update dist specific args on runtime
           setRunTimeParams   = setRunTimeParams ,
           env                = new.env(), #for resolving symbols
           names.env          = new.env(), #for storing runtime params wit class 'name'
           ## bin(trials=N)  env$trials=list(quoteValue=quote(N),computedValue=NA)
           ##                names.env$N=NA
           ##                variables in env will be evaluated at names.env
           
           
           #deriving slots, can be changed in extendd dists depends on the syntax
           link.mean     = link.mean   ,
           link.scale    = link.scale  ,
           args.mandatory=args.mandatory,
           truncated     = FALSE       ,
           censored      = FALSE       ,
           categorized   = FALSE       ,
           inflated      = FALSE       ,
           ...                         , 
           #... are specific name parameters such as include_zero in geom() disrtribution, or
           # trials=varName or trials=constant in binomial()
           
           #################################
           #Do not mess with these arguments
           runtime.pars  = runtime.pars,
           extend        = extend2     ,
           args.exposed  = args.exposed,
           invalidated   = FALSE       ,
           call          = call        ,
           nls.onargc    = nls.onargc  ,
           registeredArg = registeredArg
           )
  
  
  if(!is.null(params[[1]])&&isS4(params[[1]]))
    params0=params
  else
    params0=NULL
  
  if(!is.null(support)&&isS4(support))
    support0=support
  else
    support0=NULL
  
  fn$computed=list(params=params0,support=support0)
  
  
  fn$has.finite.discrete.support<-function(fam){
    
    if(length(fam$runtime.pars)>0)
      fam=resolve.symbols(fam)
    
    if(!is.logical(fam$categorized[1]))
      return(TRUE)
    
      dom=fam$getSupportObject(fam)
      if(dom@type=="discrete")
        return(TRUE)
      if(dom@is.integer==1){
        domm=dom@support
        if(!is.logical(fam$truncated[1]))
          domm=interval.intersect(domm,fam$truncated)
        if(!is.logical(fam$censored[1]))
          domm=interval.intersect(domm,fam$censored)
        if(domm[1]>-Inf&domm[2]<Inf)
          return(TRUE)
      } 
      
      return(FALSE)
    
  }
  
  #fam finite discrete support if it is actually finite discrete and not needs data scan
  fn$finite.discrete.support<-function(fam){
    if(length(fam$runtime.pars)>0)
      fam=resolve.symbols(fam)
    
    if(!fam$has.finite.discrete.support(fam))
      stop("gmlSEM error: Family '",fam$name,"' does not have finite discrete support")
    
     supps=fam$getSupport(fam)$vals
  }
  
  fn$is.in.support<-function(fam,val){
    
    supps=fam$getSupport(fam)
    if(val%in%supps$vals)
      return(TRUE)
    
      if(!is.na(supps$a) & (
        ((val>supps$a && val<supps$b)||
         (val==supps$a&&supps$atag=="[")||
         (val==supps$b&&supps$atag=="]")
        ))&&
        (!supps$integer || (supps$integer && val==floor(val)))
      )
        return(TRUE)
    
    
    return(FALSE)
    
  }

  #Set in package
  class(fn) = c(gmlSEMfamily.classname,"list")
  fn
  
}


extendFamily<-function(fname,...){
 res<-sapply(gmlSEMfamilies, function(g)fname%in%g$alias) 
 if(!any(res))
   stop("gmlSEM error: family name '",fname,"' unknown")
 if(sum(res)>1)
   stop("gmlSEM error: family name '",fname,"' is ambiguious")
 
 ind<-which(res)
 fam<-gmlSEMfamilies[[ind]]
 fam<-fam$extend(fam,...)
 fam
}


