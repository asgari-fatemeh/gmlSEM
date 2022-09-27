#In package


gmlSEMfamilies<-list()

family.alias<-list(
  binomial      = c("bin","binomial")     ,
  discrete   = c("categorical","cat","disc","discrete")  ,
  bernoulli     = c("bernoulli","bern")   ,
  multinomial   = c("multinomial")        ,
  gaussian      = c("gaussian","normal")  ,
  
  # halfnormal(link="log")
  halfnormal    = c("halfnormal")         ,
  skewnormal    = c("skewnormal")         ,
  
  #does it exist in stan::math?
  invnormal     = c("invgaussian","inversegaussian",
  "inverse.gaussian","inv.normal","inverse.normal","inversenormal"),
  
  chisquared    = c("chisquared","chi2","chisq","chisqr"),
  t             = c("t")                  ,
  logistic      = c("logistic")           ,
  gamma         = c("gamma")              , #link="log"
  invgamma      = c("invgamma","inversegamma","inverse.gamma")           , #link="log"
  exponential   = c("exponential","exp")  , #link="log"
  gumble        = c("gumble")             ,
  beta          = c("beta")               ,
  lognormal     = c("lognormal","logn")   ,
  weibull       = c("weibull","weib")     ,
  poisson       = c("poisson")            ,
  negbin        = c("negbin","negativebinomial"),
  geom          = c("geometric","geo","geom") , #include_zero=NA  #When NA it is decided by scanning data or TRUE otherwise
  hypergeom     = c("hypergeometric","hypergeom","hypergeo"),
  
  #Shorthands
  zib           = c("zib","bezi")         ,
  zoib          = c("zoib")               , #inflated=c(0,1)
  zib           = c("zib")                ,
  zip           = c("zip")                ,
  lceb          = c("lceb")               ,  #family=poisson
  hurdle        = c("hurdle")             ,
  tobit         = c("tobit")    #family=normal,support=NA) //A shorthand for gaussian(censored=support)
)


known.links.mean  <-c("logit","probit","cauchit","sqr","log","cloglog","identity","inverse")
known.links.scale <-c("logit","probit","cauchit","sqr","log","cloglog","identity","inverse")
# TODO check link functions
family.links<-list(
  binomial      = list(link.mean=c("logit","probit","cauchit"))     ,
  discrete   = list(link.mean=known.links.mean)  ,
  bernoulli     = list(link.mean=known.links.mean)  ,
  multinomial   = list(link.mean=known.links.mean)  ,
  gaussian      = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  
  # halfnormal(link="log")
  halfnormal    = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  skewnormal    = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  chisquared    = list(link.mean=known.links.mean)  ,
  t             = list(link.mean=known.links.mean)  ,
  logistic      = list(link.mean=known.links.mean)  ,
  gamma         = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  invgamma      = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  exponential   = list(link.mean=known.links.mean)  ,
  gumble        = list(link.mean=known.links.mean)  ,
  beta          = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,               
  lognormal     = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  weibull       = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  poisson       = list(link.mean=known.links.mean)  ,
  negbin        = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  geom          = list(link.mean=known.links.mean)  ,
  hypergeom     = list(link.mean=known.links.mean)  ,
  
  #Shorthands
  zib           = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,               
  zoib          = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,               
  zib           = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,               
  zip           = list(link.mean=known.links.mean)  ,
  lceb          = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,               
  hurdle        = list(link.mean=known.links.mean)  ,
  tobit         = list() #No link function
)


gmlSEMfamilies[["geom"]]<-newFamily(family.alias[["geom"]],
              params=c(mu=newSupport("interval",c(-Inf,Inf))),
              link.mean=family.links[["geom"]]$link.mean,
              include.zero=FALSE,
              support=function(){
                if(include.zero){
                  newSupport("interval",c(0,Inf),is.integer = TRUE)
                }else{
                  newSupport("interval",c(1,Inf),is.integer = TRUE)
                }
              })

fn2=extendFamily("geom",include.zero=TRUE)
fn3=extendFamily("geom",include.zero=FALSE)

fn2$getSupport(fn2)
fn2$getSupportObject(fn2)
fn2$getParamsObject(fn2)

fn2$getSupport(fn2)
fn3$getSupport(fn3)

gmlSEMfamilies[["copula"]]<-newFamily(c("copula","cop"),
                                    support=newSupport("interval",c(-Inf,Inf)),
                                    type="gaussian",
                                    on.args.change=function(type){
                                      #Nothing
                                      type=tolower(type)
                                      type=match.arg(type,c("gaussian","archimedean"))
                                    })

gmlSEMfamilies[["beta"]]<-newFamily(family.alias[["beta"]],
                                    support=newSupport("interval",c(0,1)),
                                    params=c(mu=newSupport("interval",c(0,1))),
                                    link.mean=family.links[["beta"]]$link.mean,
                                    link.scale=family.links[["beta"]]$link.scale,
                                    domain=c(0,1),
                                    on.args.change=function(domain){
                                      support=newSupport("interval",domain)
                                    })

f11<-extendFamily("beta",domain=c(-1,1),censored=c(0,1),inflated=0)
f11<-extendFamily("beta",domain=c(-1,1),inflated=0)
f11$getSupport(f11)
f11$is.in.support(f11,1)


gmlSEMfamilies[["NULL"]]<-newFamily("NULL",dim=0)

gmlSEMfamilies[["bernoulli"]]<-newFamily(family.alias[["bernoulli"]],
                                    support=newSupport("discrete",c(0,1)),
                                    params=c(mu=newSupport("interval",c(0,1))),
                                    link.mean=family.links[["bernoulli"]]$link.mean)

gmlSEMfamilies[["binomial"]]<-newFamily(family.alias[["binomial"]],
                                         params=c(mu=newSupport("interval",c(0,1))),
                                         link.mean=family.links[["binomial"]]$link.mean,
                                        on.args.change = function(trials){
                                          support=newSupport("discrete",0:trials)
                                        },
                                        trials=1 #Default value
                                        )

f1<-extendFamily("bin",trials=5)
f1$getSupport(f1)


f2=extendFamily("bin",trials=quote(N))
f2=f2$setRunTimeParams(f2,N=3)
f2$getSupport(f2)


gmlSEMfamilies[["gaussian"]]<-newFamily(family.alias[["gaussian"]],
                                        support=newSupport("interval",c(-Inf,Inf)),
                                        params=c(mu=newSupport("interval",c(-Inf,Inf)),
                                                 sigma2=newSupport("interval",c(0,Inf),include.lhs = TRUE)),
                                        link.mean=family.links[["gaussian"]]$link.mean,
                                        link.scale=family.links[["gaussian"]]$link.scale)


f1<-extendFamily("normal",link="identity",link.scale="log")
f1$getSupport(f1)

gmlSEMfamilies[["discrete"]]<-newFamily(family.alias[["discrete"]],
                                        link.mean=family.links[["discrete"]]$link.mean,
                                        N=1, #Default value
                                        domain=0:1,
                                        p=c(0.5,0.5),
                                        on.args.change = list(
                                          function(domain){
                                            N=length(domain)
                                            p=rep(1/N,N)
                                            support=newSupport("discrete",domain)
                                          },
                                          function(N){
                                            domain=1:N
                                            p=rep(1/N,N)
                                            support=newSupport("discrete",domain)
                                          },
                                          function(p){
                                            if(sum(p)!=1)
                                              stop("gmlSEM error: sum of probabilities must equals to one in 'discrete' family.")
                                            N=length(p)
                                            domain=1:N
                                            support=newSupport("discrete",domain)
                                          }
                                        ))

f1=extendFamily("discrete",domain=quote(dom))
f1=f1$setRunTimeParams(f1,dom=c(1,2,3,5))
f1$getSupport(f1)


gmlSEMfamilies[["ordered"]]<-newFamily("ordered",
                                     params = c(mu=family.links[["discrete"]]$link.mean),
                                     family=extendFamily("normal"),
                                     domain=c(0,1),
                                     thresholds=c("t1"),
                                     on.args.change = list(function(domain){
                                       if(length(domain)<2)
                                         stop("gmlSEM error: domain of family 'ordered' must be of length grater than one.")
                                       
                                       if(length(thresholds)!=length(domain))
                                         thresholds=paste0("t",1:(length(domain)-1))
                                       support=newSupport("discrete",domain)
                                     },
                                     function(thresholds){
                                       if(length(thresholds)!=(length(domain)-1))
                                         stop("gmlSEM error: 'threshold' of family 'ordered' must be of equal to the length of 'domain' minus one.")
                                     }))

f2=extendFamily("ordered",domain=c(1,2,3))
f2$getSupport(f2)



gmlSEMfamilies[["tobit"]]<-newFamily("tobit",
                                      params = c(mu=family.links[["discrete"]]$link.mean),
                                      dim=1,
                                      type=1,
                                      family=extendFamily("normal"),
                                      family1=extendFamily("normal"),
                                      family2=extendFamily("normal"),
                                      family3=extendFamily("normal"),
                                       on.args.change = list(function(type){
                                         if(type==1 || type==2){
                                           dim=1
                                           dim.latent=type
                                         }else if(type==3){
                                           dim=dim.latent=2
                                         }else if(type==4){
                                           dim=3
                                           dim.latent=3
                                         }else if(type==5){
                                           dim=2
                                           dim.latent=3
                                         }else{
                                           stop("gmlSEM error: type must be a value of 1,2,..,5 in Tobit model")
                                         }
                                       },
                                       function(domain){
                                         
                                       })
                                     )


gmlSEMfamilies2<-list(
  binomial      = list(link.mean=c("logit","probit","cauchit"))     ,
  categorical   = list(link.mean=known.links.mean)  ,
  bernoulli     = list(link.mean=known.links.mean)  ,
  multinomial   = list(link.mean=known.links.mean)  ,
  gaussian      = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  
  # halfnormal(link="log")
  halfnormal    = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  skewnormal    = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  chisquared    = list(link.mean=known.links.mean)  ,
  t             = list(link.mean=known.links.mean)  ,
  logistic      = list(link.mean=known.links.mean)  ,
  gamma         = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  invgamma      = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  exponential   = list(link.mean=known.links.mean)  ,
  gumble        = list(link.mean=known.links.mean)  ,
  beta          = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,               
  lognormal     = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  weibull       = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  poisson       = list(link.mean=known.links.mean)  ,
  negbin        = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,
  geometric     = list(link.mean=known.links.mean)  ,
  hypergeometric= list(link.mean=known.links.mean)  ,
  
  #Shorthands
  zib           = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,               
  zoib          = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,               
  zib           = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,               
  zip           = list(link.mean=known.links.mean)  ,
  lceb          = list(link.mean=known.links.mean,link.scale=known.links.scale)  ,               
  hurdle        = list(link.mean=known.links.mean)  ,
  tobit         = list() #No link function
)

getFamilyID<-function(fname){
  ind<-sapply(gmlSEMfamilies,function(g)fname%in%g$alias)
  if(any(ind)){
    names(which(ind))
  }else{
    NA
  }
}


