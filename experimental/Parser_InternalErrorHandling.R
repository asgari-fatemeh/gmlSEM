#Error and warnings handling


gmlSEMwarnings = list()

gmlSEMerror<-function(...){
  msg=paste0(lapply(list(...), function(m)as.character(m)),collapse = "")
  stop("\ngmlSEM error:\n",msg)
}

gmlSEMWarning<-function(msg){
  gmlSEMwarnings[[length(gmlSEMwarnings)+1]]<<-msg
}
