# NHANES Dietary Data
# Eating habits and food consumption data
# Longitudinal Survey from 2000 to 2022
# Each year includes two days food consumption data for more than 100K individuals

# Current data: One day consumption of 120K individuals
# Indicators are skewed and zero-inflated

data<-read.csv("D:\\R\\ZISEM\\Data\\NHANES Individual Food Consumption Day.csv")
labs<-read.csv("D:\\R\\ZISEM\\Data\\NHANES var names.csv")

inds=labs$Variable.Name%in%colnames(data)
labs<-labs[inds,1:2]
labs<-unique(labs)


#Caffeine
plott<-function(vname){
  lab=(labs$Variable.Description[labs$Variable.Name==vname])[1]
  x<-data[,vname]
  y<-x[x!=0]
  zeros=mean(x==0,na.rm = T)*100
  hist(y,probability = T,main=paste0(lab,"+zero(",round(zeros,1),"%)"))  
}


par(mfrow=c(2,2))
par(mar=c(4,2,2,2))

#print all variables
print(labs)
plott(vname="DR1ICAFF")
plott(vname="DR1ICALC")
plott(vname="DR1IVC")
plott(vname="DR1IVD")


data2=data[,colnames(data)%in%labs$Variable.Name]
cnames=colnames(data2)
colnames(data2)<-paste0("y",1:48)
labs2=data.frame(org.vname=cnames,new.name=colnames(data2),desc=NA)
for(i in 1:nrow(labs2)){
  labs2$desc[i]=(labs$Variable.Description[labs$Variable.Name==labs2$org.vname[i]])[1]
}

print(labs2)


model<-c('
         family: y1,...,y48 halfnormal(inflated=0)
         
         vt as Vitamin
         mn as minerals
         
         vt =~ y14+y22+y23+y25+y31+...+y35
         mn =~ y36+...+y42
         
         vt ~ age * gender
         mn ~ age * gender
         ')

res=gmlSEM(model,data)

