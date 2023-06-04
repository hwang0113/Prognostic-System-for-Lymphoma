#######################################
# Load the packages and read the data #
#######################################

rm(list=ls()) # Clear objects in the working space.
par.original <- par(no.readonly=TRUE); # Remember the orginal plotting parameters.

library(survival)
library(cluster)
library(protoclust)
library("readxl")


source("functions.R")

rawdata<-read_excel("lymphomas.xlsx",sheet = 1)

names(rawdata)<-c("id","time","delta","age","sex","race","site","site2","H","H_C","H_BG","year","3T","3N","3M","6T","6N","6M","7T","7N","7M","CSE","6S","7S","SSF1","SSF2","SSF3","SSF4","SSF5","AAS","BSym")

rawdata<-rawdata[rawdata$year>=2014 & rawdata$year<=2014,]

data<-rawdata



data$S<-""

data[data$"6S"=="IA","S"]    <-"IA     "
data[data$"6S"=="IB","S"]    <-"IB     "
data[data$"6S"=="IEA","S"]   <-"IEA    "
data[data$"6S"=="IEB","S"]   <-"IEB    "
data[data$"6S"=="ISA","S"]   <-"ISA    "
data[data$"6S"=="ISB","S"]   <-"ISB    "
data[data$"6S"=="IIA","S"]   <-"IIA    "
data[data$"6S"=="IIB","S"]   <-"IIB    "
data[data$"6S"=="IIEA","S"]  <-"IIEA   "
data[data$"6S"=="IIEB","S"]  <-"IIEB   "
data[data$"6S"=="IISA","S"]  <-"IISA   "
data[data$"6S"=="IISB","S"]  <-"IISB   "
data[data$"6S"=="IIESA","S"] <-"IIESA  "
data[data$"6S"=="IIESB","S"] <-"IIESB  "
data[data$"6S"=="IIIA","S"]  <-"IIIA   "
data[data$"6S"=="IIIB","S"]  <-"IIIB   "
data[data$"6S"=="IIIEA","S"] <-"IIIEA  "
data[data$"6S"=="IIIEB","S"] <-"IIIEB  "
data[data$"6S"=="IIISA","S"] <-"IIISA  "
data[data$"6S"=="IIISB","S"] <-"IIISB  "
data[data$"6S"=="IIIESA","S"]<-"IIIESA "
data[data$"6S"=="IIIESB","S"]<-"IIIESB "
data[data$"6S"=="IVA","S"]   <-"IVA    "
data[data$"6S"=="IVB","S"]   <-"IVB    "

table(data$S)


data$I<-""

data[data$site=="Hodgkin - Extranodal" ,"I"]<-"EN "
data[data$site=="Hodgkin - Nodal" ,"I"]<-"N  "
data[data$site=="NHL - Extranodal" ,"I"]<-"EN "
data[data$site=="NHL - Nodal" ,"I"]<-"N  "
table(data$I)

data$T<-""

data[data$site=="Hodgkin - Extranodal" ,"T"]<-"HL  "
data[data$site=="Hodgkin - Nodal" ,"T"]<-"HL  "
data[data$site=="NHL - Extranodal" ,"T"]<-"NHL "
data[data$site=="NHL - Nodal" ,"T"]<-"NHL "
table(data$T)


data$A<-""
data[data$age>=0 & data$age<13 ,"A"]<-"Y "
data[data$age>=13 & data$age<20 ,"A"]<-"O "
table(data$A)

data$X<-""
data[data$sex==1,"X"]<-"M"
data[data$sex==2,"X"]<-"F"
table(data$X)



nrow(data)


# Choose for "time"
choose_time<-function(data){
  newdata<-data[data$time!=9999,]
  return(newdata)}
data<-choose_time(data)

# Choose for "delta"
choose_delta<-function(data){
  newdata<-data[data$delta==0 | data$delta==1,]
  return(newdata)}
data<-choose_delta(data)



# Choose for "S"
choose_S<-function(data){
  newdata<- data[data$S%in%c("IA     ",
                             "IB     ",
                             "IEA    ",
                             "IEB    ",
                             "ISA    ",
                             "ISB    ",
                             "IIA    ",
                             "IIB    ",
                             "IIEA   ",
                             "IIEB   ",
                             "IISA   ",
                             "IISB   ",
                             "IIESA  ",
                             "IIESB  ",
                             "IIIA   ",
                             "IIIB   ",
                             "IIIEA  ",
                             "IIIEB  ",
                             "IIISA  ",
                             "IIISB  ",
                             "IIIESA ",
                             "IIIESB ",
                             "IVA    ",
                             "IVB    ")  ,]
  return(newdata)}
data<-choose_S(data)
table(data$S)

data$S<-factor(data$S,levels=c("IA     ",
                               "IB     ",
                               "IEA    ",
                               "IEB    ",
                               "ISA    ",
                               "ISB    ",
                               "IIA    ",
                               "IIB    ",
                               "IIEA   ",
                               "IIEB   ",
                               "IISA   ",
                               "IISB   ",
                               "IIESA  ",
                               "IIESB  ",
                               "IIIA   ",
                               "IIIB   ",
                               "IIIEA  ",
                               "IIIEB  ",
                               "IIISA  ",
                               "IIISB  ",
                               "IIIESA ",
                               "IIIESB ",
                               "IVA    ",
                               "IVB    "))

table(data$S)

# Choose for "I"
choose_I<-function(data){
  newdata<-data[data$I%in%c("EN ","N  ") ,]
  return(newdata)}
data<-choose_I(data)

table(data$I)


# Choose for "T"
choose_T<-function(data){
  newdata<-data[data$T%in%c("HL  ","NHL ") ,]
  return(newdata)}
data<-choose_T(data)

table(data$T)


# Choose for "A"
choose_A<-function(data){
  newdata<-data[data$A%in%c("Y ","O ") ,]
  return(newdata)}
data<-choose_A(data)

table(data$A)

# Choose for "X"
choose_X<-function(data){
  newdata<-data[data$X%in%c("M","F") ,]
  return(newdata)}
data<-choose_X(data)


table(data$S)
table(data$I)
table(data$T)
table(data$A)
table(data$X)



round(100*table(data$S)/nrow(data),1)
round(100*table(data$I)/nrow(data),1)
round(100*table(data$T)/nrow(data),1)
round(100*table(data$A)/nrow(data),1)
round(100*table(data$X)/nrow(data),1)





########################################################

vars<-c("S","I","T","A","X")   # Construct a vector "vars" that contains the name of the factors we are interested in. The names must be exact the same as the factor names in the inputted data.


elimSize<-0  # set a number "elimSize" so that any combination with size < "elimSize" will be eliminated.

displayFinalKM<-1 # 1--display survival curves of "boxes"; 0--don't display

varIndex=match(vars,names(data)) #Find the index of positions of variables
numVar<-length(varIndex) #Count the number of variables



mydata=cbind.data.frame(data[,varIndex],time=data$time,delta=data$delta)  # Keep only "time", "delta" and interested factors.

# Add a column "comb" to "mydata" to indicate the combinations (the reason of using if-else syntax is that "paste" function has to be used with at least two variables.)
if(numVar==1){
  comb<-data[vars]
  mydata<-cbind.data.frame(mydata,comb)
  colnames(mydata)[1]<-"comb"
}else{
  mydata$comb<-apply(data[,varIndex],1,paste,collapse="") 
}

mydata<-cbind(A=data$A,H=data$H,mydata)

# "Count" is a variable that counts number of observation for each combination.
Count<-aggregate(count ~ comb, cbind(count=1,mydata),length)

# Names of Combinations with size > elimSize will be recorded in the variable "restcomb"
restcomb<-Count[Count$count>=elimSize,]$comb

# Count how many different combinations left.
numComb<-length(restcomb)

# Construct a "testdata" that only keeps combinations in "restcomb"
data=mydata[mydata$comb %in% restcomb,]



######################
# Set the parameters #
######################

linkage<-"complete" # Choose a linkage method
vars<-c("S","I","T","A","X")    # Construct a vector "vars" that contains the name of the subset factors we are interested in. The names must be exact the same as the factor names in the inputted data.


elimSize<-0  # set a number "elimSize" so that any combination with size < "elimSize" will be eliminated.
finalGroups<-0 # Set the number of clusters in the final results (1 < "finalGroups" < total number of combinations (with size smaller than "elimSize") ) If 0 is input, then the number of groups will be decided by C-index
displayFinalKM<-1 # 1--display survival curves of "boxes"; 0--don't display
w<-60 # Set the study time by months

varIndex=match(vars,names(data)) #Find the index of positions of variables
numVar<-length(varIndex) #Count the number of variables


mydata=cbind.data.frame(data[,varIndex],time=data$time,delta=data$delta)  # Keep only "time", "delta" and interested factors.

# Add a column "comb" to "mydata" to indicate the combinations (the reason of using if-else syntax is that "paste" function has to be used with at least two variables.)
if(numVar==1){
  comb<-data[vars]
  mydata<-cbind.data.frame(mydata,comb)
  colnames(mydata)[1]<-"comb"
}else{
  mydata$comb<-apply(data[,varIndex],1,paste,collapse="") 
}


mydata<-cbind(A=data$A,H=data$H,mydata)

# "Count" is a variable that counts number of observation for each combination.
Count<-aggregate(count ~ comb, cbind(count=1,mydata),length)

# Names of Combinations with size > elimSize will be recorded in the variable "restcomb"
restcomb<-Count[Count$count>=elimSize,]$comb

# Count how many different combinations left.
numComb<-length(restcomb)

# Construct a "testdata" that only keeps combinations in "restcomb"
testdata=mydata[mydata$comb %in% restcomb,]

# Add a column "comb_int" to represent the index of "comb" in the "restcomb" for each observation (in alphabeta order).
testdata$comb_ind<-as.integer(factor(testdata$comb))

nrow(testdata)

##################################################################################




################################################
# Plot K-M survival curve for each combination #
################################################

test<-survfit(Surv(time,delta)~comb, data=testdata) # Use function "survfit" in package "survival" to plot the K-M survival curves


# The reason to use if-else syntax is that we want to add legend for combinations only if total number of combination is not too large.
if(numComb<25) #add legend for (groups of combinations) <12 for the K-M Curve
{
  plot(test,ylim=c(.0, 1.0),
       
       xlim=c(0,w),
       xlab=("Survival Time in Month"),
       ylab=("Proportion Surviving"),
       cex=1,
       cex.lab=1.5,
       cex.axis=1.5,
       col=1:6,
       lwd=3.0,
       lty=c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4),
       mark.time=FALSE)
  title(main = "Kaplan-Meier Curves for \nDifferent Combinations")
  legend("bottomleft",
         legend=restcomb,
         col=1:6,
         cex=0.7,
         lty=c(1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4),
         lwd=3.0,
         xjust=1, yjust=1)
}else
{
  plot(test,ylim=c(.0, 1.0),
       
       xlim=c(0,w),
       xlab=("Survival Time in Month"),
       ylab=("Proportion Surviving"),
       cex=1,
       cex.lab=1.5,
       cex.axis=1.5,
       col=1:6,
       lwd=3.0,
       lty=c(1,1,1,1,1,1,4,4,4,4,4,4),
       mark.time=FALSE)
  title(main = "Kaplan-Meier Curves for \nDifferent Combinations")
}


##########################################################
# Find w-month survival probability for each combination #
##########################################################


# Define a Functions "KMsurv" to calulate K-M survival probability, The input of this function is a data set with two columns "time" and "delta". 
# The output of the function is a data set with 2 columns "Ti" (distinct death time) and "Si" (corresponding K-M survival probability)

KMsurv<-function(data){
  #Count the frequences for each combination of time and delta
  if (all(data$delta==1))
  {freq<-cbind(0,table(data[,1],data[,2]))}
  else if (all(data$delta==0))
  {return(rbind.data.frame(c(1,1),c(max(data[,1]),1)))}
  else{freq<-table(data[,1],data[,2])}
  
  Ti<-as.numeric(rownames(freq))
  Di<-as.numeric(freq[,2]) # Count how many deaths before ti
  Ci<-as.numeric(freq[,1]) # Count how many censored cases before ti
  out<-cumsum(c(0,Ci+Di)) # Count how many cases dropped before ti
  nrow=length(Ti);
  Yi<-rep(nrow(data),nrow);
  Qi<-rep(0,nrow);
  Pi<-rep(0,nrow);
  
  for(s in 1:nrow){
    Yi[s]<-Yi[s]-out[s]
    Qi[s]<-1-Di[s]/Yi[s]
    
  }
  
  Si<-cumprod(Qi[1:nrow]) # Calculate K-M survival probability
  
  return(matrix(cbind(Ti,Si)[Di>=1,],ncol=2))
}



################################################
# Use C-index to give suggest number of groups #
################################################

# Define a function "Esurv" to calculated estimated survival time.
Esurv<-function(data){
  surv<-KMsurv(data)
  Es<-surv[1,1]*1
  if(nrow(surv)>1){
    for(i in 1:(nrow(surv)-1)){
      Es<-sum(Es,surv[i,2]*(surv[i+1,1]-surv[i,1]))
    }
  }
  # Add the part after the last death
  max_t<-data[which.max(data[,1]),1]
  Es<-sum(Es,surv[nrow(surv),2]*(max_t-surv[nrow(surv),1]))
  return(Es)
}

# Define a function "Esurv_c" to calculated estimated survival time.
Esurv_c<-function(data,maxt){
  surv<-KMsurv(data)
  if(maxt<surv[nrow(surv),1]){stop("overall_maxt < maxt")}
  Es<-surv[1,1]*1
  if(nrow(surv)>1){
    for(i in 1:(nrow(surv)-1)){
      Es<-sum(Es,surv[i,2]*(surv[i+1,1]-surv[i,1]))
    }
  }
  # Add the part after the last death
  Es<-sum(Es,surv[nrow(surv),2]*(maxt-surv[nrow(surv),1]))
  return(Es)
}


  
  ###############################################################
  # Define function "C_index" that will be used to find C-index #
  ###############################################################
  
  # Define a function C_index with an arguments "data" to calculate Harrell's C-index. "data" has 3 columns "time", "delta" and "E_T".
  C_index<-function(data){
    
    # Aggregate the rows in data by adding a column "count"
    data <- aggregate(count ~ ., cbind(count=1,data),length) 
    
    # Sort the data by order of time (ascending) and delta (descending)
    ord <- order(data[,1], -data[,2]) 
    data <- as.matrix(data[ord,])
    
    # Find the row positions of deaths. Minus 1 as the index of element in C starts at 0 instead of 1
    di <- which(data[,2] == 1)-1
    
    library(Rcpp) # Call "Rcpp" package to use C function
    cppFunction('NumericVector Cindex(NumericMatrix data,NumericVector di){ 
                int nrow = data.nrow();
                int ndi = di.size();
                double numerator = 0;
                double denominator = 0;
                NumericVector out(2);
                
                int ii = 0;
                int i = di[ii];
                
                NumericVector f(3);
                NumericVector g(3);
                
                while(ii<ndi){
                if(i==(nrow-1)){break;}
                
                for (int j=i+1;j<nrow;j++){
                int count = data(i,3)*data(j,3);
                
                if (data(j,0) > data(i,0)) {
                denominator += count;
                if (data(j,2) > data(i,2))
                {numerator += count; }
                else if (data(j,2) == data(i,2)) 
                {numerator += 0.5*count; }
                
                }else{
                if (data(j,1)==0){
                denominator += count;
                if (data(j,2) > data(i,2)) {
                numerator += count;}
                else if (data(j,2) == data(i,2)) 
                {numerator += 0.5*count;}
                }
                }
                }
                i = di[++ii];
                }
                out[0]=numerator;
                out[1]=denominator;
                return out;
  }
                ')
    result<-Cindex(data,di)
    return(result[[1]]/result[[2]])
  }






################
## Validation ##
################



# Construct the Cindex_validation function which is used to find C-index for the validation dataset
Cindex_validation<-function(Cindexdata_validation){
  # Add a column "group_ind" to match the groups using previous prognostic system
  Cindexdata_validation$E_T<-0 
  for(j in 1:length(ordered_tempassign)){
    Cindexdata_validation[Cindexdata_validation$comb %in% unlist(ordered_tempassign[j]),"E_T"]<-length(ordered_tempassign)+1-j
  }
  
  # Delete unmatched combinations
  Cindexdata_validation<-Cindexdata_validation[Cindexdata_validation$E_T!=0,]
  
  
  Cindexdata_validation<-Cindexdata_validation[,c("time","delta","E_T")] # Keep only 3 columns of interest

  
  validation_allC<-C_index(Cindexdata_validation)
}

load(paste("prog","_allC.Rdata",sep=""))

allC_validation<-NULL

Cindexdata_validation<-testdata

load(paste("prog","_ordered_assign",1,".Rdata",sep=""))
Cindexdata_validation<-Cindexdata_validation[Cindexdata_validation$comb%in%ordered_tempassign[[1]],]


for (i in 1:length(allC)){
  
  load(paste("prog","_ordered_assign",i,".Rdata",sep=""))
  
  tempC<-Cindex_validation(Cindexdata_validation)
  
  
  
  allC_validation<-c(allC_validation,tempC)
  
}

print("cindex_validation")
print(allC_validation)


numComb<-length(table(Cindexdata_validation$comb))

nstar<-8    


# Plot C-indices and new C-indices against number of groups 

plot(0,0,xlim = c(1,length(allC)),ylim = c(0.5,min(1,max(allC)+0.05)),type = "n",main="C-indices against Number of Groups",xlab="Number of Groups",ylab="C-index")
lines(1:length(allC), allC  , col =" gray ", lty =1, lwd =2)
lines(1:length(allC), allC_validation  , col =" blue ", lty =1, lwd =2)
abline(v=nstar,col='red')

mtext(paste("n* =",nstar),side=1,at=nstar)




# Plot C-indices against number of groups and give n*
plot(0,0,xlim = c(1,length(allC)),ylim = c(0.5,min(1,max(allC)+0.05)),cex.axis=0.4,type = "n",main="",xlab="Number of Groups",ylab="C-index",cex.lab=1)
lines(1:length(allC), allC  , col =" red ", lty =1, lwd =1)
abline(v=nstar,col='gray', lwd =0.5)
mtext(paste("n* =",nstar), cex = 0.6,adj=0.5, padj=0.2,side=1,at=nstar, col = "red")
points(nstar, allC[nstar],type="p", pch=19, col="red",cex = 0.4)
text(nstar, y = allC[nstar], labels =paste("C-index = ",formatC(allC[nstar], format='f',digit=4),sep=''), adj = c(-0.05,+1.2),
     pos = NULL, offset = 0.5, vfont = NULL,
     cex = 0.6, col = "red", font = NULL)

lines(1:length(allC_validation), allC_validation  , col =" darkblue ", lty =1, lwd =1)
points(nstar, allC_validation[nstar],type="p", pch=19, col="darkblue",cex = 0.4)
text(nstar, y = allC_validation[nstar], labels =paste("C-index = ",formatC(allC_validation[nstar], format='f',digit=4),sep=''), adj = c(-0.05,-1.2),
     pos = NULL, offset = 0.5, vfont = NULL,
     cex = 0.6, col = "darkblue", font = NULL)

legend("bottomright",
       legend=c("Training set","Validation set"),
       col=c("red","darkblue"),
       cex=0.6,
       lty=1,
       lwd=1,bg="white" ,
       xjust=1, yjust=1)


plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="")




########################################
# Plot K-M curves for the final groups #
########################################

groupdata<-Cindexdata_validation

load(paste("prog","_ordered_assign",nstar,".Rdata",sep=""))

legend<-rep(0,length(ordered_tempassign))
ES<-rep(0,length(ordered_tempassign)) # Calculate estimated survival time for each group
ESorder<-order(-ES)
groupdata$group_ind<-0

groupdata$r_group_ind<-0 # Reorder "group_ind" based on estimated survival time
for(i in 1:length(ordered_tempassign)){
  groupdata[groupdata$comb %in% unlist(ordered_tempassign[i]),]$r_group_ind<-i
  legend[i]<-paste("Group",i," n=",sum(groupdata$r_group_ind==i),sep="")
}


grouptest<-survfit(Surv(time,delta)~r_group_ind, data=groupdata)


plot(grouptest,ylim=c(0, 1.0),
     
     xlim=c(0,max(data$time)),
     xlab=("Survival Time by Months"),
     ylab=("Proportion Surviving"),
     cex=1,
     cex.lab=1.5,
     cex.axis=1.5,
     col=2:(nstar+1),
     lwd=3.0,
     lty=ceiling(((2:(nstar+1))-1)/8),
     mark.time=FALSE)
title(main = "")
legend("bottomleft",
       legend=legend,
       col=ESorder+1,
       cex=0.7,
       lty=ceiling((ESorder)/8),
       lwd=3.0,
       xjust=1, yjust=1)



library("scales")
library("ggthemes")
library("ggsci")

par(mar=c(5.1, 4.1, 4.1, 3.1))
summ<-summary(grouptest,time=60)
plot(grouptest,ylim=c(0, 1.0),
     
     xlim=c(0,60),
     xlab=("Survival Time in Months"),
     ylab=("Proportion Surviving"),
     cex=1.2,
     cex.lab=1,
     cex.axis=1,
     col=gdocs_pal()(10)[2:9],
     lwd=2.0,
     lty=c(2,2,2,2,1,1,1,1),
     mark.time=FALSE)
title(main = "")
legend("bottomleft",
       legend=legend,
       col=gdocs_pal()(10)[2:9],
       cex=0.7,
       lty=c(2,2,2,2,1,1,1,1),
       lwd=2.0,
       xjust=1, yjust=1)


rysurv<-paste(round(100*summ$surv, 1), "%", sep="")
text(66,cex = 0.6,summ$surv,rysurv,col=gdocs_pal()(10)[2:9], xpd=NA)

par(mar=c(5.1, 4.1, 4.1, 2.1))


# Construct initial dissimilarities matrix "dis0_m". 
Comp_Surv<-matrix(0, nrow = nstar, ncol=nstar)
for (i in 1:(nstar)){
  for (j in 1:nstar){	
    testi<-as.matrix(groupdata[groupdata$r_group_ind==i,c("time","delta")])
    testj<-as.matrix(groupdata[groupdata$r_group_ind==j,c("time","delta")])
    tempdata<-rbind.data.frame(cbind(testi,1),cbind(testj,2))
    names(tempdata)<-c("time","delta","group")
    tempdata$group<-as.factor(tempdata$group)
    
    Comp_Surv[i,j]<-WL(tempdata)[2]
  }
}

Comp_Surv

# Construct initial dissimilarities matrix "dis0_m". 
Comp_Surv<-matrix(0, nrow = nstar, ncol=nstar)
for (i in 1:(nstar-1)){
  j=i+1
  
  testi<-as.matrix(groupdata[groupdata$r_group_ind==i,c("time","delta")])
  testj<-as.matrix(groupdata[groupdata$r_group_ind==j,c("time","delta")])
  tempdata<-rbind.data.frame(cbind(testi,1),cbind(testj,2))
  names(tempdata)<-c("time","delta","group")
  tempdata$group<-as.factor(tempdata$group)
  
      print(WL(tempdata)[2])
  fit_coxph<-coxph(Surv(time, delta) ~ group, data = tempdata)
  print(fit_coxph)
  print(summary(fit_coxph)$coefficients)
  print(pchisq((summary(fit_coxph)$coefficients[4])^2, df=1, lower.tail = FALSE))
  print(exp(confint(fit_coxph)))
}
