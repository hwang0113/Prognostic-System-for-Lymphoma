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

rawdata<-rawdata[rawdata$year>=2004 & rawdata$year<=2013,]

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


elimSize<-25  # set a number "elimSize" so that any combination with size < "elimSize" will be eliminated.

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

linkage<-"minimax" # Choose a linkage method
vars<-c("S","I","T","A","X")    # Construct a vector "vars" that contains the name of the subset factors we are interested in. The names must be exact the same as the factor names in the inputted data.


elimSize<-25  # set a number "elimSize" so that any combination with size < "elimSize" will be eliminated.
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


##################################################################################


################################################
# Plot K-M survival curve for each combination #
################################################
legend<-rep(0,24)


legend[1] <-paste("Stage ","IA      "," n=",sum(testdata$S=="IA     "),sep="")
legend[2] <-paste("Stage ","IB      "," n=",sum(testdata$S=="IB     "),sep="")
legend[3] <-paste("Stage ","IEA    "," n=",sum(testdata$S=="IEA    "),sep="")
legend[4] <-paste("Stage ","IEB    "," n=",sum(testdata$S=="IEB    "),sep="")
legend[5] <-paste("Stage ","ISA    "," n=",sum(testdata$S=="ISA    "),sep="")
legend[6] <-paste("Stage ","ISB    "," n=",sum(testdata$S=="ISB    "),sep="")
legend[7] <-paste("Stage ","IIA     "," n=",sum(testdata$S=="IIA    "),sep="")
legend[8] <-paste("Stage ","IIB     "," n=",sum(testdata$S=="IIB    "),sep="")
legend[9] <-paste("Stage ","IIEA   "," n=",sum(testdata$S=="IIEA   "),sep="")
legend[10]<-paste("Stage ","IIEB    "," n=",sum(testdata$S=="IIEB   "),sep="")
legend[11]<-paste("Stage ","IISA    "," n=",sum(testdata$S=="IISA   "),sep="")
legend[12]<-paste("Stage ","IISB    "," n=",sum(testdata$S=="IISB   "),sep="")
legend[13]<-paste("Stage ","IIESA  "," n=",sum(testdata$S=="IIESA  "),sep="")
legend[14]<-paste("Stage ","IIESB  "," n=",sum(testdata$S=="IIESB  "),sep="")
legend[15]<-paste("Stage ","IIIA     "," n=",sum(testdata$S=="IIIA   "),sep="")
legend[16]<-paste("Stage ","IIIB     "," n=",sum(testdata$S=="IIIB   "),sep="")
legend[17]<-paste("Stage ","IIIEA   "," n=",sum(testdata$S=="IIIEA  "),sep="")
legend[18]<-paste("Stage ","IIIEB   "," n=",sum(testdata$S=="IIIEB  "),sep="")
legend[19]<-paste("Stage ","IIISA   "," n=",sum(testdata$S=="IIISA  "),sep="")
legend[20]<-paste("Stage ","IIISB   "," n=",sum(testdata$S=="IIISB  "),sep="")
legend[21]<-paste("Stage ","IIIESA "," n=",sum(testdata$S=="IIIESA "),sep="")
legend[22]<-paste("Stage ","IIIESB "," n=",sum(testdata$S=="IIIESB "),sep="")
legend[23]<-paste("Stage ","IVA     "," n=",sum(testdata$S=="IVA    "),sep="")
legend[24]<-paste("Stage ","IVB    "," n=",sum(testdata$S=="IVB    "),sep="")

test<-survfit(Surv(time,delta)~S, data=testdata) # Use function "survfit" in package "survival" to plot the K-M survival curves

color<-c(2,"orange" ,"cyan2", "#33a02c" ,"purple","#4575b4",
         2,"orange" ,"cyan2", "#33a02c" ,"purple","#4575b4","rosybrown4",9,
         2,"orange" ,"cyan2", "#33a02c" ,"purple","#4575b4","rosybrown4",9,
         2,"orange")

summ<-summary(test,time=60)

par(mar=c(5.1, 4.1, 4.1, 8.5),xpd=F) 
plot(test,ylim=c(0.5, 1.0),
     xlim=c(0,60),
     xlab=("Survival Time in Months"),
     ylab=("Proportion Surviving"),
     cex=2,
     cex.lab=1,
     cex.axis=1,
     col=color,
     lwd=2,
     lty=c(rep(1,6),rep(2,8),rep(3,8),4,4),
     mark.time=FALSE)
title(main = "")

legend("right", 
       legend=legend,  
       col=color,
       cex=0.7,
       lty=c(rep(1,6),rep(2,8),rep(3,8),4,4),
       lwd=2.0,
       xjust=1, yjust=1)





###########################################################################################


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



###########################################################################################
Cindexdata<-testdata 

table(Cindexdata$S)


maxt<-max(Cindexdata$time)
# Add a column "E_T" to record the estimated survival time for each case
Cindexdata$E_T<-0


Cindexdata$E_T[Cindexdata$S=="IA     "]<-Esurv_c(Cindexdata[Cindexdata$S=="IA     ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IB     "]<-Esurv_c(Cindexdata[Cindexdata$S=="IB     ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IEA    "]<-Esurv_c(Cindexdata[Cindexdata$S=="IEA    ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IEB    "]<-Esurv_c(Cindexdata[Cindexdata$S=="IEB    ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="ISA    "]<-Esurv_c(Cindexdata[Cindexdata$S=="ISA    ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="ISB    "]<-Esurv_c(Cindexdata[Cindexdata$S=="ISB    ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIA    "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIA    ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIB    "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIB    ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIEA   "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIEA   ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIEB   "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIEB   ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IISA   "]<-Esurv_c(Cindexdata[Cindexdata$S=="IISA   ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IISB   "]<-Esurv_c(Cindexdata[Cindexdata$S=="IISB   ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIESA  "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIESA  ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIESB  "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIESB  ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIIA   "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIIA   ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIIB   "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIIB   ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIIEA  "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIIEA  ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIIEB  "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIIEB  ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIISA  "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIISA  ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIISB  "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIISB  ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIIESA "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIIESA ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IIIESB "]<-Esurv_c(Cindexdata[Cindexdata$S=="IIIESB ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IVA    "]<-Esurv_c(Cindexdata[Cindexdata$S=="IVA    ",c("time","delta")],maxt)
Cindexdata$E_T[Cindexdata$S=="IVB    "]<-Esurv_c(Cindexdata[Cindexdata$S=="IVB    ",c("time","delta")],maxt)


C_index(Cindexdata[,c("time","delta","E_T")] )

##########################################################################################################


###########################################################################################
Cindexdata<-testdata 

table(Cindexdata$S)


maxt<-max(Cindexdata$time)
# Add a column "E_T" to record the estimated survival time for each case
Cindexdata$E_T<-0


Cindexdata$E_T[Cindexdata$S=="IA     "]<-24
Cindexdata$E_T[Cindexdata$S=="IB     "]<-23
Cindexdata$E_T[Cindexdata$S=="IEA    "]<-22
Cindexdata$E_T[Cindexdata$S=="IEB    "]<-21
Cindexdata$E_T[Cindexdata$S=="ISA    "]<-20
Cindexdata$E_T[Cindexdata$S=="ISB    "]<-19
Cindexdata$E_T[Cindexdata$S=="IIA    "]<-18
Cindexdata$E_T[Cindexdata$S=="IIB    "]<-17
Cindexdata$E_T[Cindexdata$S=="IIEA   "]<-16
Cindexdata$E_T[Cindexdata$S=="IIEB   "]<-15
Cindexdata$E_T[Cindexdata$S=="IISA   "]<-14
Cindexdata$E_T[Cindexdata$S=="IISB   "]<-13
Cindexdata$E_T[Cindexdata$S=="IIESA  "]<-12
Cindexdata$E_T[Cindexdata$S=="IIESB  "]<-11
Cindexdata$E_T[Cindexdata$S=="IIIA   "]<-10
Cindexdata$E_T[Cindexdata$S=="IIIB   "]<-9
Cindexdata$E_T[Cindexdata$S=="IIIEA  "]<-8
Cindexdata$E_T[Cindexdata$S=="IIIEB  "]<-7
Cindexdata$E_T[Cindexdata$S=="IIISA  "]<-6
Cindexdata$E_T[Cindexdata$S=="IIISB  "]<-5
Cindexdata$E_T[Cindexdata$S=="IIIESA "]<-4
Cindexdata$E_T[Cindexdata$S=="IIIESB "]<-3
Cindexdata$E_T[Cindexdata$S=="IVA    "]<-2
Cindexdata$E_T[Cindexdata$S=="IVB    "]<-1

C_index(Cindexdata[,c("time","delta","E_T")] )
