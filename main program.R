#######################################
# Load the packages and read the data #
#######################################

rm(list=ls()) # Clear objects in the working space.
par.original <- par(no.readonly=TRUE); # Remember the original plotting parameters.

library(survival)
library(cluster)
library(protoclust)
library(clValid)
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

nrow(data)

# Choose for "delta"
choose_delta<-function(data){
  newdata<-data[data$delta==0 | data$delta==1,]
  return(newdata)}
data<-choose_delta(data)

nrow(data)

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
#w<-max(mydata$time)

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


# THe reason to use if-else syntax is that we want to add legend for combinations only if total number of combination is not too large.
if(numComb<25) #add legend for (groups of combinations) <12 for the K-M Curve
{
  plot(test,ylim=c(.0, 1.0),
       
       xlim=c(0,w),
       xlab=("Survival Time in Months"),
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
       xlab=("Survival Time in Months"),
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



# Construct a vector contain w months' survival probabilities for each combination.
Size<-rep(0,numComb)
SP<-rep(0,numComb)
for (i in 1:(numComb)){
  temp<-testdata[testdata$comb_ind==i,c("time","delta")]
  KM<-KMsurv(temp)
  Size[i]<-nrow(temp)
  pos<-sum(KM[,1]<=w) # Find the position of w months in table of K-M survival probability
  if(pos==0){SP[i]=1}
  else{SP[i]<-KM[pos,2]} # Record w months survival probability to "SP"
}

# Add a label of w-month survival probability to the each leaf of the dendrogram.
label<-paste(100*round(SP,3),":",restcomb)
label2<-paste(100*round(SP,3),":",restcomb,"(n=",Size,")")



# Construct initial dissimilarities matrix "dis0_m". 
dis0_m<-matrix(0, nrow = numComb, ncol=numComb)
for (i in 1:(numComb)){
  for (j in 1:numComb){	
    testi<-as.matrix(testdata[testdata$comb_ind==i,c("time","delta")])
    testj<-as.matrix(testdata[testdata$comb_ind==j,c("time","delta")])
    tempdata<-rbind.data.frame(cbind(testi,1),cbind(testj,2))
    names(tempdata)<-c("time","delta","group")
    tempdata$group<-as.factor(tempdata$group)
    
    dis0_m[i,j]<-abs(ES_WL(tempdata)[2])
  }
}



# Put the dissimilarities in the upper triangular of the matrix "dis0_m" without ones in the diagonal into a vector dis0.
# This is because "learner" function requires an input format of distant matrix.
dis0<-dist(dis0_m)
for (i in 1:(numComb-1))
  for (j in (i+1):numComb)
  {dis0[(i-1)*numComb-i*(i-1)/2+j-i]<-dis0_m[i,j]}


#####################################
# Define functions for Standard PAM #
#####################################


# Define a function "pam.std" to apply Partitioning Around Medoids algorithm. The input of the function is the initial disimialrities
# for "numComb" subjects and a random number K for the number of partitions. "pam.std" will automatically partition subjects into 
# K clusters. The output of the function is a matrix using "TRUE" of "FALSE" to indicate whether any of two subjects are particated
# in the different clusters.
pam.std<-function(dis0,K)
{
  med<-sample(numComb,K,replace = FALSE) # Random assign "K" medoids
  classes<-pam(dis0, K, diss=T, medoids = NULL,do.swap = TRUE)$cluster # vector of group assignment. pam is a system embedded function. return K clusters. 
  dissimilarity<-c()
  for (i in 1:numComb)
  { 
    dissimilarity<-cbind(dissimilarity, classes!=classes[i])
  }
  return(dissimilarity)
}

dis_l<-NULL

ensemble<-function(x,n,Alg){
  
  dissimilarity<-matrix(rep(FALSE, n^2), ncol=n)
  for (i in 1:(n-1)) # aggregate dissimilarity matrix for different number of clusters
  {
    dissimilarity<-dissimilarity + pam.std(x,i)
  }
  
  # Define the dissimilarity matrix for maximum possible number of clusters
  maxdis<-matrix(rep(TRUE, n^2), ncol=n)
  diag(maxdis)<-FALSE
  
  dissimilarity<-dissimilarity  +maxdis
  
  dissimilarity<-dissimilarity / n # find the average of dissimilarity. comments by huan
  dis<-dist(dissimilarity) # create a distance matrix in order to use the structures
  
  k<-1
  for (i in 1:(n-1))
  {
    for (j in (i+1):n)
    {
      dis[k]<-dissimilarity[j,i]
      k      <- k + 1
    }
  }
  
  dis_l<<-dis
  
  if (Alg=="minimax"){
    hc <- protoclust(dis)
    # Run Standard PAM using initial dissimilarities
    protoclust <- list()  # initialize empty object
    # define merging pattern: 
    #    negative numbers are leaves, 
    #    positive are merged clusters (defined by row number in $merge)
    protoclust$merge <- hc$merge
    protoclust$height <- hc$height   # define merge heights
    protoclust$order <- hc$order              # order of leaves(trivial if hand-entered)
    class(protoclust) <- "hclust"   
    
    
    return(protoclust)
    
    
  }
  
  
  
  tempClust<-hclust(dis,method=Alg) # hclust: hiarchical clustering method embeded in R
  return(tempClust)
}


#########################################################################
# Apply Standard PAM to initial dissimilarities and draw the dendrogram #
#########################################################################

# Add a subtitle to show which variables are selected.
if(numVar==1){
  pastename1<-"Analyzed variable is"
  pastename2<-vars
}else{
  pastename1<-"Analyzed variables are"
  pastename<-vars[1]
  i=2
  while(i<numVar){
    pastename<-paste(pastename,vars[i],sep=", ")
    i<-i+1
  }
  pastename2<-paste(pastename,vars[numVar],sep=" and ")
}


# Run Standard PAM using initial dissimilarities
myClust<-ensemble(dis0,numComb,linkage)


plot(myClust,labels=as.character(label), xlab=NA,main="",ylab="Dissimilarity",sub=NA, cex=min(40/numComb,0.6))


################################################
# Use C-index to give suggest number of groups #
################################################



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

if (finalGroups==0){
  
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
  
  
  temporder<--myClust$order # Starting with leaves' order in the dendrogram, construct a vector "temporder" temporarily recording groups' order in the dendrogram.  
  # Note that a negative number means it is an index of combination that has not been merged to any group yet, a positive number means it is an index of new formed group.
  tempassign<-as.list(myClust$order) # Starting with leaves' order in the dendrogram, construct a list "tempassign" temporarily recording combinations' order in the dendrogram
  
  o_assign<-function(assign){
    
    # Construct an assigndata to order the assignment
    assigndata<-testdata 
    
    maxt<-max(assigndata$time)
    # Construct a vector assignE_T to store the expected survival time for each group
    assignE_T<-1:length(assign) 
    for(j in 1:length(assign)){
      assigndatai<-assigndata[assigndata$comb_ind %in% unlist(assign[j]),c("time","delta")]
      assignE_T[j]<-Esurv_c(assigndatai,maxt)
    }
    
    assignorder<-order(-assignE_T)
    assign[assignorder]
  }
  
  allC<-NULL # Construct a vector "allC" to record the C-index value of each grouping
  allDI_0<-NULL # Construct a vector "allDI" to record the Dunn Index of each grouping based on initial dissimilarities
  allDI_l<-NULL # Construct a vector "allDI" to record the Dunn Index of each grouping based on learned dissimilarities
  
  for (i in 1:(numComb)){
 
    
    ########################################
    # Save current ordered assignment "ordered_tempassign" #
    ########################################
    ordered_tempassign<-o_assign(tempassign)
    for(s in 1:length(ordered_tempassign)){
      ordered_tempassign[[s]]<-restcomb[ordered_tempassign[[s]]]
    }
    
    assign_i<-paste("prog","_assign",numComb-i+1,".Rdata",sep="")
    save(tempassign, file=assign_i)
    ordered_assign_i<-paste("prog","_ordered_assign",numComb-i+1,".Rdata",sep="")
    save(ordered_tempassign, file=ordered_assign_i)
    
    
    ###########################################################
    # Calculate C-index using current assignment "tempassign" #
    ###########################################################
    
    # Construct Cindexdata which is used to input in function "C_index"
    Cindexdata<-testdata 
    
    # Add a column "group_ind" to record the group number for each case
    Cindexdata$group_ind<-0 
    for(j in 1:length(tempassign)){
      Cindexdata[Cindexdata$comb_ind %in% unlist(tempassign[j]),]$group_ind<-j
    }
    
    maxt<-max(Cindexdata$time)
    # Add a column "E_T" to record the estimated survival time for each case
    Cindexdata$E_T<-0
    for(k in 1:length(tempassign)){
      group_k_data<-Cindexdata[Cindexdata$group_ind==k,c("time","delta")]
      Cindexdata$E_T[Cindexdata$group_ind==k]<-Esurv_c(group_k_data,maxt)
    }
    
    Cindexdata<-Cindexdata[,c("time","delta","E_T")] # Keep only 3 columns of interest
    
    allC <- cbind(allC, C_index(Cindexdata)) # Record the C-index in vector "allC"
    
    ##############################################################
    # Calculate Dunn index using current assignment "tempassign" #
    ##############################################################
    
    tempcluster<-rep(seq_along(tempassign), lapply(tempassign, length))[order(unlist(tempassign))]

    allDI_0 <- cbind(allDI_0, dunn(dis0, tempcluster))
    allDI_l <- cbind(allDI_l, dunn(dis_l, tempcluster)) 
    

    
    #######################
    # Update "tempassign" #
    #######################
    if(i==numComb) break
    
    mergedpair<-myClust$merge[i,] # find merged pair of groups in the ith merging step
    
    p<-which(mergedpair[1]==temporder) # find the position of first element of "mergedpair" in the "temporder"
    q<-which(mergedpair[2]==temporder) # find the position of second element of "mergedpair" in the "temporder"
    
    # Merge pth row and qth row in "tempassign" as a new row. (Note that according to the property of "myClust$order", pth row and qth row should be next to each other, i.e. |p-q|=1)
    tempassign[[p]]<-sort(c(unlist(tempassign[[p]]),unlist(tempassign[[q]]))) # put the merged combinations in pth row
    tempassign[[q]]<-NULL # Remove qth row
    
    # Update "temporder" (replace pth and qth element in "temporder" with a new group index)
    temporder[p]<-i # The index of new formed group in ith merging step in "temporder" is i.
    temporder<-temporder[-q] # Remove qth element in "temporder" 
    
  }
  
  allC<-rev(allC) # Reverse the order of elements in allC
  
  save(allC, file=paste("prog","_allC.Rdata",sep=""))
  
  increase<-rep(0,length(allC)-1)
  for(i in 1:(length(allC)-1)){
    increase[i]<-diff(allC)[i]/allC[i]
  }
  round(increase,5)

}



bestn<-8 

# Plot C-indices against number of groups and give n*
plot(0,0,xlim = c(1,numComb),ylim = c(0.5,min(1,max(allC)+0.05)),cex.axis=1.5,type = "n",main="",xlab="Number of Groups",ylab="C-index",cex.lab=1.5)
lines(1:numComb, allC  , col =" gray ", lty =1, lwd =2)
abline(v=bestn,col='red')
mtext(paste("n* =",bestn), cex = 1.5,adj=0.1, padj=0.2,side=1,at=bestn)
points(bestn, allC[bestn],type="p", pch=19, col="red")
text(bestn, y = allC[bestn], labels =paste("C-index = ",round(allC[bestn],4),sep=''), adj = c(-0.05,+1.2),
     pos = NULL, offset = 0.5, vfont = NULL,
     cex = 1.5, col = NULL, font = NULL)


finalGroups<-bestn


plot(myClust,labels=as.character(label), xlab=NA,main="",ylab="Dissimilarity",sub=NA, cex=min(40/numComb,0.6))
assign<-rect.hclust(myClust, k=finalGroups, border="red") # draw boxes



########################################
# Plot K-M curves for the final groups #
########################################

groupdata<-testdata
legend<-rep(0,length(assign))
ES<-rep(0,length(assign)) # Calculate estimated survival time for each group
groupdata$group_ind<-0
maxt<-max(groupdata$time)

for(i in 1:length(assign)){
  groupdata[groupdata$comb_ind %in% unlist(assign[i]),]$group_ind<-i
  ES[i]<-Esurv_c(groupdata[groupdata$group_ind==i,c("time","delta")],maxt)
}

ESorder<-order(-ES)
assign<-assign[ESorder]

groupdata$r_group_ind<-0 # Reorder "group_ind" based on estimated survival time
for(i in 1:length(assign)){
  groupdata[groupdata$comb_ind %in% unlist(assign[i]),]$r_group_ind<-i
  legend[i]<-paste("Group",i," n=",sum(groupdata$r_group_ind==i),sep="")
}


grouptest<-survfit(Surv(time,delta)~group_ind, data=groupdata)


plot(grouptest,ylim=c(0, 1.0),
     
     xlim=c(0,max(data$time)),
     xlab=("Survival Time by Months"),
     ylab=("Proportion Surviving"),
     cex=1,
     cex.lab=1.5,
     cex.axis=1.5,
     col=2:(finalGroups+1),
     lwd=3.0,
     lty=ceiling(((2:(finalGroups+1))-1)/8),
     mark.time=FALSE)
title(main = "")
legend("bottomleft",
       legend=legend,
       col=ESorder+1,
       cex=0.7,
       lty=ceiling((ESorder)/8),
       lwd=3.0,
       xjust=1, yjust=1)



plot(myClust,labels=as.character(label), xlab=NA,main="Original Dendrogram",ylab="Dissimilarity",sub='', cex=min(30/numComb,1.5))
assign<-rect.hclust(myClust, k=finalGroups, border="red") # draw boxes

# Define a function "Inv_Per" for inverse permutation
Inv_Per<-function(Per){
  inv<-rep(0,length(Per))
  for (i in 1:length(Per)){
    inv[i]<-which(ESorder==i)
  }
  return(inv)
}

# Find the inverse permutation of "ESorder"
Inv_ESorder<-Inv_Per(ESorder)


# Add group number beneath the dendrogram
x<-0
for (i in 1:length(assign)){
  x_i<-x+0.5*length(assign[[i]])+0.5
  mtext(Inv_ESorder[i],side=1,at=x_i)
  x<-x+length(assign[[i]])
}




# Rotate the groups in the dendrogram in ascending order of severity from left to right.

k<-0
for (i in 1:length(assign)){
  for (j in 1:length(assign[[i]])){
    assign[[i]][j]<-myClust$order[k+j] 
  }
  k<-k+length(assign[[i]])
}
assign<-assign[ESorder] # Reorder the rows of "assign" by E_T
myClust$order<-unlist(assign)



plot(myClust,labels=as.character(label),cex.axis=1.5, xlab=NA,main="",ylab="Dissimilarity",sub='',cex.lab=1.5, cex=min(80/numComb,1.5))
assign<-rect.hclust(myClust, k=finalGroups, border="red") # draw boxes


# Add group number beneath the dendrogram
x<-0
for (i in 1:length(assign)){
  x_i<-x+0.5*length(assign[[i]])+0.5
  mtext(i,side=1,at=x_i,cex=1.5,padj = 0.3)
  x<-x+length(assign[[i]])
}



grouptest<-survfit(Surv(time,delta)~r_group_ind, data=groupdata)

plot(grouptest,ylim=c(0, 1.0),
     
     xlim=c(0,144),
     xlab=("Survival Time by Months"),
     ylab=("Proportion Surviving"),
     cex=1,
     cex.lab=1.5,
     cex.axis=1.5,
     col=c("gray50","cyan3","orange" , 6,  3, 2, 4, 9),
     lwd=3.0,
     lty=1,
     mark.time=FALSE)
title(main = "")
legend("bottomleft",
       legend=legend,
       col=c("gray50","cyan3","orange" , 6,  3, 2, 4, 9),
       cex=1,
       lty=1,
       lwd=3.0,
       xjust=1, yjust=1)
