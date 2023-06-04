##########################
# Weighted logrank tests #
##########################

# Construct a function "WL" to compute the (weighted) logrank test statistics

WL<-function(data){
  data[,3]<-as.numeric(factor(data[,3]))
  
  #Count the frequences for each combination of time and delta
  if (all(data[,2]==1))
  {freq<-cbind(0,table(data[,1],data[,2]))
  }  else if (all(data[,2]==0))
  {freq<-cbind(table(data[,1],data[,2]),1)
  }  else{freq<-table(data[,1],data[,2])}
  
  Tj<-sort(unique(data[,1]))
  Dj<-as.numeric(freq[,2]) # Count how many deaths before t
  Cj<-as.numeric(freq[,1]) # Count how many censored cases before t
  
  nrow=length(Tj);
  Yj<-rep(0,nrow);
  
  data1<-data[data[,3]==1,]
  Dj1<-rep(0,nrow);
  Cj1<-rep(0,nrow);
  Yj1<-rep(0,nrow);
  
  n1<-nrow(data[data[,3]==1,])
  n2<-nrow(data[data[,3]==2,])
  
  Sj<-rep(0,nrow+1)
  Sj[1]<-1
  
  LNj1<-rep(0,nrow)
  GNj1<-rep(0,nrow)
  PNj1<-rep(0,nrow)
  
  LDj1<-rep(0,nrow)
  GDj1<-rep(0,nrow)
  PDj1<-rep(0,nrow)
  
  for(s in 1:nrow){ 
    Yj[s]<-sum(data[,1]>=Tj[s])
    Dj1[s]<-sum(data1[data1[,2]==1,1]==Tj[s]) # Count how many deaths in group 1 at time s
    Cj1[s]<-sum(data1[data1[,2]==0,1]==Tj[s]) # Count how many censored cases in group 1 at time s
    Yj1[s]<-sum(data1[,1]>=Tj[s])
    
    Sj[s+1]<-Sj[s]*(1-Dj[s]/(Yj[s]))# In Sj, [s+1] refers to time=s
    
    
    LNj1[s]<-(Dj1[s]-Yj1[s]*(Dj[s]/Yj[s]))
    GNj1[s]<-(Yj[s]/(n1+n2))*(Dj1[s]-Yj1[s]*(Dj[s]/Yj[s]))
    PNj1[s]<-Sj[s]*(Dj1[s]-Yj1[s]*(Dj[s]/Yj[s]))
    
    
    LDj1[s]<-Dj[s]*(Yj1[s]/Yj[s])*(1-Yj1[s]/Yj[s])*(Yj[s]-Dj[s])/max(Yj[s]-1,1)
    GDj1[s]<-(Yj[s]/(n1+n2))^2*Dj[s]*(Yj1[s]/Yj[s])*(1-Yj1[s]/Yj[s])*(Yj[s]-Dj[s])/max(Yj[s]-1,1)
    PDj1[s]<-Sj[s]^2*Dj[s]*(Yj1[s]/Yj[s])*(1-Yj1[s]/Yj[s])*(Yj[s]-Dj[s])/max(Yj[s]-1,1)
    
  }
  
  L<-sum(LNj1)/sqrt(sum(LDj1))
  G<-sum(GNj1)/sqrt(sum(GDj1))
  P<-sum(PNj1)/sqrt(sum(PDj1))
  
  PL<-pchisq(L^2,1,lower.tail = F)
  PG<-pchisq(G^2,1,lower.tail = F)
  PP<-pchisq(P^2,1,lower.tail = F)
  
  
  return(c(L,PL,G,PG,P,PP))
}

################################################
# Effect sizes based on weighted logrank tests #
################################################

# Construct a function "ES_WL" to compute the (weighted) logrank test statistics.

ES_WL<-function(data){
  data[,3]<-as.numeric(factor(data[,3]))
  
  #Count the frequences for each combination of time and delta
  if (all(data[,2]==1))
  {freq<-cbind(0,table(data[,1],data[,2]))
  }  else if (all(data[,2]==0))
  {freq<-cbind(table(data[,1],data[,2]),1)
  }  else{freq<-table(data[,1],data[,2])}
  
  Tj<-sort(unique(data[,1]))
  Dj<-as.numeric(freq[,2]) # Count how many deaths before t
  Cj<-as.numeric(freq[,1]) # Count how many censored cases before t
  
  nrow=length(Tj);
  Yj<-rep(0,nrow);
  
  data1<-data[data[,3]==1,]
  Dj1<-rep(0,nrow);
  Cj1<-rep(0,nrow);
  Yj1<-rep(0,nrow);
  
  n1<-nrow(data[data[,3]==1,])
  n2<-nrow(data[data[,3]==2,])
  
  Sj<-rep(0,nrow+1)
  Sj[1]<-1
  
  ESLj1<-rep(0,nrow)
  ESGj1<-rep(0,nrow)
  ESPj1<-rep(0,nrow)
  
  for(s in 1:nrow){ 
    Yj[s]<-sum(data[,1]>=Tj[s])
    Dj1[s]<-sum(data1[data1[,2]==1,1]==Tj[s]) # Count how many deaths in group 1 at time s
    Cj1[s]<-sum(data1[data1[,2]==0,1]==Tj[s]) # Count how many censored cases in group 1 at time s
    Yj1[s]<-sum(data1[,1]>=Tj[s])
    
    Sj[s+1]<-Sj[s]*(1-Dj[s]/(Yj[s]))# In Sj, [s+1] refers to time=s
    
    
    ESLj1[s]<-(Dj1[s]-Yj1[s]*(Dj[s]/Yj[s]))
    ESGj1[s]<-(Yj[s]/(n1+n2))*(Dj1[s]-Yj1[s]*(Dj[s]/Yj[s]))
    ESPj1[s]<-Sj[s]*(Dj1[s]-Yj1[s]*(Dj[s]/Yj[s]))
    
  }
  
  ESL<-sum(ESLj1)*(n1+n2)/(n1*n2)
  ESG<-sum(ESGj1)*(n1+n2)/(n1*n2)
  ESP<-sum(ESPj1)*(n1+n2)/(n1*n2)
  
  
  return(c(ESL,ESG,ESP))
}

##########
# ES_MWE #
##########

# Use exponential tails to get Dhat_e 
Dhat_e<-function(Data1,Data2,tau1,tau2){
  Data<-as.data.frame(rbind(cbind(Data1,1),cbind(Data2,2)))
  names(Data)<-c("time","delta","ind")
  
  if (all(Data$delta==0)) {return(0.5)}
  times<-Data$time[Data$delta==1]
  
  min<-min(tau1,tau2)
  max<-max(tau1,tau2)
  
  test<-survfit(Surv(time,delta)~ind, data=Data) 
  
  #Part A
  time1<-unique(times[times<=min])
  time1<-time1[order(time1)]
  
  mysurv1<-summary(test,times=time1)
  
  S1<-mysurv1$surv[mysurv1$strata=="ind=1"]
  S2<-mysurv1$surv[mysurv1$strata=="ind=2"]
  
  ap1<-length(time1)-length(S1)
  ap2<-length(time1)-length(S2)
  
  S1<-c(1,S1,rep(0,ap1))
  S2<-c(1,S2,rep(0,ap2))
  
  Int1<-0
  for (i in 2:length(S1)){
    Int1<-sum(Int1,S1[i]*(S2[i-1]-S2[i]))
  }
  
  #Part B
  
  # Find maximum failure time for group 1
  max1<-max(test$time[1:test$strata[1]])
  if(max1<tau1){
    Stau1<-min(test$surv[1:test$strata[1]])
  } else if(tau1<min(test$time[1:test$strata[1]])){
    Stau1<-1
  } else{
    Stau1<-test$surv[1:test$strata[1]][max(which(test$time[1:test$strata[1]]<=tau1))]
  }
  
  lambda_1<--log(Stau1)/tau1
  
  # Find maximum failure time for group 2
  max2<-max(test$time[-c(1:test$strata[1])])
  if(max2<tau2){
    Stau2<-min(test$surv[-c(1:test$strata[1])])
  }  else if(tau2<min(test$time[-c(1:test$strata[1])])){
    Stau2<-1
  }  else{
    Stau2<-test$surv[-c(1:test$strata[1])][max(which(test$time[-c(1:test$strata[1])]<=tau2))]
  }
  
  lambda_2<--log(Stau2)/tau2
  
  if (tau1==tau2){
    Int2<-0
  }
  
  if (tau1<tau2){
    if(sum(tau1<times & tau2>=times)==0){
      Int2<-0
    }
    else {
      time2<-unique(c(tau1,times[tau1<times & tau2>=times],tau2))
      time2<-time2[order(time2)]
      mysurv2<-summary(test,times=time2)
      S1<-exp(-lambda_1*time2)
      S2<-mysurv2$surv[mysurv2$strata=="ind=2"]
      Int2<-0
      for (i in 2:length(S2)){
        Int2<-sum(Int2,S1[i]*(S2[i-1]-S2[i]))
      }
    }
    
  }
  
  if (tau1>tau2){
    time2<-unique(c(tau1,times[tau2<times & tau1>=times],tau2))
    time2<-time2[order(time2)]
    mysurv2<-summary(test,times=time2)
    S1<-mysurv2$surv[mysurv2$strata=="ind=1"]
    S2<-exp(-lambda_2*time2)
    Int2<-0
    for (i in 2:length(S2)){
      Int2<-sum(Int2,S1[i-1]*(exp(-lambda_2*time2[i-1])-exp(-lambda_2*time2[i])))
    }
  }
  
  #Part C
  
  if (is.infinite(lambda_1)){ Int3<-0
  }  else if (is.infinite(lambda_2)){ Int3<-0
  }  else {Int3<-(lambda_2/(lambda_1+lambda_2))*exp(-(lambda_1+lambda_2)*max)}
  
  Int<-sum(Int1,Int2,Int3)
  
  return(Int)
}

minimax<-function(data){
  maxt_each<-tapply(data$time,data$group, max)
  return(min(maxt_each))
}

# Construct EShat_MWE. 
ES_MWE<-function(data,tau1=minimax(data),tau2=minimax(data)){
  Data1<-as.matrix(data[data$group==1,c("time","delta")])
  Data2<-as.matrix(data[data$group==2,c("time","delta")])
  return((Dhat_e(Data2,Data1,tau2,tau1)-Dhat_e(Data1,Data2,tau1,tau2)))
}

##########
# ES_MWC #
##########

Dhat_c<-function(x,y,tau){
  testData<-as.data.frame(rbind(cbind(x,1),cbind(y,2)))
  names(testData)<-c("time","delta","ind")
  #  numGroup<-2
  if (all(testData$delta==0)) {return(0.5)}
  times<-testData$time[testData$delta==1]
  tmax1<-max(x[,1])
  tmax2<-max(y[,1])

  test<-survfit(Surv(time,delta)~ind, data=testData) 
  
  #Part A
  time1<-unique(times[times<=tau])
  time1<-time1[order(time1)]
  
  mysurv1<-summary(test,times=time1)
  
  S1<-mysurv1$surv[mysurv1$strata=="ind=1"]
  S2<-mysurv1$surv[mysurv1$strata=="ind=2"]
  
  ####### append 0s 
  ap1<-length(time1)-length(S1)
  ap2<-length(time1)-length(S2)
  
  S1<-c(1,S1,rep(0,ap1))
  S2<-c(1,S2,rep(0,ap2))
  
  Int1<-0
  for (i in 2:length(S1)){
    Int1<-sum(Int1,S1[i]*(S2[i-1]-S2[i]))
  }
  
  S1min<-S1[length(S1)]
  S2min<-S2[length(S2)]
  
  Int<-Int1/(1-S1min*S2min)
  
  return(Int)
}

# Construct EShat_MWC
ES_MWC<-function(data,tau=minimax(data)){
  Data1<-as.matrix(data[data$group==1,c("time","delta")])
  Data2<-as.matrix(data[data$group==2,c("time","delta")])
  return((Dhat_c(Data2,Data1,tau)-Dhat_c(Data1,Data2,tau)))
}