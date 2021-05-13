library(MASS)
library(tidyverse)
library(dplyr)

dlong <- function(data){
  # prepare the data set for the long-format
  dat1 <- data %>%
    dplyr::select(studlab, t1, n1, r1) %>%
    rename("study" = "studlab", 
           "treat" = "t1", 
           "events" = "r1", 
           "n" = "n1")
  
  dat2 <- data %>%
    dplyr::select(studlab, t2, n2, r2) %>%
    rename("study" = "studlab", 
           "treat" = "t2", 
           "events" = "r2", 
           "n" = "n2")
  
  data <- rbind(dat1, dat2)
  data <- data %>%
    mutate(ID = paste(study, treat, sep ="")) %>%
    filter(!duplicated(ID)) %>%
    dplyr::select(-ID)
  
  data=data[order(data$study),]
  return(data)
}

dlong1 <- function(data){
  # prepare the data set for the long-format
  dat1 <- data %>%
    dplyr::select(studlab, t1, n1, r1) %>%
    rename("study" = "studlab", 
           "treatment" = "t1", 
           "responders" = "r1", 
           "sampleSize" = "n1")
  
  dat2 <- data %>%
    dplyr::select(studlab, t2, n2, r2) %>%
    rename("study" = "studlab", 
           "treatment" = "t2", 
           "responders" = "r2", 
           "sampleSize" = "n2")
  
  data <- rbind(dat1, dat2)
  data <- data %>%
    mutate(ID = paste(study, treatment, sep ="")) %>%
    filter(!duplicated(ID)) %>%
    dplyr::select(-ID)
  
  data=data[order(data$study),]
  return(data)
}




set.seed(44)
N.sim=1000


data1=list()
data2=list()
data3=list()
logOR=list()
OR=list()
NT=5 #### number of treatments in the network
NS=2 #### number of studies per comparison
tau=0 ####  heterogeneity SD
Npmin=30 #### minimum number of patients per arm
Npmax=60 #### maximum number of patients per arm

### define treatment indices
t1=c()
t2=c()
for (i in 1:(NT-1)){
  for (k in (i+1):NT){
    for(j in 1:NS){
      t1=c(t1,i)
      t2=c(t2,k)      }}}
N.stud=length(t1)

### define patients per treatment arm
for (i in 1:N.sim)
{   
  logOR[[i]]=seq(from =1/(NT-1), to = 1, by = 1/(NT-1))  ### true logOR across studies 
  OR[[i]]=c(1,exp(logOR[[i]]))
  
  
  data1[[i]]=data.frame(t1,t2)
  data1[[i]]$studlab=c(1:(N.stud))
  data1[[i]]$n1=data1[[i]]$n2=round(runif(N.stud,Npmin,Npmax))}


#### define probabilities per treatment, per study arm
for (i in 1:N.sim)
{  
  data1[[i]]$p.ref=runif(N.stud,0.03,0.05) #### study-specific probability of an event in treatment 1
  data1[[i]]$odds.ref=data1[[i]]$p.ref/(1-data1[[i]]$p.ref)
}

#### define probabilities per treatment, per study arm
for (i in 1:N.sim)
{  
  for(j in 1:N.stud){
    data1[[i]]$odds.t1[j]=data1[[i]]$odds.ref[j]*OR[[i]][data1[[i]]$t1[j]]
    data1[[i]]$odds.t2[j]=data1[[i]]$odds.ref[j]*OR[[i]][data1[[i]]$t2[j]]
    data1[[i]]$p.t1[j]=data1[[i]]$odds.t1[j]/(1+data1[[i]]$odds.t1[j])
    data1[[i]]$p.t2[j]=data1[[i]]$odds.t2[j]/(1+data1[[i]]$odds.t2[j])
  }}


#### generate the data
for (i in 1:N.sim)
{  
  for(j in 1:N.stud){
    data1[[i]]$r1[j]=rbinom(1,data1[[i]]$n1[j],data1[[i]]$p.t1[j]) 
    data1[[i]]$r2[j]=rbinom(1,data1[[i]]$n2[j],data1[[i]]$p.t2[j]) 
  }}
for (i in 1:N.sim){ data1[[i]]=data1[[i]][,-c(6:11)]
data2[[i]]=dlong(data1[[i]])
data3[[i]]=dlong1(data1[[i]])
}
