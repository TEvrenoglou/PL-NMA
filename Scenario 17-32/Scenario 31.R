############ The following data are being generated using the netmeta version 1.2-1 and meta version 4.15-1. ####################
############ If you face any issue on reproducing the datasets please make sure that you have those versions installed. ####################

library(dplyr)
library(MASS)
library(tidyverse)
library(netmeta)






set.seed(44)
N.sim=1000


data1=list()
data2=list()
data3=list()
data4=list()
logOR=list()
logOR1=list()
OR=list()
NT=5 #### number of treatments in the network
NS=8 #### number of studies per comparison
tau=0 ####  heterogeneity SD
Npmin=100 #### minimum number of patients per arm
Npmax=200 #### maximum number of patients per arm

### define treatment indices
t1=c(1:NT)
t=(rep(t1,NS))
stud=rep(1:NS, each=NT)
for (i in 1:N.sim)
{ 
  data1[[i]]=data.frame(stud,t)
}
### define patients per treatment arm
for (i in 1:N.sim)
{   
  data1[[i]]$n=rep(round(runif(NS,Npmin,Npmax)),each=NT)
}

narms=length(data1[[1]]$stud)


### define logOR
for (i in 1:N.sim)
{   
  logOR[[i]]=seq(from =1/(NT-1), to = 1, by = 1/(NT-1))  ### true logOR across studies 
  OR[[i]]=c(1,exp(logOR[[i]]))
}


#### define odds per treatment, per study arm
for (i in 1:N.sim)
{  
  data1[[i]]$p.ref=rep(runif(NS,0.005,0.1),each=NT)
  data1[[i]]$odds.ref=data1[[i]]$p.ref/(1-data1[[i]]$p.ref)
}

Sigma=matrix(c(tau^2,tau^2/2,tau^2/2,tau^2/2,tau^2/2,
               tau^2/2,tau^2,tau^2/2,tau^2/2,tau^2/2,
               tau^2/2,tau^2/2,tau^2,tau^2/2,tau^2/2,
               tau^2/2,tau^2/2,tau^2/2,tau^2,tau^2/2,
               tau^2/2,tau^2/2,tau^2/2,tau^2/2,tau^2),
             nrow=5)

for (i in 1:N.sim)
{   
  logOR1[[i]]=c(0,logOR[[i]])
  for(j in 1:narms){
    
    
    data1[[i]]$truelogOR.t[j]=logOR1[[i]][data1[[i]]$t[j]]
    #data1[[i]]$truelogOR.t2[j]=logOR1[[i]][data1[[i]]$t2[j]]
    
    test1=mvrnorm(1,c(0,0.25,0.50,0.75,1),Sigma)
    test2=rnorm(1,0.25,tau)
    test3=rnorm(1,0.5,tau)
    test4=rnorm(1,0.75,tau)
    test5=rnorm(1,1,tau)
    
    data1[[i]]$st.sp.logOR.t1[j]=test1[1]*(data1[[i]]$t[j]!=1)
    data1[[i]]$st.sp.logOR.t2[j]=test1[2]*(data1[[i]]$t[j]!=1)+test2*(data1[[i]]$t[j]==1)
    data1[[i]]$st.sp.logOR.t3[j]=test1[3]*(data1[[i]]$t[j]!=1)+test3*(data1[[i]]$t[j]==1)
    data1[[i]]$st.sp.logOR.t4[j]=test1[4]*(data1[[i]]$t[j]!=1)+test4*(data1[[i]]$t[j]==1)
    data1[[i]]$st.sp.logOR.t5[j]=test1[5]*(data1[[i]]$t[j]!=1)+test5*(data1[[i]]$t[j]==1)
    
    data1[[i]]$odds.t1[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t1[j])
    data1[[i]]$odds.t2[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t2[j])
    data1[[i]]$odds.t3[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t3[j])
    data1[[i]]$odds.t4[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t4[j])
    data1[[i]]$odds.t5[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t5[j])
    
    data1[[i]]$p.t1[j]=data1[[i]]$odds.t1[j]/(1+data1[[i]]$odds.t1[j])
    data1[[i]]$p.t2[j]=data1[[i]]$odds.t2[j]/(1+data1[[i]]$odds.t2[j])
    data1[[i]]$p.t3[j]=data1[[i]]$odds.t3[j]/(1+data1[[i]]$odds.t3[j])
    data1[[i]]$p.t4[j]=data1[[i]]$odds.t4[j]/(1+data1[[i]]$odds.t4[j])
    data1[[i]]$p.t5[j]=data1[[i]]$odds.t5[j]/(1+data1[[i]]$odds.t5[j])
  }}


#### generate the events
for (i in 1:N.sim){
  for(j in 1:narms){
    data1[[i]]$r1[j]=-100
    data1[[i]]$r1[j]=rbinom(1,data1[[i]]$n[j],data1[[i]]$p.ref[j]) *(data1[[i]]$t[j]==1)
  }
}


for (i in 1:N.sim)
{  
  for(j in 1:narms){
    data1[[i]]$r2[j]=rbinom(1,data1[[i]]$n[j],data1[[i]]$p.t2[j]) *(data1[[i]]$t[j]==2)
    data1[[i]]$r3[j]=rbinom(1,data1[[i]]$n[j],data1[[i]]$p.t3[j]) *(data1[[i]]$t[j]==3)
    data1[[i]]$r4[j]=rbinom(1,data1[[i]]$n[j],data1[[i]]$p.t4[j]) *(data1[[i]]$t[j]==4)
    data1[[i]]$r5[j]=rbinom(1,data1[[i]]$n[j],data1[[i]]$p.t5[j]) *(data1[[i]]$t[j]==5)
  }}

for (i in 1:N.sim){   data1[[i]]$r=data1[[i]]$r1+data1[[i]]$r2+data1[[i]]$r3+data1[[i]]$r4+data1[[i]]$r5}
for (i in 1:N.sim){ data1[[i]]=data1[[i]][,-c(4:26)]}
data2=data1
data4=data1
for (i in 1:N.sim){ names(data2[[i]])=c("study","treat","n","events")}
for (i in 1:N.sim){ names(data4[[i]])=c("study","treatment","sampleSize","responders")}
for(i in 1:N.sim){
  data3[[i]] =pairwise(treat=t, event=r, n=n, studlab=stud, data=data1[[i]])
}
data1=data3
for (i in 1:N.sim){ data1[[i]]=data1[[i]][,-c(2:6)]}
for (i in 1:N.sim){ data1[[i]]=data1[[i]][,-c(3,5,6,7,8)]

}
##################
