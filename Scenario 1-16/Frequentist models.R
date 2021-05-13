library(parallel)
library(beepr)
library(netmeta)
library(brglm)
library(lme4)



# Initiate cluster
no_cores <- detectCores()-3
cl <- makeCluster(no_cores)


X1=list()
for (i in 1:N.sim){ X1[[i]]=list("data"=data1[[i]],"logOR"=logOR[[i]],"tau"=tau[[1]]) }

################################################################
###### Mantel Haenszel NMA ###########
MH=function(X)
{
  biasMH<-c()
  coverageMH<-c()
  mseMH=c()
  MH1<-netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "MH", data=X$data)
  MH.res<-data.frame(mean=MH1$TE.fixed[2:NT,1], lowerCI=MH1$lower.fixed[2:NT,1],upperCI=MH1$upper.fixed[2:NT,1],se=MH1$seTE.fixed[2:NT,1])
  biasMH<-c(biasMH, (MH.res$mean-X$logOR))
  MH.res$cover<-(MH.res$lowerCI<X$logOR)&(MH.res$upperCI>X$logOR)
  coverageMH<-c(coverageMH, MH.res$cover)
  mseMH=c(mseMH,biasMH*biasMH+MH.res$se*MH.res$se)
  
  return(list("bias"=biasMH,"cov"=coverageMH,"mseMH"=mseMH))
  
}

clusterExport(cl,"X1")
clusterExport(cl,"NT")
clusterExport(cl,"N.sim")
clusterExport(cl, "MH")
clusterEvalQ(cl, {library(netmeta)})
l2=parLapply(cl,1:N.sim, function(x) MH(X1[[x]]))

biasMH=c()
for (i in 1:N.sim){  biasMH=c(biasMH, l2[[i]]$bias)}
mean(biasMH)

coverageMH=c()
for (i in 1:N.sim){  coverageMH=c(coverageMH, l2[[i]]$cov)}
mean(coverageMH)

mseMH=c()
for (i in 1:N.sim){  mseMH=c(mseMH, l2[[i]]$mseMH)}
mean(mseMH)
################################################################


################################################################
###### NCH NMA ###########
NCH=function(X)
{
  biasNCH<-c()
  coverageNCH<-c()
  mseNCH=c()
  
  NCH1<-netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "NCH", data=X$data)
  NCH.res<-data.frame(mean=NCH1$TE.fixed[2:NT,1], lowerCI=NCH1$lower.fixed[2:NT,1],upperCI=NCH1$upper.fixed[2:NT,1],se=NCH1$seTE.fixed[2:NT,1])
  
  biasNCH<-c(biasNCH, (NCH.res$mean-X$logOR))
  NCH.res$cover<-(NCH.res$lowerCI<X$logOR)&(NCH.res$upperCI>X$logOR)
  coverageNCH<-c(coverageNCH, NCH.res$cover)
  mseNCH=c(mseNCH,biasNCH*biasNCH+NCH.res$se*NCH.res$se)
  
  return(list("bias"=biasNCH,"cov"=coverageNCH,"mseNCH"=mseNCH))
  
}



clusterExport(cl,"X1")
clusterExport(cl,"NT")
clusterExport(cl,"N.sim")
clusterExport(cl, "NCH")
clusterEvalQ(cl, {library(netmeta)})
l3=parLapply(cl,1:N.sim, function(x) NCH(X1[[x]]))


biasNCH=c()
for (i in 1:N.sim){  biasNCH=c(biasNCH, l3[[i]]$bias)}
mean(biasNCH)

coverageNCH=c()
for (i in 1:N.sim){  coverageNCH=c(coverageNCH, l3[[i]]$cov)}
mean(coverageNCH)


mseNCH=c()
for (i in 1:N.sim){  mseNCH=c(mseNCH, l3[[i]]$mseNCH)}



IV=function(X)
{
  biasIV.FE<-c()
  coverageIV.FE<-c()
  biasIV.RE<-c()
  coverageIV.RE<-c()
  mseIV.FE=c()
  mseIV.RE=c()
  bias_tau=c()
  
  IV1<-netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "Inverse", incr=0.5 ,cc.pooled=T, allstudies=T ,sm="OR",data=X$data)
  
  IV.FE.res<-data.frame(mean=IV1$TE.fixed[2:NT,1], lowerCI=IV1$lower.fixed[2:NT,1],upperCI=IV1$upper.fixed[2:NT,1],seFE=IV1$seTE.fixed[2:NT,1])
  biasIV.FE<-c(biasIV.FE, (IV.FE.res$mean-X$logOR))
  IV.FE.res$cover<-(IV.FE.res$lowerCI<X$logOR)&(IV.FE.res$upperCI>X$logOR)
  coverageIV.FE<-c(coverageIV.FE, IV.FE.res$cover)
  mseIV.FE=c(mseIV.FE,biasIV.FE*biasIV.FE+IV.FE.res$seFE*IV.FE.res$seFE)
  
  IV.RE.res<-data.frame(mean=IV1$TE.random[2:NT,1], lowerCI=IV1$lower.random[2:NT,1],upperCI=IV1$upper.random[2:NT,1],seRE=IV1$seTE.random[2:NT,1])
  biasIV.RE<-c(biasIV.RE, (IV.RE.res$mean-X$logOR))
  IV.RE.res$cover<-(IV.RE.res$lowerCI<X$logOR)&(IV.RE.res$upperCI>X$logOR)
  coverageIV.RE<-c(coverageIV.RE, IV.RE.res$cover)
  mseIV.RE=c(mseIV.RE,biasIV.RE*biasIV.RE+IV.RE.res$seRE*IV.RE.res$seRE)
  
  tau1=IV1$tau
  bias_tau<-c(bias_tau,(tau1-X$tau))
  
  return(list("biasFE"=biasIV.FE,"biasRE"=biasIV.RE,"covFE"=coverageIV.FE, "covRE"=coverageIV.RE,"mseIV.FE"=mseIV.FE,"mseIV.RE"=mseIV.RE,"tau"=tau1,"bias_tau"=bias_tau))
}
clusterExport(cl,"X1")
clusterExport(cl,"NT")
clusterExport(cl,"N.sim")
clusterExport(cl, "IV")
clusterEvalQ(cl, {library(netmeta)})
l1=parLapply(cl,1:N.sim, function(x) IV(X1[[x]]))


biasIV.FE=c()
for (i in 1:N.sim){  biasIV.FE=c(biasIV.FE, l1[[i]]$biasFE)}
mean(biasIV.FE)

coverageIV.FE=c()
for (i in 1:N.sim){  coverageIV.FE=c(coverageIV.FE, l1[[i]]$covFE)}
mean(coverageIV.FE)

biasIV.RE=c()
for (i in 1:N.sim){  biasIV.RE=c(biasIV.RE, l1[[i]]$biasRE)}
mean(biasIV.RE)


bias_tau_IV=c()
for (i in 1:N.sim){bias_tau_IV=c(bias_tau_IV, l1[[i]]$bias_tau)}
mean(bias_tau_IV)

coverageIV.RE=c()
for (i in 1:N.sim){  coverageIV.RE=c(coverageIV.RE, l1[[i]]$covRE)}
mean(coverageIV.RE)

mseIV.FE=c()
for (i in 1:N.sim){  mseIV.FE=c(mseIV.FE, l1[[i]]$mseIV.FE)}
mean(mseIV.FE)

mseIV.RE=c()
for (i in 1:N.sim){  mseIV.RE=c(mseIV.RE, l1[[i]]$mseIV.RE)}
mean(mseIV.RE)

tau_IV=c()
for (i in 1:N.sim){tau_IV=c(tau_IV, l1[[i]]$tau)}

stopCluster(cl)





################################################################


################################################################
###### PL-NMA ###########

no_cores1 <- detectCores()-3
cl1 <- makeCluster(no_cores1)

X2=list()
for (i in 1:N.sim){ X2[[i]]=list("data"=data2[[i]],"logOR"=logOR[[i]]) }
FL_F=function(X)
{
  biasFL_F=c()
  coverageFL_F=c()
  coverageFL_Fr=c()
  pcoverageFL_F=c()
  mseFL_F=c()
  mseFL_Fr=c()

  
  prof=function(data.=X,pos=0){
    assign("data.",data.,pos)
    FL_F1=brglm(cbind(events,n-events)~factor(treat)+factor(study),data=data.)
    prof1=profile(FL_F1)
    prof2=confint(prof1)
    prof2=prof2[grep("treat",row.names(prof2)),]
    colnames(prof2)=c("pci.lower","pci.upper")
    prof2=as.data.frame(prof2)
    return(prof2)}
  
  dispersion_Fletcher=function(model){
    require(AICcmodavg)
    phi1=c_hat(model,method = "fletcher")
    if(phi1>1){
      phi=phi1}
    else if((phi1<1)||(phi1=1)){
      phi=1}
    return(phi)}
  
  FL_F=brglm(cbind(events,n-events)~factor(treat)+factor(study),data=X$data)
  p2=prof(data.=X$data,pos=1)
  phi=dispersion_Fletcher(FL_F)
  ests=summary(FL_F)$coefficients
  ests=ests[grep("treat",rownames(ests)),]
  ests=ests[,-c(3,4)]
  ci.lower=ests[,1]-1.96*ests[,2]
  ci.upper=ests[,1]+1.96*ests[,2]
  
  ci.rlower=ests[,1]-1.96*(ests[,2]*sqrt(phi))
  ci.rupper=ests[,1]+1.96*(ests[,2]*sqrt(phi))
  
  n.se=ests[,2]*sqrt(phi)
  
  pci.lower=p2[,1]
  pci.upper=p2[,2]
  ests=cbind(ests,ci.lower,ci.upper,pci.lower,pci.upper)
  ests=as.data.frame(ests)
  
  FL_F.res=data.frame(mean=ests[,1],lowerCI=ci.lower,upperCI=ci.upper,lowerPCI=pci.lower,upperPCI=pci.upper,random_lowerCI=ci.rlower,
                      random_upperCI=ci.rupper)
  
  biasFL_F=c(biasFL_F,(FL_F.res$mean-X$logOR))
  mseFL_F=biasFL_F*biasFL_F+ests[,2]*ests[,2]
  mseFL_Fr=biasFL_F*biasFL_F+n.se*n.se
  FL_F.res$cover=(FL_F.res$lowerCI<X$logOR)&(FL_F.res$upperCI>X$logOR)
  FL_F.res$rcover=((FL_F.res$random_lowerCI<X$logOR)&(FL_F.res$random_upperCI>X$logOR))
  FL_F.res$pcover=(FL_F.res$lowerPCI<X$logOR)&(FL_F.res$upperPCI>X$logOR)
  counter=(phi>1)
  
  coverageFL_F=c(coverageFL_F, FL_F.res$cover)
  coverageFL_Fr=c(coverageFL_Fr,FL_F.res$rcover)
  pcoverageFL_F=c(pcoverageFL_F,FL_F.res$pcover)
  dispersion=phi
  return(list("bias"=biasFL_F,"cov"=coverageFL_F,"rcov"=coverageFL_Fr,"pcov"=pcoverageFL_F,"mse"=mseFL_F,"mse1"=mseFL_Fr,
              "count"=counter,"dispersion"=dispersion))
  
}

clusterExport(cl1,"X2")
clusterExport(cl1,"NT")
clusterExport(cl1,"N.sim")
clusterExport(cl1, "FL_F")
clusterEvalQ(cl1, {library(brglm)})
l8=parLapply(cl1,1:N.sim, function(x) FL_F(X2[[x]]))

biasFL_F=c()
for (i in 1:N.sim){  biasFL_F=c(biasFL_F, l8[[i]]$bias)}
mean(biasFL_F)

coverageFL_F=c()
for (i in 1:N.sim){  coverageFL_F=c(coverageFL_F, l8[[i]]$cov)}
mean(coverageFL_F)


coverageFL_Fr=c()
for (i in 1:N.sim){  coverageFL_Fr=c(coverageFL_Fr, l8[[i]]$rcov)}
mean(coverageFL_Fr)


pcoverageFL_F=c()
for (i in 1:N.sim){  pcoverageFL_F=c(pcoverageFL_F, l8[[i]]$pcov)}
mean(pcoverageFL_F)

mseFL_F=c()
for (i in 1:N.sim){mseFL_F=c(mseFL_F, l8[[i]]$mse)}
mean(mseFL_F)


mseFL_Fr=c()
for (i in 1:N.sim){mseFL_Fr=c(mseFL_Fr, l8[[i]]$mse1)}
mean(mseFL_Fr)


counter=c()
for (i in 1:N.sim){counter=c(counter, l8[[i]]$count)}
sum(counter)


dispersion=c()
for (i in 1:N.sim){dispersion=c(dispersion, l8[[i]]$dispersion)}
phi=dispersion




no_cores2 <- detectCores()-3
cl2 <- makeCluster(no_cores2)


X3=list()
for (i in 1:N.sim){ X3[[i]]=list("data"=data2[[i]],"logOR"=logOR[[i]]) }

################################################################
###### Binomial-Normal ###########
BN=function(X)
{
  biasBN<-c()
  coverageBN<-c()
  mseBN=c()
  
  BN1=glmer(cbind(events,n-events)~factor(treat)+factor(study)+(treat-1|study),family =binomial,data=X$data)
  
  ests=summary(BN1)$coefficients
  ests=ests[grep("treat",rownames(ests)),]
  ests=ests[,-c(3,4)]
  ci.lower=ests[,1]-1.96*ests[,2]
  ci.upper=ests[,1]+1.96*ests[,2]
  BN.res=data.frame(mean=ests[,1],lowerCI=ci.lower,upperCI=ci.upper)
  biasBN<-c(biasBN, (BN.res$mean-X$logOR))
  BN.res$cover<-(BN.res$lowerCI<X$logOR)&(BN.res$upperCI>X$logOR)
  coverageBN<-c(coverageBN, BN.res$cover)
  mseBN=c(mseBN,biasBN*biasBN+ests[,2]*ests[,2])
  
  heterogeneity=data.frame(VarCorr(BN1))[,4]
  counterBN=(heterogeneity>0)
  
  return(list("bias"=biasBN,"cov"=coverageBN,"mseBN"=mseBN,"counterBN"=counterBN,"heterogeneity"=heterogeneity))
  
}

clusterExport(cl2,"X3")
clusterExport(cl2,"NT")
clusterExport(cl2,"N.sim")
clusterExport(cl2, "BN")
clusterEvalQ(cl2, {library(lme4)})
l9=parLapply(cl2,1:N.sim, function(x) BN(X3[[x]]))

biasBN=c()
for (i in 1:N.sim){  biasBN=c(biasBN, l9[[i]]$bias)}
mean(biasBN)

coverageBN=c()
for (i in 1:N.sim){  coverageBN=c(coverageBN, l9[[i]]$cov)}
mean(coverageBN)


mseBN=c()
for (i in 1:N.sim){  mseBN=c(mseBN, l9[[i]]$mseBN)}
mean(mseBN)


counterBN=c()
for (i in 1:N.sim){counterBN=c(counterBN, l9[[i]]$counterBN)}
sum(counterBN)


heterogeneity=c()
for (i in 1:N.sim){heterogeneity=c(heterogeneity, l9[[i]]$heterogeneity)}
tau2_BN=heterogeneity



#############################

################################### RESULTS ##################################################


round(mean(biasIV.FE), digits=2)
round(mean(abs(biasIV.FE)), digits=2)
round(mean(coverageIV.FE), digits=3)*100
round(mean(mseIV.FE),digits = 2)
# IV.RE
round(mean(biasIV.RE), digits=2)
round(mean(abs(biasIV.RE)), digits=2)
round(mean(coverageIV.RE), digits=3)*100
round(mean(mseIV.RE),digits = 2)

avg_bias_tau_IV=c(round(mean(bias_tau_IV), digits=2))
avg_bias_tau_IV=data.frame(avg_bias_tau_IV)
names(avg_bias_tau_IV)=c("mean_bias_tau")
# MH
round(mean(biasMH), digits=2)
round(mean(abs(biasMH)), digits=2)
round(mean(coverageMH), digits=3)*100
round(mean(mseMH),digits = 2)
# NCH
round(mean(biasNCH), digits=2)
round(mean(abs(biasNCH)), digits=2)
round(mean(coverageNCH), digits=3)*100
round(mean(mseNCH),digits = 2)
# FL_F
round(mean(biasFL_F), digits=2)
round(mean(abs(biasFL_F)), digits=2)
round(mean(coverageFL_F), digits=3)*100
round(mean(coverageFL_Fr), digits=3)*100
round(mean(pcoverageFL_F), digits=3)*100
round(mean(mseFL_F),digits = 2)
round(mean(mseFL_Fr),digits = 2)
#BN
round(mean(biasBN), digits=2)
round(mean(abs(biasBN)), digits=2)
round(mean(coverageBN), digits=3)*100
round(mean(mseBN),digits = 2)





bias=c(round(mean(biasIV.FE), digits=2),round(mean(biasIV.RE), digits=2),round(mean(biasMH), digits=2),round(mean(biasNCH), digits=2),
       round(mean(biasFL_F), digits=2),round(mean(biasFL_F), digits=2),round(mean(biasBN), digits=2))





absolute_bias=c(round(mean(abs(biasIV.FE)), digits=2),round(mean(abs(biasIV.RE)), digits=2),round(mean(abs(biasMH)), digits=2),
                round(mean(abs(biasNCH)), digits=2),round(mean(abs(biasFL_F)), digits=2),round(mean(abs(biasFL_F)), digits=2)
                ,round(mean(abs(biasBN)), digits=2))

coverage=c(round(mean(coverageIV.FE), digits=3)*100,round(mean(coverageIV.RE), digits=3)*100,round(mean(coverageMH), digits=3)*100,
           round(mean(coverageNCH), digits=3)*100,round(mean(coverageFL_F), digits=3)*100,round(mean(coverageFL_Fr), digits=3)*100,
           round(mean(coverageBN), digits=3)*100)

coverage1=c(round(mean(coverageIV.FE), digits=3)*100,round(mean(coverageIV.RE), digits=3)*100,round(mean(coverageMH), digits=3)*100,
            round(mean(coverageNCH), digits=3)*100,round(mean(pcoverageFL_F), digits=3)*100,round(mean(coverageFL_Fr), digits=3)*100,
            round(mean(coverageBN), digits=3)*100)

mse=c(round(mean(mseIV.FE),digits = 2),round(mean(mseIV.RE),digits = 2),round(mean(mseMH),digits = 2),round(mean(mseNCH),digits = 2),
      round(mean(mseFL_F),digits = 2),round(mean(mseFL_Fr),digits = 2),round(mean(mseBN),digits = 2))

results=cbind(bias,absolute_bias,mse,coverage)
results=data.frame(results)
row.names(results)=c("IV_FE","IV_RE","MH","NCH","FL_F","FL_F(R)","B-N")
results

results1=cbind(bias,absolute_bias,mse,coverage1)
results1=data.frame(results1)
row.names(results1)=c("IV_FE","IV_RE","MH","NCH","FL_F(P)","FL_F(R)","B-N")
results1


het1=cbind(counter,dispersion)
het1=as.data.frame(het1)
names(het1)=c("counter_Fletcher","dispersion")
het2=cbind(counterBN,tau2_BN)
het2=as.data.frame(het2)

het=data.frame(het1,het2)


tau_info_IV=data.frame(bias_tau_IV,tau_IV)
names(tau_info_IV)=c("bias_tau_IV","tau_IV")


