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
  seMH=c()
  length_MH=c()
  lower_MH=c()
  upper_MH=c()
  estimate_MH=c()
  
  
  MH1<-netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "MH", data=X$data)
  MH.res<-data.frame(mean=MH1$TE.fixed[2:NT,1], lowerCI=MH1$lower.fixed[2:NT,1],upperCI=MH1$upper.fixed[2:NT,1],se=MH1$seTE.fixed[2:NT,1])
  
  biasMH<-c(biasMH, (MH.res$mean-X$logOR))
  MH.res$cover<-(MH.res$lowerCI<X$logOR)&(MH.res$upperCI>X$logOR)
  coverageMH<-c(coverageMH, MH.res$cover)
  mseMH=c(mseMH,(biasMH*biasMH)+(MH.res$se*MH.res$se))
  
  seMH=c(seMH,MH.res$se)
  length_MH=c(length_MH,(MH.res$upperCI-MH.res$lowerCI))
  
  lower_MH=c(lower_MH,MH.res$lowerCI)
  upper_MH=c(upper_MH,MH.res$upperCI)
  estimate_MH=c(estimate_MH,MH.res$mean)
  
  return(list("bias"=biasMH,"cov"=coverageMH,"mseMH"=mseMH,
              "seMH"=seMH,"length_MH"=length_MH,
              "lower_MH"=lower_MH,
              "upper_MH"=upper_MH,
              "estimate_MH"=estimate_MH
              
  ))
  
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
#beep(sound=3)

mseMH=c()
for (i in 1:N.sim){  mseMH=c(mseMH, l2[[i]]$mseMH)}
mean(mseMH)

seMH=c()
for (i in 1:N.sim){  seMH=c(seMH, l2[[i]]$seMH)}
mean(seMH)

length_MH=c()
for (i in 1:N.sim){  length_MH=c(length_MH, l2[[i]]$length_MH)}
mean(length_MH)

upper_MH=c()
for (i in 1:N.sim){  upper_MH=c(upper_MH, l2[[i]]$upper_MH)}

lower_MH=c()
for (i in 1:N.sim){  lower_MH=c(lower_MH, l2[[i]]$lower_MH)}


estimate_MH=c()
for (i in 1:N.sim){  estimate_MH=c(estimate_MH, l2[[i]]$estimate_MH)}



#beep(sound=3)
################################################################


################################################################
###### NCH NMA ###########
NCH=function(X)
{
  biasNCH<-c()
  coverageNCH<-c()
  mseNCH=c()
  se_NCH=c()
  length_NCH=c()
  
  upper_NCH=c()
  lower_NCH=c()
  estimate_NCH=c()
  
  
  NCH1<-netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "NCH", data=X$data)
  NCH.res<-data.frame(mean=NCH1$TE.fixed[2:NT,1], lowerCI=NCH1$lower.fixed[2:NT,1],upperCI=NCH1$upper.fixed[2:NT,1],se=NCH1$seTE.fixed[2:NT,1])
  
  
  biasNCH<-c(biasNCH, (NCH.res$mean-X$logOR))
  NCH.res$cover<-(NCH.res$lowerCI<X$logOR)&(NCH.res$upperCI>X$logOR)
  coverageNCH<-c(coverageNCH, NCH.res$cover)
  mseNCH=c(mseNCH,((biasNCH*biasNCH)+(NCH.res$se*NCH.res$se)))
  
  se_NCH=c(se_NCH,NCH.res$se)
  
  length_NCH=c(length_NCH,NCH.res$upperCI-NCH.res$lowerCI)
  
  upper_NCH=c(upper_NCH,NCH.res$upperCI)
  
  lower_NCH=c(lower_NCH,NCH.res$lowerCI)
  
  estimate_NCH=c(estimate_NCH,NCH.res$mean)
  
  return(list("bias"=biasNCH,"cov"=coverageNCH,"mseNCH"=mseNCH,
              "se_NCH"=se_NCH,"length_NCH"=length_NCH,
              "upper_NCH"=upper_NCH,
              "lower_NCH"=lower_NCH,
              "estimate_NCH"=estimate_NCH
              
  ))
  
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
#beep(sound=3)

mseNCH=c()
for (i in 1:N.sim){  mseNCH=c(mseNCH, l3[[i]]$mseNCH)}
mean(mseNCH)

se_NCH=c()
for (i in 1:N.sim){  se_NCH=c(se_NCH, l3[[i]]$se_NCH)}
mean(se_NCH)

length_NCH=c()
for (i in 1:N.sim){  length_NCH=c(length_NCH, l3[[i]]$length_NCH)}
mean(length_NCH)

upper_NCH=c()
for (i in 1:N.sim){  upper_NCH=c(upper_NCH, l3[[i]]$upper_NCH)}


lower_NCH=c()
for (i in 1:N.sim){  lower_NCH=c(lower_NCH, l3[[i]]$lower_NCH)}


estimate_NCH=c()
for (i in 1:N.sim){  estimate_NCH=c(estimate_NCH, l3[[i]]$estimate_NCH)}

################################################################


################################################################

IV=function(X)
{
  biasIV.FE<-c()
  coverageIV.FE<-c()
  biasIV.RE<-c()
  coverageIV.RE<-c()
  mseIV.FE=c()
  mseIV.RE=c()
  bias_tau=c()
  
  se_IV.FE=c()
  se_IV.RE=c()
  
  
  length_IV.FE=c()
  length_IV.RE=c()
  
  upper_IV.FE=c()
  upper_IV.RE=c()
  
  lower_IV.FE=c()
  lower_IV.RE=c()
  
  estimate_IV.FE=c()
  estimate_IV.RE=c()
  
  
  IV1<-netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "Inverse", incr=0.5 ,cc.pooled=T, allstudies=T ,sm="OR",data=X$data)
  
  IV.FE.res<-data.frame(mean=IV1$TE.fixed[2:NT,1], lowerCI=IV1$lower.fixed[2:NT,1],upperCI=IV1$upper.fixed[2:NT,1],seFE=IV1$seTE.fixed[2:NT,1])
  biasIV.FE<-c(biasIV.FE, (IV.FE.res$mean-X$logOR))
  IV.FE.res$cover<-(IV.FE.res$lowerCI<X$logOR)&(IV.FE.res$upperCI>X$logOR)
  coverageIV.FE<-c(coverageIV.FE, IV.FE.res$cover)
  mseIV.FE=c(mseIV.FE,biasIV.FE*biasIV.FE+IV.FE.res$seFE*IV.FE.res$seFE)
  
  se_IV.FE=c(se_IV.FE,IV.FE.res$seFE)
  length_IV.FE=c(length_IV.FE,(IV.FE.res$upperCI-IV.FE.res$lowerCI))
  
  lower_IV.FE=c(lower_IV.FE,IV.FE.res$lowerCI)
  upper_IV.FE=c(upper_IV.FE,IV.FE.res$upperCI)
  estimate_IV.FE=c(estimate_IV.FE,IV.FE.res$mean)
  
  IV.RE.res<-data.frame(mean=IV1$TE.random[2:NT,1], lowerCI=IV1$lower.random[2:NT,1],upperCI=IV1$upper.random[2:NT,1],seRE=IV1$seTE.random[2:NT,1])
  biasIV.RE<-c(biasIV.RE, (IV.RE.res$mean-X$logOR))
  IV.RE.res$cover<-(IV.RE.res$lowerCI<X$logOR)&(IV.RE.res$upperCI>X$logOR)
  coverageIV.RE<-c(coverageIV.RE, IV.RE.res$cover)
  mseIV.RE=c(mseIV.RE,biasIV.RE*biasIV.RE+IV.RE.res$seRE*IV.RE.res$seRE)
  
  se_IV.RE=c(se_IV.RE,IV.RE.res$seRE)
  length_IV.RE=c(length_IV.RE,(IV.RE.res$upperCI-IV.RE.res$lowerCI))
  
  
  lower_IV.RE=c(lower_IV.RE,IV.RE.res$lowerCI)
  upper_IV.RE=c(upper_IV.RE,IV.RE.res$upperCI)
  estimate_IV.RE=c(estimate_IV.RE,IV.RE.res$mean)
  
  tau1=IV1$tau
  bias_tau<-c(bias_tau,(tau1-X$tau))
  
  return(list("biasFE"=biasIV.FE,"biasRE"=biasIV.RE,"covFE"=coverageIV.FE, "covRE"=coverageIV.RE,
              "mseIV.FE"=mseIV.FE,"mseIV.RE"=mseIV.RE,"tau"=tau1,"bias_tau"=bias_tau,
              "se_IV_FE"=se_IV.FE,"length_IV_FE"=length_IV.FE,
              "lower_IV_FE"=lower_IV.FE,"upper_IV_FE"=upper_IV.FE,
              "estimate_IV_FE"=estimate_IV.FE,
              
              "se_IV_RE"=se_IV.RE, "length_IV_RE"=length_IV.RE,
              "lower_IV_RE"=lower_IV.RE,"upper_IV_RE"=upper_IV.RE,
              "estimate_IV_RE"=estimate_IV.RE
              
              
  ))
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

length_IV_FE=c()
for (i in 1:N.sim){length_IV_FE=c(length_IV_FE, l1[[i]]$length_IV_FE)}
mean(length_IV_FE)


length_IV_RE=c()
for (i in 1:N.sim){length_IV_RE=c(length_IV_RE, l1[[i]]$length_IV_RE)} 
mean(length_IV_RE)

se_IV_FE=c()
for (i in 1:N.sim){se_IV_FE=c(se_IV_FE, l1[[i]]$se_IV_FE)}

se_IV_RE=c()
for (i in 1:N.sim){se_IV_RE=c(se_IV_RE, l1[[i]]$se_IV_RE)}

lower_IV_FE=c()
for (i in 1:N.sim){lower_IV_FE=c(lower_IV_FE, l1[[i]]$lower_IV_FE)}

lower_IV_RE=c()
for (i in 1:N.sim){lower_IV_RE=c(lower_IV_RE, l1[[i]]$lower_IV_RE)}


upper_IV_FE=c()
for (i in 1:N.sim){upper_IV_FE=c(upper_IV_FE, l1[[i]]$upper_IV_FE)}

upper_IV_RE=c()
for (i in 1:N.sim){upper_IV_RE=c(upper_IV_RE, l1[[i]]$upper_IV_RE)}

estimate_IV_FE=c()
for (i in 1:N.sim){estimate_IV_FE=c(estimate_IV_FE, l1[[i]]$estimate_IV_FE)}

estimate_IV_RE=c()
for (i in 1:N.sim){estimate_IV_RE=c(estimate_IV_RE, l1[[i]]$estimate_IV_RE)}



stopCluster(cl)



rm(MH)
rm(NCH)
rm(IV)
rm(X1)
rm(data1)
rm(l1)
rm(l2)
rm(l3)

################################################################


################################################################
###### PL-NMA ###########

no_cores1 <- detectCores()-2
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
  
  se_fixed=c()
  se_random=c()
  
  estimatef=c()
  
  upperf=c()
  lowerf=c()
  
  upperr=c()
  lowerr=c()
  
  upperp=c()
  lowerp=c()
  
  lengthf=c()
  lengthr=c()
  lengthp=c()
  
  #counter=c()
  #dispersion=c()
  
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
  #ests[,2]=ests[,2]*sqrt(phi)
  ci.lower=ests[,1]-1.96*ests[,2]
  ci.upper=ests[,1]+1.96*ests[,2]
  
  ci.rlower=ests[,1]-1.96*(ests[,2]*sqrt(phi))
  ci.rupper=ests[,1]+1.96*(ests[,2]*sqrt(phi))
  
  n.se=ests[,2]*sqrt(phi)
  n.se=unname(n.se)
  
  pci.lower=p2[,1]
  pci.upper=p2[,2]
  ests=cbind(ests,ci.lower,ci.upper,pci.lower,pci.upper)
  ests=as.data.frame(ests)
  
  FL_F.res=data.frame(mean=ests[,1],lowerCI=ci.lower,upperCI=ci.upper,lowerPCI=pci.lower,upperPCI=pci.upper,random_lowerCI=ci.rlower,
                      random_upperCI=ci.rupper,se=ests[,2],se_random=n.se)
  
  biasFL_F=c(biasFL_F,(FL_F.res$mean-X$logOR))
  mseFL_F=c(mseFL_F,(biasFL_F*biasFL_F)+(ests[,2]*ests[,2]))
  mseFL_Fr=c(mseFL_Fr,(biasFL_F*biasFL_F)+(n.se*n.se))
  FL_F.res$cover=(FL_F.res$lowerCI<X$logOR)&(FL_F.res$upperCI>X$logOR)
  FL_F.res$rcover=((FL_F.res$random_lowerCI<X$logOR)&(FL_F.res$random_upperCI>X$logOR))
  FL_F.res$pcover=(FL_F.res$lowerPCI<X$logOR)&(FL_F.res$upperPCI>X$logOR)
  counter=(phi>1)
  
  coverageFL_F=c(coverageFL_F, FL_F.res$cover)
  coverageFL_Fr=c(coverageFL_Fr,FL_F.res$rcover)
  pcoverageFL_F=c(pcoverageFL_F,FL_F.res$pcover)
  dispersion=phi
  
  estimatef=c(estimatef,FL_F.res$mean)
  
  se_fixed=c(se_fixed,FL_F.res$se)
  se_random=c(se_random,FL_F.res$se_random)
  
  upperf=c(upperf,FL_F.res$upperCI)
  lowerf=c(lowerf,FL_F.res$lowerCI)
  
  lengthf=c(lengthf,FL_F.res$upperCI-FL_F.res$lowerCI)
  
  
  upperp=c(upperp,FL_F.res$upperPCI)
  lowerp=c(lowerp,FL_F.res$lowerPCI)
  
  lengthp=c(lengthp,FL_F.res$upperPCI-FL_F.res$lowerPCI)
  
  upperr=c(upperr,FL_F.res$random_upperCI)
  lowerr=c(lowerr,FL_F.res$random_lowerCI)
  
  lengthr=c(lengthr,FL_F.res$random_upperCI-FL_F.res$random_lowerCI)
  
  
  return(list("bias"=biasFL_F,"cov"=coverageFL_F,"rcov"=coverageFL_Fr,"pcov"=pcoverageFL_F,"mse"=mseFL_F,"mse1"=mseFL_Fr,
              "count"=counter,"dispersion"=dispersion,
              "estimatef"=estimatef,
              "se_fixed"=se_fixed,"se_random"=se_random,"upperf"=upperf,"lowerf"=lowerf,"lengthf"=lengthf,
              "upperp"=upperp,"lowerp"=lowerp,"lengthp"=lengthp, "upperr"=upperr,"lowerr"=lowerr,"lengthr"=lengthr
              
  ))
  
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
#beep(sound=3)

coverageFL_Fr=c()
for (i in 1:N.sim){  coverageFL_Fr=c(coverageFL_Fr, l8[[i]]$rcov)}
mean(coverageFL_Fr)
#beep(sound=3)

pcoverageFL_F=c()
for (i in 1:N.sim){  pcoverageFL_F=c(pcoverageFL_F, l8[[i]]$pcov)}
mean(pcoverageFL_F)
#beep(sound=3)

mseFL_F=c()
for (i in 1:N.sim){mseFL_F=c(mseFL_F, l8[[i]]$mse)}
mean(mseFL_F)
#beep(sound=3)

mseFL_Fr=c()
for (i in 1:N.sim){mseFL_Fr=c(mseFL_Fr, l8[[i]]$mse1)}
mean(mseFL_Fr)
#beep(sound=3)

counter=c()
for (i in 1:N.sim){counter=c(counter, l8[[i]]$count)}
sum(counter)
#beep(sound=3)

dispersion=c()
for (i in 1:N.sim){dispersion=c(dispersion, l8[[i]]$dispersion)}
phi=dispersion

estimatef=c()
for (i in 1:N.sim){estimatef=c(estimatef, l8[[i]]$estimatef)}

se_fixed=c()
for (i in 1:N.sim){se_fixed=c(se_fixed, l8[[i]]$se_fixed)}

se_random=c()
for (i in 1:N.sim){se_random=c(se_random, l8[[i]]$se_random)}

lower_fixed=c()
for (i in 1:N.sim){lower_fixed=c(lower_fixed, l8[[i]]$lowerf)}

upper_fixed=c()
for (i in 1:N.sim){upper_fixed=c(upper_fixed, l8[[i]]$upperf)}

length_fixed=c()
for (i in 1:N.sim){length_fixed=c(length_fixed, l8[[i]]$lengthf)}
mean(length_fixed)

lower_random=c()
for (i in 1:N.sim){lower_random=c(lower_random, l8[[i]]$lowerr)}

upper_random=c()
for (i in 1:N.sim){upper_random=c(upper_random, l8[[i]]$upperr)}


length_random=c()
for (i in 1:N.sim){length_random=c(length_random, l8[[i]]$lengthr)}
mean(length_random)


lower_profile=c()
for (i in 1:N.sim){lower_profile=c(lower_profile, l8[[i]]$lowerp)}


upper_profile=c()
for (i in 1:N.sim){upper_profile=c(upper_profile, l8[[i]]$upperp)}

length_profile=c()
for (i in 1:N.sim){length_profile=c(length_profile, l8[[i]]$lengthp)}
mean(length_profile)

#beep(sound=3)

########################################################
########################################################


rm(FL_F)
rm(l8)
rm(X2)
#### Binomial-Normal model

no_cores2 <- detectCores()#-3
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
  se_BN=c()
  length_BN=c()
  lower_BN=c()
  upper_BN=c()
  estimate_BN=c()
  count_fail_fixed=c()
  
  biasBN_random<-c()
  coverageBN_random<-c()
  mseBN_random=c()
  se_BN_random=c()
  length_BN_random=c()
  lower_BN_random=c()
  upper_BN_random=c()
  estimate_BN_random=c()
  count_fail_random=c()
  
  biasBN_glm<-c()
  coverageBN_glm<-c()
  mseBN_glm=c()
  se_BN_glm=c()
  length_BN_glm=c()
  lower_glm=c()
  upper_glm=c()
  estimate_glm=c()
  
  BN1=glmer(cbind(events,n-events)~factor(treat)+factor(study)+(treat-1|study),family =binomial,data=X$data)
  
  BN2=glmer(cbind(events,n-events)~factor(treat)+(treat-1|study)+(1|study),family =binomial,data=X$data)
  
  BN3=glm(cbind(events,n-events)~factor(treat)+factor(study),family =binomial,data=X$data)
  
  ests=summary(BN1)$coefficients
  ests=ests[grep("treat",rownames(ests)),]
  ests=ests[,-c(3,4)]
  #ests[,2]=ests[,2]*sqrt(phi)
  ci.lower=ests[,1]-1.96*ests[,2]
  ci.upper=ests[,1]+1.96*ests[,2]
  BN.res=data.frame(mean=ests[,1],lowerCI=ci.lower,upperCI=ci.upper)
  biasBN<-c(biasBN, (BN.res$mean-X$logOR))
  BN.res$cover<-(BN.res$lowerCI<X$logOR)&(BN.res$upperCI>X$logOR)
  coverageBN<-c(coverageBN, BN.res$cover)
  mseBN=c(mseBN,biasBN*biasBN+ests[,2]*ests[,2])
  
  se_BN=c(se_BN,ests[,2])
  
  length_BN=c(length_BN,ci.upper-ci.lower)
  
  heterogeneity=data.frame(VarCorr(BN1))[,4]
  counterBN=(heterogeneity>0)
  
  lower_BN=c(lower_BN,BN.res$lowerCI)
  upper_BN=c(upper_BN,BN.res$upperCI)
  estimate_BN=c(estimate_BN,BN.res$mean)
  
  
  m_fixed <- tryCatch({
    withCallingHandlers({
      error <- FALSE
      list(model = glmer(cbind(events,n-events)~factor(treat)+factor(study)+(treat-1|study),family =binomial,data=X$data),
           error = error)
    },warning = function(w) {
      if(grepl('failed to converge', w$message)) error <<- TRUE
    }
    )})
  
  count_fail_fixed=c(count_fail_fixed,m_fixed$error)
  
  
  
  
  
  #####################################################################
  
  
  ests_random=summary(BN2)$coefficients
  ests_random=ests_random[grep("treat",rownames(ests_random)),]
  ests_random=ests_random[,-c(3,4)]
  #ests[,2]=ests[,2]*sqrt(phi)
  ci.lower_random=ests_random[,1]-1.96*ests_random[,2]
  ci.upper_random=ests_random[,1]+1.96*ests_random[,2]
  
  BN.res_random=data.frame(mean=ests_random[,1],lowerCI=ci.lower_random,upperCI=ci.upper_random)
  
  biasBN_random<-c(biasBN_random, (BN.res_random$mean-X$logOR))
  BN.res_random$cover<-(BN.res_random$lowerCI<X$logOR)&(BN.res_random$upperCI>X$logOR)
  coverageBN_random<-c(coverageBN_random, BN.res_random$cover)
  mseBN_random=c(mseBN_random,biasBN_random*biasBN_random+ests_random[,2]*ests_random[,2])
  se_BN_random=c(se_BN_random,ests_random[,2])
  length_BN_random=c(length_BN_random,ci.upper_random-ci.lower_random)
  
  heterogeneity_random=data.frame(VarCorr(BN2))[,4]
  counterBN_random=(heterogeneity>0)
  
  lower_BN_random=c(lower_BN_random,BN.res_random$lowerCI)
  upper_BN_random=c(upper_BN_random,BN.res_random$upperCI)
  estimate_BN_random=c(estimate_BN_random,BN.res_random$mean)
  
  heterogeneity_random=VarCorr(BN2)
  C=as.data.frame(heterogeneity_random)
  
  C=subset(C,C$var1=="treat")
  heterogeneity_random=C$sdcor
  
  
  m_random <- tryCatch({
    withCallingHandlers({
      error <- FALSE
      list(model = glmer(cbind(events,n-events)~factor(treat)+(treat-1|study)+(1|study),family =binomial,data=X$data),
           error = error)
    },warning = function(w) {
      if(grepl('failed to converge', w$message)) error <<- TRUE
    }
    )})
  
  count_fail_random=c(count_fail_random,m_random$error)
  
  
  
  #####################################################################
  
  
  ests_glm=summary(BN3)$coefficients
  ests_glm=ests_glm[grep("treat",rownames(ests_glm)),]
  ests_glm=ests_glm[,-c(3,4)]
  #ests[,2]=ests[,2]*sqrt(phi)
  ci.lower_glm=ests_glm[,1]-1.96*ests_glm[,2]
  ci.upper_glm=ests_glm[,1]+1.96*ests_glm[,2]
  
  BN.res_glm=data.frame(mean=ests_glm[,1],lowerCI=ci.lower_glm,upperCI=ci.upper_glm)
  
  biasBN_glm<-c(biasBN_glm, (BN.res_glm$mean-X$logOR))
  BN.res_glm$cover<-(BN.res_glm$lowerCI<X$logOR)&(BN.res_glm$upperCI>X$logOR)
  coverageBN_glm<-c(coverageBN_glm, BN.res_glm$cover)
  mseBN_glm=c(mseBN_glm,biasBN_glm*biasBN_glm+ests_glm[,2]*ests_glm[,2])
  se_BN_glm=c(se_BN_glm,ests_glm[,2])
  length_BN_glm=c(length_BN_glm,ci.upper_glm-ci.lower_glm)
  
  lower_glm=c(lower_glm,BN.res_glm$lowerCI)
  upper_glm=c(upper_glm,BN.res_glm$upperCI)
  estimate_glm=c(estimate_glm,BN.res_glm$mean)
  
  #####################################################################
  
  return(list("bias"=biasBN,"cov"=coverageBN,"mseBN"=mseBN,"counterBN"=counterBN,"heterogeneity"=heterogeneity,
              "se_BN"=se_BN,"length_BN"=length_BN,
              "count_fail_fixed"=count_fail_fixed,
              "bias_random"=biasBN_random,"cov_random"=coverageBN_random,"mseBN_random"=mseBN_random,
              "counterBN_random"=counterBN_random,"heterogeneity_random"=heterogeneity_random,
              "se_BN_random"=se_BN_random,
              "length_BN_random"=length_BN_random,
              "count_fail_random"=count_fail_random,
              "bias_glm"=biasBN_glm,"cov_glm"=coverageBN_glm,"mse_glm"=mseBN_glm,
              "se_glm"=se_BN_glm,"length_glm"=length_BN_glm,
              "lower_BN"=lower_BN,"upper_BN"=upper_BN,"estimate_BN"=estimate_BN,
              "lower_BN_random"=lower_BN_random,"upper_BN_random"=upper_BN_random,"estimate_BN_random"=estimate_BN_random,
              "lower_glm"=lower_glm,"upper_glm"=upper_glm,"estimate_glm"=estimate_glm
              
  )
  )
  
}

clusterExport(cl2,"X3")
clusterExport(cl2,"NT")
clusterExport(cl2,"N.sim")
clusterExport(cl2, "BN")
clusterEvalQ(cl2, {library(lme4)})
l9=parLapply(cl2,1:N.sim, function(x) BN(X3[[x]]))

biasBN=c()
for (i in 1:N.sim){  biasBN=c(biasBN, l9[[i]]$bias)}

biasBN_random=c()
for (i in 1:N.sim){  biasBN_random=c(biasBN_random, l9[[i]]$bias_random)}

bias_glm=c()
for (i in 1:N.sim){  bias_glm=c(bias_glm, l9[[i]]$bias_glm)}


mean(biasBN)
mean(biasBN_random)
mean(bias_glm)

se_BN=c()
for (i in 1:N.sim){  se_BN=c(se_BN, l9[[i]]$se_BN)}

length_BN=c()
for (i in 1:N.sim){  length_BN=c(length_BN, l9[[i]]$length_BN)}
mean(length_BN)


se_BN_random=c()
for (i in 1:N.sim){  se_BN_random=c(se_BN_random, l9[[i]]$se_BN_random)}

length_BN_random=c()
for (i in 1:N.sim){  length_BN_random=c(length_BN_random, l9[[i]]$length_BN_random)}
mean(length_BN_random)


se_glm=c()
for (i in 1:N.sim){  se_glm=c(se_glm, l9[[i]]$se_glm)}

length_glm=c()
for (i in 1:N.sim){  length_glm=c(length_glm, l9[[i]]$length_glm)}
mean(length_glm)


coverageBN=c()
for (i in 1:N.sim){  coverageBN=c(coverageBN, l9[[i]]$cov)}
mean(coverageBN)


coverageBN_random=c()
for (i in 1:N.sim){  coverageBN_random=c(coverageBN_random, l9[[i]]$cov_random)}
mean(coverageBN_random)

coverage_glm=c()
for (i in 1:N.sim){  coverage_glm=c(coverage_glm, l9[[i]]$cov_glm)}
mean(coverage_glm)

mseBN=c()
for (i in 1:N.sim){  mseBN=c(mseBN, l9[[i]]$mseBN)}
mean(mseBN)


mseBN_random=c()
for (i in 1:N.sim){  mseBN_random=c(mseBN_random, l9[[i]]$mseBN_random)}
mean(mseBN_random)


mse_glm=c()
for (i in 1:N.sim){  mse_glm=c(mse_glm, l9[[i]]$mse_glm)}
mean(mse_glm)


counterBN=c()
for (i in 1:N.sim){counterBN=c(counterBN, l9[[i]]$counterBN)}
sum(counterBN)


counterBN_random=c()
for (i in 1:N.sim){counterBN_random=c(counterBN_random, l9[[i]]$counterBN_random)}
sum(counterBN_random)


heterogeneity=c()
for (i in 1:N.sim){heterogeneity=c(heterogeneity, l9[[i]]$heterogeneity)}
tau2_BN=heterogeneity


heterogeneity_random=c()
for (i in 1:N.sim){heterogeneity_random=c(heterogeneity_random, l9[[i]]$heterogeneity_random)}
tau2_BN_random=heterogeneity_random


lower_BN=c()
for (i in 1:N.sim){lower_BN=c(lower_BN, l9[[i]]$lower_BN)}

upper_BN=c()
for (i in 1:N.sim){upper_BN=c(upper_BN, l9[[i]]$upper_BN)}

estimate_BN=c()
for (i in 1:N.sim){estimate_BN=c(estimate_BN, l9[[i]]$estimate_BN)}

lower_BN_random=c()
for (i in 1:N.sim){lower_BN_random=c(lower_BN_random, l9[[i]]$lower_BN_random)}

upper_BN_random=c()
for (i in 1:N.sim){upper_BN_random=c(upper_BN_random, l9[[i]]$upper_BN_random)}

estimate_BN_random=c()
for (i in 1:N.sim){estimate_BN_random=c(estimate_BN_random, l9[[i]]$estimate_BN_random)}

lower_glm=c()
for (i in 1:N.sim){lower_glm=c(lower_glm, l9[[i]]$lower_glm)}

upper_glm=c()
for (i in 1:N.sim){upper_glm=c(upper_glm, l9[[i]]$upper_glm)}

estimate_glm=c()
for (i in 1:N.sim){estimate_glm=c(estimate_glm, l9[[i]]$estimate_glm)}

count_fail_fixed=c()
for (i in 1:N.sim){count_fail_fixed=c(count_fail_fixed, l9[[i]]$count_fail_fixed)}


count_fail_random=c()
for (i in 1:N.sim){count_fail_random=c(count_fail_random, l9[[i]]$count_fail_random)}

count_fail=cbind.data.frame(count_fail_fixed,count_fail_random)

summaries_fail=data.frame(sum(count_fail[,1]),sum(count_fail[,2]))
names(summaries_fail)=c("fail fixed intercept","fail random intercept")

new_het=cbind.data.frame(heterogeneity,heterogeneity_random)

bias_taulme4_fixed=new_het[,1]-tau
bias_taulme4_random=new_het[,2]-tau

bias_taulme4=cbind.data.frame(bias_taulme4_fixed,bias_taulme4_random)
names(bias_taulme4)=c("fixed intercept","random intercept")

bias_taulme4_mean=cbind.data.frame(mean(bias_taulme4[,1]),mean(bias_taulme4[,2]))
names(bias_taulme4_mean)=c("fixed intercept","random intercept")

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
## BN random


round(mean(biasBN_random), digits=2)
round(mean(abs(biasBN_random)), digits=2)
round(mean(coverageBN_random), digits=3)*100
round(mean(mseBN_random),digits = 2)





############

bias=c(round(mean(biasIV.FE), digits=2),round(mean(biasIV.RE), digits=2),round(mean(biasMH), digits=2),round(mean(biasNCH), digits=2),
       round(mean(biasFL_F), digits=2),round(mean(biasFL_F), digits=2),round(mean(biasBN), digits=2),
       round(mean(biasBN_random), digits=2),
       round(mean(bias_glm), digits=2)
       
)





absolute_bias=c(round(mean(abs(biasIV.FE)), digits=2),round(mean(abs(biasIV.RE)), digits=2),round(mean(abs(biasMH)), digits=2),
                round(mean(abs(biasNCH)), digits=2),round(mean(abs(biasFL_F)), digits=2),round(mean(abs(biasFL_F)), digits=2)
                ,round(mean(abs(biasBN)), digits=2),
                round(mean(abs(biasBN_random)), digits=2),
                round(mean(abs(bias_glm)), digits=2)
                
)

coverage=c(round(mean(coverageIV.FE), digits=3)*100,round(mean(coverageIV.RE), digits=3)*100,round(mean(coverageMH), digits=3)*100,
           round(mean(coverageNCH), digits=3)*100,round(mean(coverageFL_F), digits=3)*100,round(mean(coverageFL_Fr), digits=3)*100,
           round(mean(coverageBN), digits=3)*100,
           round(mean(coverageBN_random), digits=3)*100,
           round(mean(coverage_glm), digits=3)*100
)

coverage1=c(round(mean(coverageIV.FE), digits=3)*100,round(mean(coverageIV.RE), digits=3)*100,round(mean(coverageMH), digits=3)*100,
            round(mean(coverageNCH), digits=3)*100,round(mean(pcoverageFL_F), digits=3)*100,round(mean(coverageFL_Fr), digits=3)*100,
            round(mean(coverageBN), digits=3)*100,
            
            round(mean(coverageBN_random), digits=3)*100,
            round(mean(coverage_glm), digits=3)*100
)

mse=c(round(mean(mseIV.FE),digits = 2),round(mean(mseIV.RE),digits = 2),round(mean(mseMH),digits = 2),round(mean(mseNCH),digits = 2),
      round(mean(mseFL_F),digits = 2),round(mean(mseFL_Fr),digits = 2),round(mean(mseBN),digits = 2),
      round(mean(mseBN_random),digits = 2),round(mean(mse_glm),digits = 2)
      
)

results=cbind(bias,absolute_bias,mse,coverage)
results=data.frame(results)
row.names(results)=c("IV_FE","IV_RE","MH","NCH","FL_F","FL_F(R)","B-N","B-N random intercept","Common logistic")
results

results1=cbind(bias,absolute_bias,mse,coverage1)
results1=data.frame(results1)
row.names(results1)=c("IV_FE","IV_RE","MH","NCH","FL_F(P)","FL_F(R)","B-N","B-N random intercept","Common logistic")
results1


het1=cbind(counter,dispersion)
het1=as.data.frame(het1)
names(het1)=c("counter_Fletcher","dispersion")
het2=cbind(counterBN,tau2_BN,counterBN_random,tau2_BN_random)
het2=as.data.frame(het2)

het=data.frame(het1,het2)


tau_info_IV=data.frame(bias_tau_IV,tau_IV)
names(tau_info_IV)=c("bias_tau_IV","tau_IV")

countIV=ifelse(tau_IV>0,1,0)


