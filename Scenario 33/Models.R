library(parallel)
library(beepr)
library(netmeta)
library(brglm)
library(lme4)



####################### PLNMA ##########################################################
no_cores1 <- detectCores()-2
cl1 <- makeCluster(no_cores1)

X2=list()

for (i in 1:N.sim){ X2[[i]]=list("data"=data2[[i]],"logOR"=logOR[[i]]) } ### runs it in the dataset with all-zero events studies

#for (i in 1:N.sim){ X2[[i]]=list("data"=data20[[i]],"logOR"=logOR[[i]]) } ### runs it in the dataset without all-zero events studies

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

counter_PLNMA=c()
for (i in 1:N.sim){counter_PLNMA=c(counter_PLNMA, l8[[i]]$count)}
sum(counter_PLNMA)
#beep(sound=3)

dispersion_PLNMA=c()
for (i in 1:N.sim){dispersion_PLNMA=c(dispersion_PLNMA, l8[[i]]$dispersion)}
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


######################################################################################################


no_cores2 <- detectCores()
cl2 <- makeCluster(no_cores2)


X3=list()

for (i in 1:N.sim){ X3[[i]]=list("data"=data2[[i]],"logOR"=logOR[[i]]) } ### runs it in the dataset with all-zero events studies

#for (i in 1:N.sim){ X2[[i]]=list("data"=data20[[i]],"logOR"=logOR[[i]]) } ### runs it in the dataset without all-zero events studies

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
  
  # heterogeneity_random=data.frame(VarCorr(BN2))[,4]
  # counterBN_random=(heterogeneity>0)
  
  lower_BN_random=c(lower_BN_random,BN.res_random$lowerCI)
  upper_BN_random=c(upper_BN_random,BN.res_random$upperCI)
  estimate_BN_random=c(estimate_BN_random,BN.res_random$mean)
  
  heterogeneity_random=VarCorr(BN2)
  C=as.data.frame(heterogeneity_random)
  
  C=subset(C,C$var1=="treat")
  heterogeneity_random=C$sdcor
  counterBN_random=(heterogeneity_random>0)
  
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

########### Bayesian models ################################

no_cores1 <- detectCores()-1
cl11 <- makeCluster(no_cores1)

X4=list()

for (i in 1:N.sim){ X4[[i]]=list("data"=data3[[i]],"logOR"=logOR[[i]],"tau"=tau[[1]])} ### runs in the dataset with all zero events studies
#for (i in 1:N.sim){ X4[[i]]=list("data"=data30[[i]],"logOR"=logOR[[i]],"tau"=tau[[1]])} ### runs in the dataset without all zero events studies

BL_NMA=function(X)
{
  
  biasBL_NMA=c()
  coverageBL_NMA=c()
  mseBL_NMA=c()
  se_BL_NMA=c()
  length_BL_NMA=c()
  lower_BL_NMA=c()
  upper_BL_NMA=c()
  estimate_BL_NMA=c()
  # conv_ests=c()
  # conv_tau=c()
  
  
  biasBL_NMA_Random=c()
  biasBL_NMA_Random_2=c()
  bias_tau=c()
  bias_tau_2=c()
  coverageBL_NMA_Random=c()
  coverageBL_NMA_Random_2=c()
  mseBL_NMA_Random=c()
  mseBL_NMA_Random_2=c()
  se_BL_NMA_Random=c()
  se_BL_NMA_Random_2=c()
  length_BL_NMA_Random=c()
  length_BL_NMA_Random_2=c()
  lower_BL_NMA_Random=c()
  lower_BL_NMA_Random_2=c()
  upper_BL_NMA_Random=c()
  upper_BL_NMA_Random_2=c()
  estimate_BL_NMA_Random=c()
  estimate_BL_NMA_Random_2=c()
  
  
  
  lower_tau=c()
  lower_tau_2=c()
  upper_tau=c()
  upper_tau_2=c()
  length_tau=c()
  length_tau_2=c()
  se_tau=c()
  se_tau_2=c()
  
  conv_ests=c()
  conv_ests_random=c()
  conv_ests_random_informative=c()
  name=c()
  
  conv_tau=c()
  conv_tau_informative=c()
  
  mtc.network <- mtc.network(data.ab =X$data, description = "Network")
  
  # specify the estimation parameters 
  mtc.model1 <-mtc.model(mtc.network, type ="consistency", om.scale = 2,
                         hy.prior=mtc.hy.prior("std.dev","dunif", 0, 2),linearModel = "fixed",re.prior.sd = 100,n.chain = 2)
  
  mtc.model1$inits[[1]]$.RNG.name <- "base::Mersenne-Twister"
  mtc.model1$inits[[2]]$.RNG.name <- "base::Mersenne-Twister"
  #mtc.model1$inits[[3]]$.RNG.name <- "base::Mersenne-Twister"
  #mtc.model1$inits[[4]]$.RNG.name <- "base::Mersenne-Twister"
  
  mtc.model1$inits[[1]]$.RNG.seed <- 58982
  mtc.model1$inits[[2]]$.RNG.seed <- 58983
  #mtc.model1$inits[[3]]$.RNG.seed <- 58984
  #mtc.model1$inits[[4]]$.RNG.seed <- 58985
  mtc.run1 <- mtc.run(mtc.model1, thin = 1,n.iter = 50000,n.adapt = 10000)
  
  conv_fixed=gelman.diag(mtc.run1)
  conv_fixed=conv_fixed$psrf[,1]
  l_fixed=length(conv_fixed)
  conv_ests=unname(conv_fixed[1:l_fixed])
  name=names(conv_fixed)
  # conv_ests1=cbind.data.frame(conv_ests,name)
  #conv_tau=conv_fixed[l_fixed]
  
  s1=summary(mtc.run1)$summaries[1]
  s2=summary(mtc.run1)$summaries[2]
  s1=as.data.frame(s1)
  s2=as.data.frame(s2)
  res=cbind(s1[,1],s1[,2],s2[,1],s2[,5])
  
  BL_res=as.data.frame(res)
  names(BL_res)=c("mean","se","lowerCR","upperCR")
  
  biasBL_NMA=c(biasBL_NMA,(BL_res$mean-X$logOR))
  BL_res$cover=(BL_res$lowerCR<X$logOR)&(BL_res$upperCR>X$logOR)
  
  mseBL_NMA=c(mseBL_NMA,biasBL_NMA*biasBL_NMA+BL_res[,2]*BL_res[,2])
  
  coverageBL_NMA=c(coverageBL_NMA, BL_res$cover)
  
  se_BL_NMA=c(se_BL_NMA,BL_res$se)
  lower_BL_NMA=c(lower_BL_NMA,BL_res$lowerCR)
  upper_BL_NMA=c(upper_BL_NMA,BL_res$upperCR)
  length_BL_NMA=c(length_BL_NMA,BL_res$upperCR-BL_res$lowerCR)
  estimate_BL_NMA=c(estimate_BL_NMA,BL_res$mean)
  ########################################################################################
  ########################################################################################
  
  mtc.model11 <-mtc.model(mtc.network, type ="consistency", om.scale = 2,
                          hy.prior=mtc.hy.prior("std.dev","dunif", 0, 2),linearModel = "random",re.prior.sd = 100,n.chain = 2)
  
  hy.prior <- mtc.hy.prior(type="std.dev", distr="dhnorm", 0, 9.77)
  
  mtc.model12 <-mtc.model(mtc.network, type ="consistency",
                          hy.prior=hy.prior,linearModel = "random",re.prior.sd = 10,n.chain = 2)
  
  
  mtc.model11$inits[[1]]$.RNG.name <- "base::Mersenne-Twister"
  mtc.model11$inits[[2]]$.RNG.name <- "base::Mersenne-Twister"
  
  mtc.model12$inits[[1]]$.RNG.name <- "base::Mersenne-Twister"
  mtc.model12$inits[[2]]$.RNG.name <- "base::Mersenne-Twister"
  
  mtc.model11$inits[[1]]$.RNG.seed <- 58962
  mtc.model11$inits[[2]]$.RNG.seed <- 58963
  
  mtc.model12$inits[[1]]$.RNG.seed <- 58964
  mtc.model12$inits[[2]]$.RNG.seed <- 58965
  
  mtc.run11 <- mtc.run(mtc.model11, thin = 1,n.iter = 50000,n.adapt = 10000)
  
  conv_random=gelman.diag(mtc.run11)
  conv_random=conv_random$psrf[,1]
  l_random=length(conv_random)
  conv_ests_random=unname(conv_random[1:l_random-1])
  # name_random=names(conv_ests_random)
  # conv_ests_random_1=cbind.data.frame(conv_ests_random,name_random)
  conv_tau=conv_random[l_random]
  
  
  mtc.run12 <- mtc.run(mtc.model12, thin = 1,n.iter = 50000,n.adapt = 10000)
  
  conv_random_informative=gelman.diag(mtc.run12)
  conv_random_informative=conv_random_informative$psrf[,1]
  l_random_informative=length(conv_random_informative)
  conv_ests_random_informative=unname(conv_random_informative[1:l_random_informative-1])
  # name_random_informative=names(conv_ests_random_informative)
  # conv_ests_random_informative_1=cbind.data.frame(conv_ests_random_informative,name_random_informative)
  conv_tau_informative=conv_random_informative[l_random_informative]
  
  #convergence_estimates=rbind.data.frame(conv_ests,conv_ests_random,conv_ests_random_informative,name)
  
  s11=summary(mtc.run11)$summaries[1]
  s21=summary(mtc.run11)$summaries[2]
  s11=as.data.frame(s11)
  s21=as.data.frame(s21)
  res1=cbind(s11[,1],s11[,2],s21[,1],s21[,5])
  
  
  s11_2=summary(mtc.run12)$summaries[1]
  s21_2=summary(mtc.run12)$summaries[2]
  s11_2=as.data.frame(s11_2)
  s21_2=as.data.frame(s21_2)
  res12=cbind(s11_2[,1],s11_2[,2],s21_2[,1],s21_2[,5])
  
  
  tau1=res1[length(res1[,1]),1]
  tau12=res12[length(res12[,1]),1]
  
  lower_tau=c(lower_tau,res1[length(res1[,1]),3])
  upper_tau=c(upper_tau,res1[length(res1[,1]),4])
  se_tau=c(se_tau,res1[length(res1[,1]),2])
  length_tau=c(length_tau,res1[length(res1[,1]),4]-res1[length(res1[,1]),3])
  
  
  lower_tau_2=c(lower_tau_2,res12[length(res12[,1]),3])
  upper_tau_2=c(upper_tau_2,res12[length(res12[,1]),4])
  se_tau_2=c(se_tau_2,res1[length(res12[,1]),2])
  length_tau_2=c(length_tau_2,res12[length(res12[,1]),4]-res12[length(res12[,1]),3])
  
  tau_test=res1[length(res1[,1]),1]
  tau_test_2=res12[length(res12[,1]),1]
  
  if(tau1<0 | tau1==0){
    tau1=0}
  if(tau1>0){
    tau1=tau1
  }
  res1=res1[-c(length(res1[,1])),]
  
  if(tau12<0 | tau12==0){
    tau12=0}
  if(tau12>0){
    tau12=tau12
  }
  res12=res12[-c(length(res12[,1])),]
  
  
  BL_res1=as.data.frame(res1)
  names(BL_res1)=c("mean","se","lowerCR","upperCR")
  
  BL_res12=as.data.frame(res12)
  names(BL_res12)=c("mean","se","lowerCR","upperCR")
  
  biasBL_NMA_Random=c(biasBL_NMA_Random,(BL_res1$mean-X$logOR))
  biasBL_NMA_Random_2=c(biasBL_NMA_Random_2,(BL_res12$mean-X$logOR))
  
  bias_tau=c(bias_tau,(tau1-X$tau))
  bias_tau_2=c(bias_tau_2,(tau12-X$tau))
  
  BL_res1$cover=(BL_res1$lowerCR<X$logOR)&(BL_res1$upperCR>X$logOR)
  BL_res12$cover=(BL_res12$lowerCR<X$logOR)&(BL_res12$upperCR>X$logOR)
  
  mseBL_NMA_Random=c(mseBL_NMA_Random,biasBL_NMA_Random*biasBL_NMA_Random+BL_res1[,2]*BL_res1[,2])
  mseBL_NMA_Random_2=c(mseBL_NMA_Random_2,biasBL_NMA_Random_2*biasBL_NMA_Random_2+BL_res12[,2]*BL_res12[,2])
  
  coverageBL_NMA_Random=c(coverageBL_NMA_Random, BL_res1$cover)
  coverageBL_NMA_Random_2=c(coverageBL_NMA_Random_2, BL_res12$cover)
  
  estimate_BL_NMA_Random=c(estimate_BL_NMA_Random,BL_res1$mean)
  estimate_BL_NMA_Random_2=c(estimate_BL_NMA_Random_2,BL_res12$mean)
  
  se_BL_NMA_Random=c(se_BL_NMA_Random,BL_res1$se)
  se_BL_NMA_Random_2=c(se_BL_NMA_Random_2,BL_res12$se)
  
  lower_BL_NMA_Random=c(lower_BL_NMA_Random,BL_res1$lowerCR)
  lower_BL_NMA_Random_2=c(lower_BL_NMA_Random_2,BL_res12$lowerCR)
  
  upper_BL_NMA_Random=c(upper_BL_NMA_Random,BL_res1$upperCR)
  upper_BL_NMA_Random_2=c(upper_BL_NMA_Random_2,BL_res12$upperCR)
  
  
  length_BL_NMA_Random=c(length_BL_NMA_Random,BL_res1$upperCR-BL_res1$lowerCR)
  length_BL_NMA_Random_2=c(length_BL_NMA_Random_2,BL_res12$upperCR-BL_res12$lowerCR)
  
  
  
  return(list("bias_fixed"=biasBL_NMA,"cov_fixed"=coverageBL_NMA,"mse_fixed"=mseBL_NMA,
              "estimate_fixed_bayes"=estimate_BL_NMA,"se_fixed_bayes"=se_BL_NMA,"lower_fixed_bayes"=lower_BL_NMA,"upper_fixed_bayes"=upper_BL_NMA,
              "length_fixed_bayes"=length_BL_NMA,
              "bias_random"=biasBL_NMA_Random,"cov_random"=coverageBL_NMA_Random,"mse_random"=mseBL_NMA_Random,"bias_tau"=bias_tau,
              "tau"=tau1,"tau_test"=tau_test,
              "estimate_random_bayes"=estimate_BL_NMA_Random,"se_random_bayes"=se_BL_NMA_Random,"lower_random_bayes"=lower_BL_NMA_Random,
              "upper_random_bayes"=upper_BL_NMA_Random,"length_random_bayes"=length_BL_NMA_Random,
              "lower_tau"=lower_tau,"upper_tau"=upper_tau,"se_tau"=se_tau,"length_tau"=length_tau,
              
              "bias_random_2"=biasBL_NMA_Random_2,"cov_random_2"=coverageBL_NMA_Random_2,"mse_random_2"=mseBL_NMA_Random_2,"bias_tau_2"=bias_tau_2,
              "tau_2"=tau12,"tau_test_2"=tau_test_2,
              "estimate_random_bayes_2"=estimate_BL_NMA_Random_2,"se_random_bayes_2"=se_BL_NMA_Random_2,"lower_random_bayes_2"=lower_BL_NMA_Random_2,
              "upper_random_bayes_2"=upper_BL_NMA_Random_2,"length_random_bayes_2"=length_BL_NMA_Random_2,
              "lower_tau_2"=lower_tau_2,"upper_tau_2"=upper_tau_2,"se_tau_2"=se_tau_2,"length_tau_2"=length_tau_2,
              
              "conv_ests_fixed"=conv_ests,"conv_ests_random"=conv_ests_random,"conv_ests_random_informative"=conv_ests_random_informative,
              "names"=name,"conv_tau"=conv_tau,"conv_tau_informative"=conv_tau_informative
  ))
}




clusterExport(cl11,"X4")
clusterExport(cl11,"NT")
clusterExport(cl11,"N.sim")
clusterExport(cl11, "BL_NMA")
clusterEvalQ(cl11, {library(gemtc)})
l10=parLapply(cl11,1:N.sim, function(x) BL_NMA(X4[[x]]))

convergence_fixed=c()
convergence_random=c()
convergence_random_informative=c()
convergence_tau_random=c()
convergence_tau_random_informative=c()
name=c()
for (i in 1:N.sim){  convergence_fixed=c(convergence_fixed, l10[[i]]$conv_ests_fixed)}
for (i in 1:N.sim){  convergence_random=c(convergence_random, l10[[i]]$conv_ests_random)}
for (i in 1:N.sim){  convergence_random_informative=c(convergence_random_informative, l10[[i]]$conv_ests_random_informative)}
for (i in 1:N.sim){  name=c(name, l10[[i]]$names)}
for (i in 1:N.sim){  convergence_tau_random=c(convergence_tau_random, l10[[i]]$conv_tau)}
for (i in 1:N.sim){  convergence_tau_random_informative=c(convergence_tau_random_informative, l10[[i]]$conv_tau_informative)}

convergence_summary=cbind.data.frame(convergence_fixed,convergence_random,convergence_random_informative,name)


convergence_tau=cbind.data.frame(convergence_tau_random,convergence_tau_random_informative)

names(convergence_summary)=c("Common effect","Random-effects (d~N(0,10000),tau~U(0,2))","Random-effects (d~N(0,100),tau~HN(1))","parameter")
names(convergence_tau)=c("Uniform(0,2)","Half-Normal")

results_convergence=c("Convergence times-Common effect"=length(which(convergence_summary[,1]<1.1)),
                      "Convergence times-Random effects (d~N(0,10000),tau~U(0,2))"=length(which(convergence_summary[,2]<1.1)),
                      "Convergence times-Random effects (d~N(0,100),tau~HN(1))"=length(which(convergence_summary[,3]<1.1)),
                      "Convergence times-tau (d~N(0,10000),tau~U(0,2))"=length(which(convergence_tau[,1]<1.1)),
                      "Convergence times-tau (d~N(0,100),tau~HN(1))"=length(which(convergence_tau[,2]<1.1))
)

results_convergence1=matrix(results_convergence,nrow=1)
results_convergence1=as.data.frame(results_convergence1)
names(results_convergence1)=names(results_convergence)




biasBL_NMA_fixed=c()
for (i in 1:N.sim){  biasBL_NMA_fixed=c(biasBL_NMA_fixed, l10[[i]]$bias_fixed)}
mean(biasBL_NMA_fixed)


coverageBL_NMA_fixed=c()
for (i in 1:N.sim){  coverageBL_NMA_fixed=c(coverageBL_NMA_fixed, l10[[i]]$cov_fixed)}
mean(coverageBL_NMA_fixed)

mseBL_NMA_fixed=c()
for (i in 1:N.sim){  mseBL_NMA_fixed=c(mseBL_NMA_fixed, l10[[i]]$mse_fixed)}
mean(mseBL_NMA_fixed)


biasBL_NMA_random=c()
for (i in 1:N.sim){  biasBL_NMA_random=c(biasBL_NMA_random, l10[[i]]$bias_random)}
mean(biasBL_NMA_random)

biasBL_NMA_random_2=c()
for (i in 1:N.sim){  biasBL_NMA_random_2=c(biasBL_NMA_random_2, l10[[i]]$bias_random_2)}
mean(biasBL_NMA_random_2)

bias_tau=c()
for (i in 1:N.sim){  bias_tau=c(bias_tau, l10[[i]]$bias_tau)}
mean(bias_tau)

bias_tau_2=c()
for (i in 1:N.sim){  bias_tau_2=c(bias_tau_2, l10[[i]]$bias_tau_2)}
mean(bias_tau_2)

tau=c()
for (i in 1:N.sim){tau=c(tau, l10[[i]]$tau)}

tau_2=c()
for (i in 1:N.sim){tau_2=c(tau_2, l10[[i]]$tau_2)}

se_tau=c()
for (i in 1:N.sim){se_tau=c(se_tau, l10[[i]]$se_tau)}

se_tau_2=c()
for (i in 1:N.sim){se_tau_2=c(se_tau_2, l10[[i]]$se_tau_2)}

upper_tau=c()
for (i in 1:N.sim){upper_tau=c(upper_tau, l10[[i]]$upper_tau)}

upper_tau_2=c()
for (i in 1:N.sim){upper_tau_2=c(upper_tau_2, l10[[i]]$upper_tau_2)}

lower_tau=c()
for (i in 1:N.sim){lower_tau=c(lower_tau, l10[[i]]$lower_tau)}

lower_tau_2=c()
for (i in 1:N.sim){lower_tau_2=c(lower_tau_2, l10[[i]]$lower_tau_2)}

length_tau=c()
for (i in 1:N.sim){length_tau=c(length_tau, l10[[i]]$length_tau)}
mean(length_tau)

length_tau_2=c()
for (i in 1:N.sim){length_tau_2=c(length_tau_2, l10[[i]]$length_tau_2)}
mean(length_tau_2)

tau_test=c()
for (i in 1:N.sim){tau_test=c(tau_test, l10[[i]]$tau_test)}

tau_test_2=c()
for (i in 1:N.sim){tau_test_2=c(tau_test_2, l10[[i]]$tau_test_2)}



coverageBL_NMA_random=c()
for (i in 1:N.sim){  coverageBL_NMA_random=c(coverageBL_NMA_random, l10[[i]]$cov_random)}
mean(coverageBL_NMA_random)

coverageBL_NMA_random_2=c()
for (i in 1:N.sim){  coverageBL_NMA_random_2=c(coverageBL_NMA_random_2, l10[[i]]$cov_random_2)}
mean(coverageBL_NMA_random_2)


mseBL_NMA_random=c()
for (i in 1:N.sim){  mseBL_NMA_random=c(mseBL_NMA_random, l10[[i]]$mse_random)}
mean(mseBL_NMA_random)

mseBL_NMA_random_2=c()
for (i in 1:N.sim){  mseBL_NMA_random_2=c(mseBL_NMA_random_2, l10[[i]]$mse_random_2)}
mean(mseBL_NMA_random)


se_fixed_bayes=c()
for (i in 1:N.sim){se_fixed_bayes=c(se_fixed_bayes, l10[[i]]$se_fixed_bayes)}

lower_fixed_bayes=c()
for (i in 1:N.sim){lower_fixed_bayes=c(lower_fixed_bayes, l10[[i]]$lower_fixed_bayes)}

upper_fixed_bayes=c()
for (i in 1:N.sim){upper_fixed_bayes=c(upper_fixed_bayes, l10[[i]]$upper_fixed_bayes)}

length_fixed_bayes=c()
for (i in 1:N.sim){length_fixed_bayes=c(length_fixed_bayes, l10[[i]]$length_fixed_bayes)}
mean(length_fixed_bayes)

estimate_fixed_bayes=c()
for (i in 1:N.sim){estimate_fixed_bayes=c(estimate_fixed_bayes, l10[[i]]$estimate_fixed_bayes)}




se_random_bayes=c()
for (i in 1:N.sim){se_random_bayes=c(se_random_bayes, l10[[i]]$se_random_bayes)}

se_random_bayes_2=c()
for (i in 1:N.sim){se_random_bayes_2=c(se_random_bayes_2, l10[[i]]$se_random_bayes_2)}

lower_random_bayes=c()
for (i in 1:N.sim){lower_random_bayes=c(lower_random_bayes, l10[[i]]$lower_random_bayes)}

lower_random_bayes_2=c()
for (i in 1:N.sim){lower_random_bayes_2=c(lower_random_bayes_2, l10[[i]]$lower_random_bayes_2)}

upper_random_bayes=c()
for (i in 1:N.sim){upper_random_bayes=c(upper_random_bayes, l10[[i]]$upper_random_bayes)}

upper_random_bayes_2=c()
for (i in 1:N.sim){upper_random_bayes_2=c(upper_random_bayes_2, l10[[i]]$upper_random_bayes_2)}

length_random_bayes=c()
for (i in 1:N.sim){length_random_bayes=c(length_random_bayes, l10[[i]]$length_random_bayes)}
mean(length_random_bayes)

length_random_bayes_2=c()
for (i in 1:N.sim){length_random_bayes_2=c(length_random_bayes_2, l10[[i]]$length_random_bayes_2)}
mean(length_random_bayes_2)

estimate_random_bayes=c()
for (i in 1:N.sim){estimate_random_bayes=c(estimate_random_bayes, l10[[i]]$estimate_random_bayes)}

estimate_random_bayes_2=c()
for (i in 1:N.sim){estimate_random_bayes_2=c(estimate_random_bayes_2, l10[[i]]$estimate_random_bayes_2)}



####################### Results #####################################

rows=c("PLNMA","PLNMA profile","PLNMA random","Binomial-Normal (fixed-intercept)","Binomial-Normal (random interecepts)",
       "Logistic Common-Standard Likelihood","Bayesian logistic-Common","Bayesian Random ((d~N(0,10000),tau~U(0,2)))","Bayesian Random ((d~N(0,100),tau~HN(1)))")


results_bias=cbind.data.frame("PLNMA"=biasFL_F,"PLNMA profile"=biasFL_F,"PLNMA random"=biasFL_F,
                              biasBN,biasBN_random,bias_glm,biasBL_NMA_fixed,biasBL_NMA_random,biasBL_NMA_random_2)

colnames(results_bias)=rows

mean_bias=as.data.frame(sapply(results_bias,mean))


mean_absolute_bias=as.data.frame(sapply(abs(results_bias),mean))

results_coverage=cbind.data.frame(coverageFL_F,pcoverageFL_F,coverageFL_Fr,
                                  coverageBN,coverageBN_random,coverage_glm,
                                  coverageBL_NMA_fixed,coverageBL_NMA_random,coverageBL_NMA_random_2)

colnames(results_coverage)=rows


mean_coverage=as.data.frame(sapply(results_coverage,mean))
mean_coverage=mean_coverage*100



results_mse=cbind.data.frame(mseFL_F,mseFL_F,mseFL_Fr,mseBN,mseBN_random,mse_glm,mseBL_NMA_fixed,mseBL_NMA_random,mseBL_NMA_random_2)
colnames(results_mse)=rows
mean_mse=as.data.frame(sapply(results_mse,mean))


results_length=cbind.data.frame(length_fixed,length_profile,length_random,
                                length_BN,length_BN_random,length_glm,
                                length_fixed_bayes,length_random_bayes,length_random_bayes_2)


colnames(results_length)=rows
mean_length=as.data.frame(sapply(results_length,mean))

results=cbind.data.frame(mean_bias,mean_absolute_bias,mean_coverage,mean_mse,mean_length)
results$model=rownames(results)
results=results[,c(6,1,2,3,4,5)]
rownames(results)=NULL
names(results)=c("model","bias","absolute bias","coverage","mse","length")



estimates=cbind.data.frame(estimatef,estimatef,estimatef,estimate_BN,estimate_BN_random,estimate_glm,estimate_fixed_bayes,estimate_random_bayes,estimate_random_bayes_2)
names(estimates)=rows

se=cbind.data.frame(se_fixed,se_fixed,se_random,se_BN,se_BN_random,se_glm,se_fixed_bayes,se_random_bayes,se_random_bayes_2)
names(se)=rows

bias_tau=cbind.data.frame(bias_tau,bias_tau_2)
names(bias_tau)=c("Uniform","Half_normal")
tau_info=data.frame(bias_tau$Uniform,tau,tau_test,bias_tau$Half_normal,tau_2,tau_test_2)
names(tau_info)=c("bias_tau","tau","tau_test","bias_tau_informative","tau_informative","tau_test_informative")
interval=cbind.data.frame(lower_fixed,upper_fixed,
                          lower_profile,upper_profile,
                          lower_random,upper_random,
                          lower_BN,upper_BN,
                          lower_BN_random,upper_BN_random,
                          lower_glm,upper_glm,
                          lower_fixed_bayes,upper_fixed_bayes,
                          lower_random_bayes,upper_random_bayes,
                          lower_random_bayes_2,upper_random_bayes_2)
                          
het1=cbind(counter_PLNMA,dispersion_PLNMA,counterBN,tau2_BN,counterBN_random,tau2_BN_random)
het1=as.data.frame(het1)
names(het1)=c("counter_Fletcher","dispersion","counter BN","tau BN","counter BN random intercepts","tau random intercept")                      



#####################################################################################################################################
write.csv(se,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\standard errors.csv",row.names = F)
write.csv(estimates,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\estimates.csv",row.names = F)
write.csv(results,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\main_results.csv",row.names = F)
write.csv(results_bias,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\bias.csv",row.names = F)
write.csv(results_coverage,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\coverage.csv",row.names = F)
write.csv(results_mse,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\mse.csv",row.names = F)
write.csv(results_length,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\length.csv",row.names = F)

##############################################################
###################### heterogeneity-convergence for binomial normal ##########################

write.csv(convergence_summary,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\convergence_estimates_bayes.csv",row.names = F)
write.csv(convergence_tau,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\convergence_tau_bayes.csv",row.names = F)
write.csv(results_convergence1,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\overall_results_convergence_bayes.csv",row.names =F)

##############################################################
######## heterogeneity-convergence for binomial normal ##########################
write.csv(count_fail,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\convergence fail.csv")
write.csv(summaries_fail,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\convergence fail_summaries.csv")
write.csv(new_het,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\heterogeneity.csv")
write.csv(bias_taulme4,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\bias_taulme4.csv")
write.csv(bias_taulme4_mean,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\mean_bias_taulme4.csv")
#######################################
write.csv(tau_info,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\tau_info_bayes.csv",row.names =F)
write.csv(interval,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\intervals.csv",row.names =F)
write.csv(het1,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Paper\\Simulation NMA\\Scenaria 1-16\\More zeroes\\t=0\\Included\\heterogeneity_info.csv",row.names =F)
