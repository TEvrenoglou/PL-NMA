library(tidyverse)
library(brglm)
library(lme4)
library(netmeta)
library(AICcmodavg)
library(gemtc)

data("Dong2013") #### The dataset is already embedded in netmeta package


PLNMA=function(data,events,n,study,treat,ref.group){
  require(brglm)
  require(AICcmodavg)
  
  treat1=as.factor(data$treat)
  treat1=relevel(treat1,ref = ref.group)
  
  dispersion_Fletcher=function(model){
    phi1=c_hat(model,method = "fletcher")
    if(phi1>1){
      phi=phi1}
    else if((phi1<1)||(phi1=1)){
      phi=1}
    return(phi)}
  
  PL=brglm(cbind(events,n-events)~factor(treat1)+factor(study),data=data)
  phi=dispersion_Fletcher(PL)
  ests=summary(PL)$coefficients
  ests=ests[grep("treat",rownames(ests)),]
  ests=ests[,-c(3,4)]
  ci.lower=ests[,1]-1.96*ests[,2]
  ci.upper=ests[,1]+1.96*ests[,2]
  
  ci.rlower=ests[,1]-1.96*(ests[,2]*sqrt(phi))
  ci.rupper=ests[,1]+1.96*(ests[,2]*sqrt(phi))
  
  n.se=ests[,2]*sqrt(phi)
  
  ests=cbind(ests,ci.lower,ci.upper)
  ests=as.data.frame(ests)
  
  
  
  PL_res=data.frame(OR=exp(ests[,1]),lowerCI=exp(ci.lower),upperCI=exp(ci.upper),
                    random_lowerCI=exp(ci.rlower),random_upperCI=exp(ci.rupper))
  PL_res=round(PL_res,digits = 2)
  PL_res=as.matrix(PL_res)
  
  row.names(PL_res)=paste(levels(treat1)[-1],"vs",ref.group)
  
  PL_res=list(main_results=PL_res,overdispersion=phi)
  return(PL_res)
}

p1 = pairwise(treatment, death, randomized, studlab = id,
              data = Dong2013, sm = "OR")


################## MH-NMA model #####################

MH1 = netmetabin(p1, ref = "Placebo")

#############################################

#################### NCH-NMA model #####################

NCH1 = netmetabin(p1, ref = "Placebo", method = "NCH")

##############################################


modelPLNMA=PLNMA(data = Dong2013,events = Dong2013$death,n=Dong2013$randomized,treat = Dong2013$treatment,ref.group = "Placebo",study = Dong2013$id)

print(modelPLNMA)

#### Calculating of profile likelihood confidence intervals

Dong2013$treatment=as.factor(Dong2013$treatment)
Dong2013$treatment=relevel(Dong2013$treatment,ref = "Placebo")


model_profile=brglm(cbind(death,randomized-death)~factor(treatment)+factor(id),data=Dong2013)

### Calculating profile likelihood confidence intervals for OR
A=profileModel(model_profile,objective = "modifiedScoreStatistic",X = model.matrix(model_profile)[,!is.na(model_profile$coefficients)],
               quantile=qchisq(0.95, 1))
A1=profConfint(A)
prof_lik_conf1=round(exp(A1),digits = 2)
round(prof_lik_conf1,digits = 2)



############ Common-effect NMA #######################
model_glm=glm(cbind(death,randomized-death)~factor(treatment)+factor(id),data=Dong2013,
                family = binomial)


ests_glm=summary(model_glm)$coefficients
ests_glm=ests_glm[grep("treat",rownames(ests_glm)),]
OR_glm=exp(ests_glm[,1])
ests_glm=cbind(ests_glm,OR_glm)

ci.lower_glm=ests_glm[,1]-1.96*ests_glm[,2]
ci.upper_glm=ests_glm[,1]+1.96*ests_glm[,2]

ci.lower_glm=exp(ci.lower_glm)
ci.upper_glm=exp(ci.upper_glm)

ests_glm=cbind(ests_glm,ci.lower_glm,ci.upper_glm)
ests_glm=round(ests_glm,digits = 2)
print(ests_glm)

########################### Preparing data for Bayesian models #########################

treat=matrix(ncol=1,nrow=length(Dong2013$treatment))


for(i in 1:length(Dong2013$treatment)){
  
  
  if(Dong2013$treatment[i]=="Placebo"){
    treat[i,1]=1}
  else if (Dong2013$treatment[i]=="ICS"){
    treat[i,1]=2}
  else if (Dong2013$treatment[i]=="LABA"){
    treat[i,1]=3}
  else if (Dong2013$treatment[i]=="LABA-ICS"){
    treat[i,1]=4}
  else if (Dong2013$treatment[i]=="TIO-HH"){
    treat[i,1]=5}
  else if (Dong2013$treatment[i]=="TIO-SMI"){
    treat[i,1]=6}
}

Dong2013$treat=treat

###### Bayesian models ######################3
data_bayes=Dong2013[,-2]
names(data_bayes)=c("study","responders","sampleSize","treatment")



mtc.network <- mtc.network(data.ab =data_bayes, description = "Network")


##### Random-effects model (d~N(0,10000), tau~U(0,2)) 

mtc.model1 <-mtc.model(mtc.network, type ="consistency", om.scale = 2,
                       hy.prior=mtc.hy.prior("std.dev","dunif", 0, 2),linearModel = "random",n.chain = 2,re.prior.sd = 100)


mtc.model1$inits[[1]]$.RNG.name <- "base::Mersenne-Twister"
mtc.model1$inits[[2]]$.RNG.name <- "base::Mersenne-Twister"

mtc.model1$inits[[1]]$.RNG.seed <- 58982
mtc.model1$inits[[2]]$.RNG.seed <- 58983
mtc.run1 <- mtc.run(mtc.model1, thin = 1,n.iter =50000,n.adapt = 10000)

s1=summary(mtc.run1)$summaries[1]
s2=summary(mtc.run1)$summaries[2]
s1=as.data.frame(s1)
s2=as.data.frame(s2)

###### Printing results in OR scale
print(round((exp(s1[1:nrow(s1)-1,1])),digits = 2)) ### treatments effects
print(round((exp(s2[1:nrow(s2)-1,1])),digits = 2)) ## lower bounds of credible interval
print(round((exp(s2[1:nrow(s2)-1,5])),digits = 2)) ## upper bounds of credible interval


tau=round((s1[nrow(s1),1]),digits = 2) ## exporting tau
tau_lb=(round((s2[nrow(s2),1]),digits = 2)) ## lower bound for tau
tau_ub=(round((s2[nrow(s2),5]),digits = 2)) ## upper bound for tau

print(tau)
print(tau_lb)
print(tau_ub)

A=gelman.diag(mtc.run1)
############# Random-effects model d~N(0,100), tau~HN(1) #############
hy.prior <- mtc.hy.prior(type="std.dev", distr="dhnorm", 0, 1)

mtc.model11 <-mtc.model(mtc.network, type ="consistency", om.scale = 2,
                       hy.prior=hy.prior,linearModel = "random",n.chain = 2,re.prior.sd = 10)


mtc.model11$inits[[1]]$.RNG.name <- "base::Mersenne-Twister"
mtc.model11$inits[[2]]$.RNG.name <- "base::Mersenne-Twister"

mtc.model11$inits[[1]]$.RNG.seed <- 58982
mtc.model11$inits[[2]]$.RNG.seed <- 58983
mtc.run11 <- mtc.run(mtc.model11, thin = 1,n.iter =50000,n.adapt = 10000)

s11=summary(mtc.run11)$summaries[1]
s21=summary(mtc.run11)$summaries[2]
s11=as.data.frame(s11)
s21=as.data.frame(s21)

###### Printing results in OR scale
print(round((exp(s11[1:nrow(s11)-1,1])),digits = 2)) ### treatments effects
print(round((exp(s21[1:nrow(s21)-1,1])),digits = 2)) ## lower bounds of credible interval
print(round((exp(s21[1:nrow(s21)-1,5])),digits = 2)) ## upper bounds of credible interval


tau1=round((s11[nrow(s11),1]),digits = 2) ## exporting tau
tau_lb1=(round((s21[nrow(s21),1]),digits = 2)) ## lower bound for tau
tau_ub1=(round((s21[nrow(s21),5]),digits = 2)) ## upper bound for tau

print(tau1)
print(tau_lb1)
print(tau_ub1)


A1=gelman.diag(mtc.run11)
#######################################################

############### Common-effect (d~N(0,10000))  model
mtc.model2 <-mtc.model(mtc.network, type ="consistency", om.scale = 2,
                       hy.prior=mtc.hy.prior("std.dev","dunif", 0, 2),linearModel = "fixed",n.chain = 2,re.prior.sd = 100)


mtc.model2$inits[[1]]$.RNG.name <- "base::Mersenne-Twister"
mtc.model2$inits[[2]]$.RNG.name <- "base::Mersenne-Twister"


mtc.model2$inits[[1]]$.RNG.seed <- 58982
mtc.model2$inits[[2]]$.RNG.seed <- 58983
mtc.run2 <- mtc.run(mtc.model2, thin = 1,n.iter =50000,n.adapt = 10000)

s12=summary(mtc.run2)$summaries[1]
s22=summary(mtc.run2)$summaries[2]
s12=as.data.frame(s12)
s22=as.data.frame(s22)


###### Printing results in OR scale
print(round((exp(s12[,1])),digits = 2)) ### treatments effects
print(round((exp(s22)),digits = 2)[,1]) ### lower bounds of credible interval
print(round((exp(s22)),digits = 2)[,5]) ### upper bounds of credible interval



##### Brooks-Gelman plots

#### Random effects model
gelman.plot(mtc.run1)

### Common effects model
gelman.plot(mtc.run2)
