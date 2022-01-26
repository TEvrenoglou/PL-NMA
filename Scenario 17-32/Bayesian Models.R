library(gemtc)
library(parallel)

################################################################


################################################################
###### Bayesian models ###########

no_cores1 <- detectCores()-1
cl11 <- makeCluster(no_cores1)

X4=list()

for (i in 1:N.sim){ X4[[i]]=list("data"=data4[[i]],"logOR"=logOR[[i]],"tau"=tau[[1]])}

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

names(convergence_summary)=c("Common effect","Random-effects with vague priors","Random-effects with informative priors","parameter")
names(convergence_tau)=c("Uniform(0,2)","Half-Normal")

results_convergence=c("Convergence times-Common effect"=length(which(convergence_summary[,1]<1.1)),
                      "Convergence times-Random effects-Vague priors"=length(which(convergence_summary[,2]<1.1)),
                      "Convergence times-Random effects-Informative priors"=length(which(convergence_summary[,3]<1.1)),
                      "Convergence times-tau-Vague priors"=length(which(convergence_tau[,1]<1.1)),
                      "Convergence times-tau-Informative priors"=length(which(convergence_tau[,2]<1.1))
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





#############################

################################### RESULTS ##################################################



bias=c(round(mean(biasBL_NMA_fixed), digits=2),round(mean(biasBL_NMA_random), digits=2),round(mean(biasBL_NMA_random_2), digits=2))
absolute_bias=c(round(mean(abs(biasBL_NMA_fixed)), digits=2),round(mean(abs(biasBL_NMA_random)), digits=2),round(mean(abs(biasBL_NMA_random_2)), digits=2))
coverage=c(round(mean(coverageBL_NMA_fixed), digits=3)*100,round(mean(coverageBL_NMA_random), digits=3)*100,round(mean(coverageBL_NMA_random_2), digits=3)*100)
mse=c(round(mean(mseBL_NMA_fixed),digits = 2),round(mean(mseBL_NMA_random),digits = 2),round(mean(mseBL_NMA_random_2),digits = 2))
bias_tau=cbind.data.frame(bias_tau,bias_tau_2)
names(bias_tau)=c("Uniform","Half_normal")


results=cbind(bias,absolute_bias,mse,coverage)
results=data.frame(results)
row.names(results)=c("BL_NMA_common","BL_NMA_random (d~N(0,10000), tau~U(0,2))","BL_NMA_random (d~N(0,100), tau~HN(1))")
results




bias1=data.frame(biasBL_NMA_fixed,biasBL_NMA_random,biasBL_NMA_random_2)
names(bias1)=c("biasBL_NMA_common","biasBL_NMA_random (d~N(0,10000), tau~U(0,2))","biasBL_NMA_random (d~N(0,100), tau~HN(1))")


avg_bias_tau=c(round(mean(bias_tau$Uniform), digits=2),round(mean(bias_tau$Half_normal), digits=2))
avg_bias_tau=matrix(avg_bias_tau,ncol = 2)
avg_bias_tau=data.frame(avg_bias_tau)
names(avg_bias_tau)=c("Uniform","Half_normal")



mse1=data.frame(mseBL_NMA_fixed,mseBL_NMA_random,mseBL_NMA_random_2)
names(mse1)=c("mseBL_NMA_common","mseBL_NMA_random","mseBL_NMA_random_informative")
#############################################################################################################################################
coverage11=data.frame(coverageBL_NMA_fixed,coverageBL_NMA_random,coverageBL_NMA_random_2)
names(coverage11)=c("coverage.BL_NMA_common","coverageBL_NMA_random (d~N(0,10000), tau~U(0,2))","coverageBL_NMA_random (d~N(0,100), tau~HN(1))")


tau_info=data.frame(bias_tau$Uniform,tau,tau_test,bias_tau$Half_normal,tau_2,tau_test_2)
names(tau_info)=c("bias_tauBL_NMA_random (d~N(0,10000), tau~U(0,2))","tauBL_NMA_random (d~N(0,10000), tau~U(0,2))","tau_testBL_NMA_random (d~N(0,10000), tau~U(0,2))",
                  "bias_tau_BL_NMA_random (d~N(0,100), tau~HN(1))","tau_BL_NMA_random (d~N(0,100), tau~HN(1))","tau_BL_NMA_random (d~N(0,100), tau~HN(1))")



lengths_bayes=data.frame(length_fixed_bayes,length_random_bayes,length_random_bayes_2)
names(lengths_bayes)=c("length.BL_NMA_common","lengthBL_NMA_random (d~N(0,10000), tau~U(0,2))","lengthBL_NMA_random (d~N(0,100), tau~HN(1))")


mean_lengths_bayes=data.frame(mean(length_fixed_bayes),mean(length_random_bayes),mean(length_random_bayes_2))
names(mean_lengths_bayes)=c("length.BL_NMA_common","lengthBL_NMA_random (d~N(0,10000), tau~U(0,2))","lengthBL_NMA_random (d~N(0,100), tau~HN(1))")


credible_intervals=data.frame(lower_fixed_bayes,upper_fixed_bayes,lower_random_bayes,upper_random_bayes,lower_random_bayes_2,upper_random_bayes_2,
                              lower_tau,upper_tau,lower_tau_2,upper_tau_2)
names(credible_intervals)=c("lower.BL_NMA_common","upper.BL_NMA_common",
                            "lower.BL_NMA_random (d~N(0,10000), tau~U(0,2))","upper.BL_NMA_random (d~N(0,10000), tau~U(0,2))",
                            "lower.BL_NMA_random (d~N(0,100), tau~HN(1))","upper.BL_NMA_random (d~N(0,100), tau~HN(1))",
                            "lower.tau_BL_NMA_random (d~N(0,10000), tau~U(0,2))","upper.tau_BL_NMA_random (d~N(0,10000), tau~U(0,2))",
                            "lower.tau_BL_NMA_random (d~N(0,100), tau~HN(1))","upper.tau_BL_NMA_random (d~N(0,100), tau~HN(1))")



standard_errors_bayes=data.frame(se_fixed_bayes,se_random_bayes,se_random_bayes_2,se_tau,se_tau_2)
names(standard_errors_bayes)=c("se.BL_NMA_common","se.BL_NMA_random (d~N(0,10000), tau~U(0,2))","se.BL_NMA_random (d~N(0,100), tau~HN(1))",
                               "se.tau.BL_NMA_random (d~N(0,10000), tau~U(0,2))","se.tau_BL_NMA_random (d~N(0,100), tau~HN(1))")



estimates_bayes=data.frame(estimate_fixed_bayes,estimate_random_bayes,estimate_random_bayes_2)
names(estimates_bayes)=c("logOR.BL_NMA_common","logOR.BL_NMA_random (d~N(0,10000), tau~U(0,2))","logOR.BL_NMA_random (d~N(0,100), tau~HN(1))")
