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
  
  
  biasBL_NMA_Random=c()
  bias_tau=c()
  coverageBL_NMA_Random=c()
  mseBL_NMA_Random=c()
  se_BL_NMA_Random=c()
  length_BL_NMA_Random=c()
  lower_BL_NMA_Random=c()
  upper_BL_NMA_Random=c()
  estimate_BL_NMA_Random=c()
  
  
  
  lower_tau=c()
  upper_tau=c()
  length_tau=c()
  se_tau=c()
  
  
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
  
  
  mtc.model11$inits[[1]]$.RNG.name <- "base::Mersenne-Twister"
  mtc.model11$inits[[2]]$.RNG.name <- "base::Mersenne-Twister"
  #mtc.model11$inits[[3]]$.RNG.name <- "base::Mersenne-Twister"
  #mtc.model11$inits[[4]]$.RNG.name <- "base::Mersenne-Twister"
  
  mtc.model11$inits[[1]]$.RNG.seed <- 58962
  mtc.model11$inits[[2]]$.RNG.seed <- 58963
  #mtc.model11$inits[[3]]$.RNG.seed <- 58964
  #mtc.model11$inits[[4]]$.RNG.seed <- 58965
  mtc.run11 <- mtc.run(mtc.model11, thin = 1,n.iter = 50000,n.adapt = 10000)
  
  s11=summary(mtc.run11)$summaries[1]
  s21=summary(mtc.run11)$summaries[2]
  s11=as.data.frame(s11)
  s21=as.data.frame(s21)
  res1=cbind(s11[,1],s11[,2],s21[,1],s21[,5])
  
  
  
  tau1=res1[length(res1[,1]),1]
  
  lower_tau=c(lower_tau,res1[length(res1[,1]),3])
  upper_tau=c(upper_tau,res1[length(res1[,1]),4])
  se_tau=c(se_tau,res1[length(res1[,1]),2])
  length_tau=c(length_tau,res1[length(res1[,1]),4]-res1[length(res1[,1]),3])
  
  tau_test=res1[length(res1[,1]),1]
  
  if(tau1<0 | tau1==0){
    tau1=0}
  if(tau1>0){
    tau1=tau1
  }
  res1=res1[-c(length(res1[,1])),]
  
  BL_res1=as.data.frame(res1)
  names(BL_res1)=c("mean","se","lowerCR","upperCR")
  
  biasBL_NMA_Random=c(biasBL_NMA_Random,(BL_res1$mean-X$logOR))
  bias_tau=c(bias_tau,(tau1-X$tau))
  BL_res1$cover=(BL_res1$lowerCR<X$logOR)&(BL_res1$upperCR>X$logOR)
  
  mseBL_NMA_Random=c(mseBL_NMA_Random,biasBL_NMA_Random*biasBL_NMA_Random+BL_res1[,2]*BL_res1[,2])
  
  coverageBL_NMA_Random=c(coverageBL_NMA_Random, BL_res1$cover)
  
  estimate_BL_NMA_Random=c(estimate_BL_NMA_Random,BL_res1$mean)
  
  se_BL_NMA_Random=c(se_BL_NMA_Random,BL_res1$se)
  
  lower_BL_NMA_Random=c(lower_BL_NMA_Random,BL_res1$lowerCR)
  upper_BL_NMA_Random=c(upper_BL_NMA_Random,BL_res1$upperCR)
  length_BL_NMA_Random=c(length_BL_NMA_Random,BL_res1$upperCR-BL_res1$lowerCR)
  
  
  
  return(list("bias_fixed"=biasBL_NMA,"cov_fixed"=coverageBL_NMA,"mse_fixed"=mseBL_NMA,
              "estimate_fixed_bayes"=estimate_BL_NMA,"se_fixed_bayes"=se_BL_NMA,"lower_fixed_bayes"=lower_BL_NMA,"upper_fixed_bayes"=upper_BL_NMA,
              "length_fixed_bayes"=length_BL_NMA,
              "bias_random"=biasBL_NMA_Random,"cov_random"=coverageBL_NMA_Random,"mse_random"=mseBL_NMA_Random,"bias_tau"=bias_tau,
              "tau"=tau1,"tau_test"=tau_test,
              "estimate_random_bayes"=estimate_BL_NMA_Random,"se_random_bayes"=se_BL_NMA_Random,"lower_random_bayes"=lower_BL_NMA_Random,
              "upper_random_bayes"=upper_BL_NMA_Random,"length_random_bayes"=length_BL_NMA_Random,
              "lower_tau"=lower_tau,"upper_tau"=upper_tau,"se_tau"=se_tau,"length_tau"=length_tau
  ))
}




clusterExport(cl11,"X4")
clusterExport(cl11,"NT")
clusterExport(cl11,"N.sim")
clusterExport(cl11, "BL_NMA")
clusterEvalQ(cl11, {library(gemtc)})
l10=parLapply(cl11,1:N.sim, function(x) BL_NMA(X4[[x]]))


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

bias_tau=c()
for (i in 1:N.sim){  bias_tau=c(bias_tau, l10[[i]]$bias_tau)}
mean(bias_tau)

tau=c()
for (i in 1:N.sim){tau=c(tau, l10[[i]]$tau)}

se_tau=c()
for (i in 1:N.sim){se_tau=c(se_tau, l10[[i]]$se_tau)}

upper_tau=c()
for (i in 1:N.sim){upper_tau=c(upper_tau, l10[[i]]$upper_tau)}

lower_tau=c()
for (i in 1:N.sim){lower_tau=c(lower_tau, l10[[i]]$lower_tau)}

length_tau=c()
for (i in 1:N.sim){length_tau=c(length_tau, l10[[i]]$length_tau)}
mean(length_tau)

tau_test=c()
for (i in 1:N.sim){tau_test=c(tau_test, l10[[i]]$tau_test)}



coverageBL_NMA_random=c()
for (i in 1:N.sim){  coverageBL_NMA_random=c(coverageBL_NMA_random, l10[[i]]$cov_random)}
mean(coverageBL_NMA_random)

mseBL_NMA_random=c()
for (i in 1:N.sim){  mseBL_NMA_random=c(mseBL_NMA_random, l10[[i]]$mse_random)}
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

lower_random_bayes=c()
for (i in 1:N.sim){lower_random_bayes=c(lower_random_bayes, l10[[i]]$lower_random_bayes)}

upper_random_bayes=c()
for (i in 1:N.sim){upper_random_bayes=c(upper_random_bayes, l10[[i]]$upper_random_bayes)}

length_random_bayes=c()
for (i in 1:N.sim){length_random_bayes=c(length_random_bayes, l10[[i]]$length_random_bayes)}
mean(length_random_bayes)

estimate_random_bayes=c()
for (i in 1:N.sim){estimate_random_bayes=c(estimate_random_bayes, l10[[i]]$estimate_random_bayes)}







#############################

################################### RESULTS ##################################################


bias=c(round(mean(biasBL_NMA_fixed), digits=2),round(mean(biasBL_NMA_random), digits=2))
absolute_bias=c(round(mean(abs(biasBL_NMA_fixed)), digits=2),round(mean(abs(biasBL_NMA_random)), digits=2))
coverage=c(round(mean(coverageBL_NMA_fixed), digits=3)*100,round(mean(coverageBL_NMA_random), digits=3)*100)
mse=c(round(mean(mseBL_NMA_fixed),digits = 2),round(mean(mseBL_NMA_random),digits = 2))
bias_tau=c(bias_tau)


results=cbind(bias,absolute_bias,mse,coverage)
results=data.frame(results)
row.names(results)=c("BL_NMA_fixed","BL_NMA_random")
results



results.csv=write.csv(results,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Simulation new\\Scenario 31\\Bayesian\\main_results_Bayesian.csv",row.names = TRUE)


bias1=data.frame(biasBL_NMA_fixed,biasBL_NMA_random)
names(bias1)=c("biasBL_NMA_fixed","biasBL_NMA_random")

write.csv(bias1,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Simulation new\\Scenario 31\\Bayesian\\bias_Bayesian.csv",row.names = FALSE)

avg_bias_tau=c(round(mean(bias_tau), digits=2))
avg_bias_tau=data.frame(avg_bias_tau)
names(avg_bias_tau)=c("mean_bias_tau")

write.csv(avg_bias_tau,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Simulation new\\Scenario 31\\Bayesian\\mean_bias_tau.csv",row.names = FALSE)


mse1=data.frame(mseBL_NMA_fixed,mseBL_NMA_random)
names(mse1)=c("mseBL_NMA_fixed","mseBL_NMA_random")
write.csv(mse1,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Simulation new\\Scenario 31\\Bayesian\\mse_Bayesian.csv",row.names = FALSE) 
#############################################################################################################################################
coverage11=data.frame(coverageBL_NMA_fixed,coverageBL_NMA_random)
names(coverage11)=c("coverage.BL_NMA_fixed","coverageBL_NMA_random")
write.csv(coverage11,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Simulation new\\Scenario 31\\Bayesian\\coverage_Bayesian.csv",row.names = FALSE)

tau_info=data.frame(bias_tau,tau,tau_test)
names(tau_info)=c("bias_tau","tau","tau_test")
write.csv(tau_info,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Simulation new\\Scenario 31\\Bayesian\\tau_info_Bayesian.csv",row.names = FALSE)


lengths_bayes=data.frame(length_fixed_bayes,length_random_bayes)
names(lengths_bayes)=c("length.BL_NMA_fixed","lengthBL_NMA_random")
write.csv(lengths_bayes,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Simulation new\\Scenario 31\\Bayesian\\lengths_Bayesian.csv",row.names = FALSE)



mean_lengths_bayes=data.frame(mean(length_fixed_bayes),mean(length_random_bayes))
names(lengths_bayes)=c("length.BL_NMA_fixed","lengthBL_NMA_random")
write.csv(mean_lengths_bayes,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Simulation new\\Scenario 31\\Bayesian\\mean_lengths_Bayesian.csv",row.names = FALSE)


credible_intervals=data.frame(lower_fixed_bayes,upper_fixed_bayes,lower_random_bayes,upper_random_bayes,lower_tau,upper_tau)
names(credible_intervals)=c("lower.BL_NMA_fixed","upper.BL_NMA_fixed","lower.BL_NMA_random","upper.BL_NMA_random",
                            "lower.tau","upper.tau")

write.csv(credible_intervals,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Simulation new\\Scenario 31\\Bayesian\\credible_intervals.csv",row.names = FALSE)





standard_errors_bayes=data.frame(se_fixed_bayes,se_random_bayes,se_tau)
names(standard_errors_bayes)=c("se.BL_NMA_fixed","se.BL_NMA_random","se.tau")

write.csv(standard_errors_bayes,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Simulation new\\Scenario 31\\Bayesian\\standard.errors_Bayesian.csv",row.names = FALSE)


estimates_bayes=data.frame(estimate_fixed_bayes,estimate_random_bayes)
names(estimates_bayes)=c("logOR.BL_NMA_fixed","logOR.BL_NMA_random")

write.csv(estimates_bayes,"C:\\Users\\Theodoros Evrenoglou\\Desktop\\Simulation new\\Scenario 31\\Bayesian\\logOR_Bayesian.csv",row.names = FALSE)


