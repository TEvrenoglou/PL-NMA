library(gemtc)
library(parallel)

################################################################


################################################################
###### Bayesian models ###########

no_cores1 <- detectCores()-1
cl11 <- makeCluster(no_cores1)

X4=list()
for (i in 1:N.sim){ X4[[i]]=list("data"=data3[[i]],"logOR"=logOR[[i]],"tau"=tau[[1]])}

BL_NMA=function(X)
{
  
  biasBL_NMA=c()
  coverageBL_NMA=c()
  mseBL_NMA=c()
  
  
  biasBL_NMA_Random=c()
  bias_tau=c()
  coverageBL_NMA_Random=c()
  mseBL_NMA_Random=c()
  
  
  
  mtc.network <- mtc.network(data.ab =X$data, description = "Network")
  
  # specify the estimation parameters 
  mtc.model1 <-mtc.model(mtc.network, type ="consistency", om.scale = 2,
                         hy.prior=mtc.hy.prior("std.dev","dunif", 0, 2),linearModel = "fixed",re.prior.sd = 100,n.chain = 2)
  
  
  mtc.model1$inits[[1]]$.RNG.name <- "base::Mersenne-Twister"
  mtc.model1$inits[[2]]$.RNG.name <- "base::Mersenne-Twister"
  
  
  mtc.model1$inits[[1]]$.RNG.seed <- 58982
  mtc.model1$inits[[2]]$.RNG.seed <- 58983
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
  
  mseBL_NMA=biasBL_NMA*biasBL_NMA+BL_res[,2]*BL_res[,2]
  
  coverageBL_NMA=c(coverageBL_NMA, BL_res$cover)
  
  ########################################################################################
  ########################################################################################
  
  mtc.model11 <-mtc.model(mtc.network, type ="consistency", om.scale = 2,
                          hy.prior=mtc.hy.prior("std.dev","dunif", 0, 2),linearModel = "random",re.prior.sd = 100,n.chain = 2)
  
  
  mtc.model11$inits[[1]]$.RNG.name <- "base::Mersenne-Twister"
  mtc.model11$inits[[2]]$.RNG.name <- "base::Mersenne-Twister"
  
  
  mtc.model11$inits[[1]]$.RNG.seed <- 58962
  mtc.model11$inits[[2]]$.RNG.seed <- 58963
  mtc.run11 <- mtc.run(mtc.model11, thin = 1,n.iter = 50000,n.adapt = 10000)
  
  s11=summary(mtc.run11)$summaries[1]
  s21=summary(mtc.run11)$summaries[2]
  s11=as.data.frame(s11)
  s21=as.data.frame(s21)
  res1=cbind(s11[,1],s11[,2],s21[,1],s21[,5])
  
  
  
  tau1=res1[length(res1[,1]),1]
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
  
  mseBL_NMA_Random=biasBL_NMA_Random*biasBL_NMA_Random+BL_res1[,2]*BL_res1[,2]
  
  coverageBL_NMA_Random=c(coverageBL_NMA_Random, BL_res1$cover)
  
  
  
  
  return(list("bias_fixed"=biasBL_NMA,"cov_fixed"=coverageBL_NMA,"mse_fixed"=mseBL_NMA,
              "bias_random"=biasBL_NMA_Random,"cov_random"=coverageBL_NMA_Random,"mse_random"=mseBL_NMA_Random,"bias_tau"=bias_tau,
              "tau"=tau1,"tau_test"=tau_test))
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

tau_test=c()
for (i in 1:N.sim){tau_test=c(tau_test, l10[[i]]$tau_test)}


coverageBL_NMA_random=c()
for (i in 1:N.sim){  coverageBL_NMA_random=c(coverageBL_NMA_random, l10[[i]]$cov_random)}
mean(coverageBL_NMA_random)

mseBL_NMA_random=c()
for (i in 1:N.sim){  mseBL_NMA_random=c(mseBL_NMA_random, l10[[i]]$mse_random)}
mean(mseBL_NMA_random)

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

