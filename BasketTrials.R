#-------------------------------------------------------------------------#
#               Analysing the BRAF-V600 basket trial                      #
#-------------------------------------------------------------------------#
# please modify the working directory to locate the folder on your PC 
# where you have saved your JAGS model scripts
setwd("~/Documents/BayesianClinicalTrials/Code")

set.seed(684324)
library(rjags)

#---------------- Approach 1: Standard hierarchical modelling ----------------#
data_list <- list(
  nMod = 6,
  r = c(8, 1, 6, 0, 2, 6),
  n = c(20, 8, 18, 10, 7, 32),
  Prior.mu.theta = c(0, 10),
  Prior.tau.HN = 0.125,
  p.cut = 0.25
)

inits <- list(
    mu.theta = 0
  )

iter <- 10000

jags <- jags.model(file = "StandardHM.txt", data = data_list, 
                   n.chains = 1,  n.adapt = iter)

update(jags, iter, progress.bar = "none")
mcmc.sampling <- jags.samples(jags, c("p", "success"), iter, progress.bar = "none")

samples.p <- mcmc.sampling$p[1,,]    # change the value in the bracket from 1 to K 
samples.success <- mcmc.sampling$success[1,,]
mean(samples.p); sd(samples.p)
mean(samples.success)

#Or alternatively
posterior_sample <- coda.samples(jags, data = data_list, 
                                 n.chains = 1,  variable.names = c("p", "success"), n.iter = iter)
summary(posterior_sample)




#-------------------- Approach 2: EXNEX --------------------#
data_list <- list(
  nMod = 6,
  r = c(8, 1, 6, 0, 2, 6),
  n = c(20, 8, 18, 10, 7, 32),
  Prior.mu.theta = c(0, 10),
  Prior.tau.HN = 0.5,
  p.cut = 0.125,
  pMix = c(0.5, 0.5),
  nex.theta = 0,
  nex.sig = 10
)

inits <- list(
  mu.theta = 0
)

iter <- 10000

jags <- jags.model(file = "EXNEX.txt", data = data_list, 
                   n.chains = 1,  n.adapt = iter)

update(jags, iter, progress.bar = "none")
mcmc.sampling <- jags.samples(jags, c("p", "success"), iter, progress.bar = "none")

samples.p <- mcmc.sampling$p[1,,]# change the value in the bracket from 1 to K 
samples.success <- mcmc.sampling$success[1,,]

mean(samples.p); sd(samples.p)
mean(samples.success)


#-------------------- Approach 3: Stand-alone analysis --------------------#
data_list <- list(
  nMod = 6,
  r = c(8, 1, 6, 0, 2, 6),
  n = c(20, 8, 18, 10, 7, 32),
  nex.theta = 0,
  nex.sig = 10,
  p.cut = 0.25
)


inits <- list(
  theta = rep(0, 6)
)

iter <- 10000
jags <- jags.model(file = "Stand-alone.txt", data = data_list, 
                   n.chains = 1,  n.adapt = iter)

update(jags, iter, progress.bar = "none")
mcmc.sampling <- jags.samples(jags, c("p", "success"), iter, progress.bar = "none")

samples.p <- mcmc.sampling$p[1,,]# change the value in the bracket from 1 to K 
samples.success <- mcmc.sampling$success[1,,]

mean(samples.p); sd(samples.p)
mean(samples.success)

# Alternatively
mcmc.sampling0 <- coda.samples(jags, data = data_list, n.chains = 1,
                               variable.names = c("p", "success"), n.iter = 1e4)

mcmc_stats0 <- summary(mcmc.sampling0)
mcmc_stats0$statistics[1:6,1]
