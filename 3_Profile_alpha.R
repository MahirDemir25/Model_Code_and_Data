#CORES<-2 ## 14 update for flux require(doParallel) registerDoParallel(CORES)
require(doParallel)
#registerDoParallel(CORES)
registerDoParallel(cores=as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")[1]))

source("Daphnia.R")

estpars <- setdiff(names(params.init),c("alpha","prey.0","noise.0","error_count.0"))
theta.t <- partrans(Daphnia_pomp,params.init,"toEst")
theta.t.hi <- theta.t.lo <- theta.t 
theta.t.lo[estpars] <- theta.t[estpars]-log(2) 
theta.t.hi[estpars] <- theta.t[estpars]+log(2)

profile_design( alpha=seq(from=-3.5,to=-1,length=25), # range of param we get profile likelihood
               lower=theta.t.lo, upper=theta.t.hi,nprof=100  # nprof is number of jobs 100
) -> pd
dim(pd)
pd <- as.data.frame(t(partrans(Daphnia_pomp,t(pd),"fromEst"))) 

bake("alpha-profile1.rds",{

  foreach (p=iter(pd,"row"), 
           .combine=rbind,
           .errorhandling="remove",
           .inorder=FALSE, .options.mpi=list(chunkSize=1,seed=1598260027L,info=TRUE)
  ) %dopar% {

    require(magrittr) 
    require(plyr) 
    require(reshape2) 
    require(pomp)
    
    options(stringsAsFactors=FALSE) 
    dat<-subset(zoopsubset,select=c("Bosmina_longirostris","day"))
  #   p=pd[1,]
    dat %>% pomp(
      times="day",
      t0=0,
      params=params.init,
      rprocess = euler(step.fun = Daphnia_rprocess, delta.t=1), 
      rmeasure= Daphnia_rmeasure,
      dmeasure = Daphnia_dmeasure,
      covar=covartable,
      obsnames = c("Bosmina_longirostris"),
      accumvars = c("error_count"),
      statenames = c("prey","noise","error_count"),
      paramnames = c("taudem","tauenv","taudem2","tauenv2","log.betaz1","alpha","FF","IP2","IP","cc","mu","eta","sdbeta","meanp","sigmap","prey.0","noise.0","error_count.0"),
      covarnames = c("seas_1","byth","nocoup","pulse","JD","Ave_Temp_diff","S_45","S_110"),
      partrans=parameter_trans( toEst=Daphnia_untrans, fromEst=Daphnia_trans)
    ) %>%
      mif2(params = unlist(p),
           Nmif = 75, # 75
           rw.sd = rw.sd(taudem=0.02, tauenv=0.02,taudem2=0.02,tauenv2=0.02,log.betaz1=0.02,log.betaz2=0.02,log.betaz3=0.02,log.betaz4=0.02,log.betaz5=0.02, alpha=0.0, IP=0.02,IP2=0.02,FF=0.005,mu=0.005, cc=0.02,eta=0.02,sdbeta=0.02,meanp=0.02,sigmap=0.02),
           Np = 4000, # 4000
           cooling.type = "geometric", 
           cooling.fraction.50 = 0.1, 
           #max.fail=250,
           ) %>%
      mif2() -> mf
    ## Runs 10 particle filters to assess Monte Carlo error in likelihood 
    pf <- replicate(10, pfilter(mf, Np = 4000))  # rep 10 and Np=4000
    ll <- sapply(pf,logLik)
    ll <- logmeanexp(ll, se = TRUE)

    data.frame(as.list(coef(mf)), 
               logLik = ll[1],
               logLik.se = ll[2]
)
  }
}) -> alpha_prof


