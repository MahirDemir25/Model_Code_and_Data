#CORES<-10 ##update for flux 10
JOBS<-100 ##update for flux 100

require(doParallel)
#registerDoParallel(cores=1)
registerDoParallel(cores=as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")[1]))


tic <- Sys.time()

mpar <- foreach(
  i=1:JOBS,
  .packages=c('pomp'),
  .inorder=FALSE) %dopar% {
    Sys.sleep(i*.1)

    NMIF<-200 ## update for flux 200
    NP<-40000 ## Number of particles 40000
    METHOD="mif2"
    source("Bosmina_in_SSM.R")
    param.tab <- read.table("params.csv", sep=",",row.names=1, header=TRUE)
    LV.pars <- c("taudem","tauenv","taudem2","tauenv2","log.betaz1","log.betaz2","log.betaz3","log.betaz4","log.betaz5","IP","alpha","FF","IP2","mu","cc","eta","sdbeta","meanp","sigmap")
    LV.ivps <- c("prey.0")
    LV.rw.sd<- rw.sd(taudem=0.02, tauenv=0.02,taudem2=0.02,tauenv2=0.02,log.betaz1=0.02,log.betaz2=0.02,log.betaz3=0.02,log.betaz4=0.02,log.betaz5=0.02, alpha=0.005, IP=0.02,IP2=0.02,FF=0.005,mu=0.02, cc=0.02,eta=0.02,sdbeta=0.02,meanp=0.02,sigmap=0.02)

    LV.hyperparams <-
      list(min=unlist(param.tab["lower.bound",]),max=unlist(param.tab["upper.bound",]))
    
    LV.rprior <- function(hyperparams, ...)
    {
      r<-runif(length(hyperparams$min),min=hyperparams$min,max=hyperparams$max)
      names(r) <- names(hyperparams$min)
      return(r)
    }
    set.seed(8100+i)
    Sys.sleep(i*0.1)
    th.draw <-LV.rprior(LV.hyperparams)
    m<-try(mif2(Daphnia_pomp,
               Nmif=NMIF,
               params=th.draw,
               rw.sd=LV.rw.sd,
               Np=NP,
               cooling.type='geometric',
               cooling.fraction= 0.3,
               
    ))
    list(pomp=m,params=th.draw)
  }

m.out <- rbind(
      pf.lik = sapply(mpar,function(x){
        if(class(x$pomp)=="mif2d_pomp") logLik(x$pomp) else NA
      }),
      sapply(mpar,function(x) {
        if(class(x$pomp)=="mif2d_pomp") coef(x$pomp) else rep(NA,length(coef(Daphnia_pomp)))
      }),
      sapply(mpar,function(x)x$params)
    )

    save(m.out,mpar,file="out.rda")
toc <- Sys.time()
print(toc-tic)
print(m.out[1,])
