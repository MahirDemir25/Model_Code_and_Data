require(magrittr)
require(plyr)
require(ggplot2)

#eta profile
alpha_prof<-readRDS("alpha-profile1.rds")
dim(alpha_prof)

alpha_prof %<>%
 # subset(nfail.max==0) %>%
  mutate(alpha=signif(alpha,5)) %>%
  ddply(~alpha,subset,rank(-logLik)<=2)
dim(alpha_prof)

#mu_prof<-mu_prof[1:50,]
dim(alpha_prof)
alpha_prof %>%
  ggplot(aes(x=log(alpha),y=logLik))+
  geom_point()+
  geom_smooth(method="loess") ->md

md + theme_classic() +
theme(axis.text=element_text(size=14),
       axis.title=element_text(size=14,face="bold"))+
labs(title = "Daphnia 110m: NCE Model- Value of alpha=0.0825 (-2.5 in log scale)") +
scale_x_continuous(name="log(alpha)", breaks=c(-5,-3.8,-3,-2.71,-1.95,-1),limits=c(-7,0.5)) +
scale_y_continuous(name="Log-Likelihood", limits=c(-3560, -3546)) +
geom_vline(xintercept=-3.8,linetype="dashed", 
           color = "red", size=1)+
geom_vline(xintercept=-1.95,linetype="dashed", 
           color = "red", size=1)+
geom_vline(xintercept=-2.71,linetype="dashed", 
             color = "grey", size=1)

#scale_y_continuous(breaks=c(0.5,1,2,3,4,5),limits=c(0.5,6))

pairs(~logLik+eta+mu + cc + alpha,
      data=alpha_prof,subset=logLik>max(logLik)-100)

alpha_prof$alpha.1
alpha_prof$logLik
max(alpha_prof$logLik)
max(alpha_prof$logLik)-1.92
