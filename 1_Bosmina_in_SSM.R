
###### Nonlinear temperature effect on NCE, function(Temp)=B*exp(- (Temp-q)*(Temp-q)/d*d )######
#install.packages("~/Downloads/pomp_2.2.tar",repos = NULL, type = "source")

require(pomp) 
require(lubridate)
require(timeSeries)
Bosmina_rprocess <- Csnippet("
                            double births;
                            double deaths;
                            double betaz;
                            double dw;

                            if (error_count != 0.0) return;
                            
                            betaz = exp(dot_product(5,&seas_1,&log_betaz1));

                            dw = rnorm(0,sqrt(dt));  // white noise for prey birth

                            births = betaz*(exp(-eta*byth*(1/(1+ exp(-1*(IP2- JD)))) ))*dt*prey*(1-prey/cc) + prey*sdbeta*dw + pulse*rlnorm(meanp,sigmap);	// prey births
                            deaths = prey*(alpha*byth*(1/(1+ exp(-1*(IP- JD))))*(exp( -byth*(FF*Ave_Temp_diff )) ) + mu )*dt;	// prey deaths due to predation

                            prey += births - deaths;
                            noise += dw;

                            // check for violations of positivity constraints
                            // nonzero error_count variable signals violation
                            if (prey <= 0.0)  {
                            prey = 0;
                            }
                            if (nocoup < 1)  {
                            prey = 0;
                             }
                            ")

params.init <- c(taudem= 31.23594,tauenv=0.3414941,taudem2= 31.23594,tauenv2=0.3414941,log.betaz1=-3.946131,log.betaz2=-2.085238,log.betaz3=-1.244633,log.betaz4=-1.254834,log.betaz5=-1.254834,alpha=0.01,FF=0.01,eta=0.1,IP=320, IP2=240, mu= 0.1282482,meanp=6.447619,sigmap= 1.693872,cc=451366.3,sdbeta=0.2263157,prey.0=0,noise.0=0,error_count.0=0)


# Import combined 45m and 110m data set: (1)We consider data between years 1998-2018 at 45m site and between 1994-2018 for 110m site. (2) We remove year 2000 in both site since this year has extremely high Bosmina abundance.
zoop<-read.csv("zoop45m_plus_110m.csv",header=TRUE)

zoop<-subset(zoop, Year!=2000)
zoop<-subset(zoop, day!=4500);zoop<-subset(zoop, day!=4509);zoop<-subset(zoop, day!=4520) # we have 3 zero observation after this first observation in year 2010
zoop<-subset(zoop, day!=14621) #remove anomalous January observation for 110m data set
zoop<-subset(zoop, day!=15113) #remove anomalous observation in year 2014 at 110m data set
zoop<-subset(zoop, Year!=2025) # year 2000 at 110m data set
#zs <-subset(zoop, select=c("DOY","Year","day", "Bythotrephes_longimanus", "Daphnia_galeata_mendotae", "Bosmina_longirostris"))

# rescale Bythotrephes Biomass (#ind/m2) divding by its sd
zoop$Byth_rs<-zoop$Bythotrephes_longimanus/372.7
sort(zoop$Byth_rs)
#generate Bythotrephes covariate
Byth<-rep(c(rep(0,50),rep(NA,315)),46); 
Byth[zoop$day]<-zoop$Byth_rs; Byth[1]<-0;Byth[16790]<-0
Byth<-interpNA(Byth, method = "linear")
sort(Byth)
#calculate 45 day moving average
Byth_corrected<-Byth+0.00340706
BythMA45<-rep(0.00340706,16790)

for (i in 23:16767) {
  BythMA45[i]<-(Byth_corrected[i-22]*Byth_corrected[i-21]*Byth_corrected[i-20]*Byth_corrected[i-19]*Byth_corrected[i-18]*Byth_corrected[i-17]*
                  Byth_corrected[i-16]*Byth_corrected[i-15]*Byth_corrected[i-14]*Byth_corrected[i-13]*Byth_corrected[i-12]*Byth_corrected[i-11]*
                  Byth_corrected[i-10]*Byth_corrected[i-9]*Byth_corrected[i-8]*
                  Byth_corrected[i-7]*Byth_corrected[i-6]*Byth_corrected[i-5]*Byth_corrected[i-4]*Byth_corrected[i-3]*Byth_corrected[i-2]*
                  Byth_corrected[i-1]*Byth_corrected[i]*Byth_corrected[i+1]*Byth_corrected[i+2]*Byth_corrected[i+3]*Byth_corrected[i+4]*
                  Byth_corrected[i+5]*Byth_corrected[i+6]*Byth_corrected[i+7]*Byth_corrected[i+8]*Byth_corrected[i+9]*
                  Byth_corrected[i+10]*Byth_corrected[i+11]*Byth_corrected[i+12]*Byth_corrected[i+13]*Byth_corrected[i+14]*
                  Byth_corrected[i+15]*Byth_corrected[i+16]*Byth_corrected[i+17]*Byth_corrected[i+18]*Byth_corrected[i+19]*Byth_corrected[i+20]*Byth_corrected[i+21]*Byth_corrected[i+22])^(1/45)
}

#need to back-transform
BythMA45<-BythMA45-0.00340706;BythMA45<-round(BythMA45,10)

## Julian days as a covariate
jdays<-1:365
JD<-rep(jdays,46)

#setting population zero at the end of each season
nocoup<-rep(c(rep(1,364),0),46)


#import temperature data for #### 45m #####
Tempdata_merge<-read.csv("Tempdata_final_full_45m.csv",head=T)
Tempdata_merge<-subset(Tempdata_merge,Year!=1994)
Tempdata_merge<-subset(Tempdata_merge,Year!=1995)
Tempdata_merge<-subset(Tempdata_merge,Year!=1996)
Tempdata_merge<-subset(Tempdata_merge,Year!=1997)
Tempdata_merge<-subset(Tempdata_merge,Year!=2000)
Tempdata_merge$day<-as.numeric(1+(as.Date(Tempdata_merge$Date, format="%m/%d/%y")-rep(as.Date("1998-01-01"),length(Tempdata_merge$Date)) ))

#generate surface temperature covariate
Temp<-rep(c(rep(0,100),rep(NA,265)),21);
Temp[Tempdata_merge$day]<-Tempdata_merge$WTMP; Temp[1]<-0;Temp[7665]<-0
Temp<-interpNA(Temp, method = "linear")
Temp<-Temp[1:7665]


## Temp at 40m in Lake Michigan (LM) for 45m
L_Temp<-read.csv("Temp_at40m_98-18_for45m.csv",header=TRUE)
L_Temp<-subset(L_Temp,Year!=2000)
L_Temp$day<-as.numeric((as.Date( L_Temp$Date, format="%m/%d/%y")-rep(as.Date("1998-01-01"),length(L_Temp$Date)) ))

## covariate of Temp at 40m
l_temp<-rep(c(rep(0,86),rep(NA,279)),21);
l_temp[L_Temp$day]<-L_Temp$Temp_at40m; l_temp[1]<-0;l_temp[7665]<-0
l_temp<-interpNA(l_temp, method = "linear")

Temp_diff45=Temp-l_temp

#replace any negative value with 0
for(i in 1:nrow(Temp_diff45))
{
  if (Temp_diff45[i,1]<0)
    Temp_diff45[i,1]=0
}


#import surface temperature data for #### 110m ######
Tempdata_merge<-read.csv("Tempdata_final_94-18.csv",head=T,sep=";", dec=",")
Tempdata_merge$day<-as.numeric(1+(as.Date(Tempdata_merge$Date)-rep(as.Date("1994-01-01"),length(Tempdata_merge$Date))))

#generate surface temperature covariate
Temp<-rep(c(rep(0,100),rep(NA,265)),25); 
Temp[Tempdata_merge$day]<-Tempdata_merge$WTMP; Temp[1]<-0;Temp[9125]<-0
Temp<-interpNA(Temp, method = "linear")

## Temp at 40m in Lake Michigan (LM) for 110m
L_Temp<-read.csv("Temp_at40m_94-18_full.csv",header=TRUE,sep=",")
L_Temp$day<-as.numeric((as.Date( L_Temp$Date, format="%m/%d/%y")-rep(as.Date("1994-01-01"),length(L_Temp$Date)) ))
## covariate of Temp at 40m
l_temp<-rep(c(rep(0,86),rep(NA,279)),25);
l_temp[L_Temp$day]<-L_Temp$Temp_at40m; l_temp[1]<-0;l_temp[9125]<-0
l_temp<-interpNA(l_temp, method = "linear")
#L_Temp$Date<-as.Date(L_Temp$DOY, origin=as.Date("1994-01-01"))

## Temp difference (Temperature gradient)
Temp_diff110m=Temp-l_temp
#replace any negative value with 0
for(i in 1:nrow(Temp_diff110m))
{
  if (Temp_diff110m[i,1]<0)
    Temp_diff110m[i,1]=0
}

Temp_diff<-c(Temp_diff45,Temp_diff110m)
Tempp=Temp_diff
# obtain the avaraged julian year for Temp data
Temp1998<-Tempp[1:365]
Temp1999<-Tempp[366:730]
Temp2000<-Tempp[731:1095]
Temp2001<-Tempp[1096:1460]
Temp2002<-Tempp[1461:1825]
Temp2003<-Tempp[1826:2190]
Temp2004<-Tempp[2191:2555]
Temp2005<-Tempp[2556:2920]
Temp2006<-Tempp[2921:3285]
Temp2007<-Tempp[3286:3650]

Temp2008<-Tempp[3651:4015]
Temp2009<-Tempp[4016:4380]
Temp2010<-Tempp[4381:4745]

Temp2011<-Tempp[4746:5110]
Temp2012<-Tempp[5111:5475]
Temp2013<-Tempp[5476:5840]
Temp2014<-Tempp[5841:6205]
Temp2015<-Tempp[6206:6570]
Temp2016<-Tempp[6571:6935]
Temp2017<-Tempp[6936:7300]
Temp2018<-Tempp[7301:7665]

# obtain the avaraged for Byth data vs julian year for 110 m site
Temp2019<-Tempp[7666:8030]     # corresponds year 1994 in 110m site
Temp2020<-Tempp[8031:8395]     # corresponds year 1995 in 110m site
Temp2021<-Tempp[8396:8760]     # corresponds year 1996 in 110m site
Temp2022<-Tempp[8761:9125]     # corresponds year 1997 in 110m site
Temp2023<-Tempp[9126:9490]     # corresponds year 1998 in 110m site
Temp2024<-Tempp[9491:9855]     # corresponds year 1999 in 110m site
Temp2025<-Tempp[9856:10220]    # corresponds year 2000 in 110m site 
Temp2026<-Tempp[10221:10585]   # corresponds year 2001 in 110m site
Temp2027<-Tempp[10586:10950]   # corresponds year 2002 in 110m site
Temp2028<-Tempp[10951:11315]   # corresponds year 2003 in 110m site

Temp2029<-Tempp[11316:11680]   # corresponds year 2004 in 110m site
Temp2030<-Tempp[11681:12045]   # corresponds year 2005 in 110m site
Temp2031<-Tempp[12046:12410]   # corresponds year 2006 in 110m site

Temp2032<-Tempp[12411:12775]   # corresponds year 2007 in 110m site
Temp2033<-Tempp[12776:13140]   # corresponds year 2008 in 110m site
Temp2034<-Tempp[13141:13505]   # corresponds year 2009 in 110m site
Temp2035<-Tempp[13506:13870]   # corresponds year 2010 in 110m site
Temp2036<-Tempp[13871:14235]   # corresponds year 2011 in 110m site
Temp2037<-Tempp[14236:14600]   # corresponds year 2012 in 110m site
Temp2038<-Tempp[14601:14965]   # corresponds year 2013 in 110m site
Temp2039<-Tempp[14966:15330]   # corresponds year 2014 in 110m site
Temp2040<-Tempp[15331:15695]   # corresponds year 2015 in 110m site
Temp2041<-Tempp[15696:16060]   # corresponds year 2016 in 110m site
Temp2042<-Tempp[16061:16425]   # corresponds year 2017 in 110m site
Temp2043<-Tempp[16426:16790]   # corresponds year 2018 in 110m site


# mean of temperature gradient for 38 years with excluding the years 2000, 2004-2006 from both sites(45m and 110m)
Temp_diff_m<- rowMeans(cbind(Temp1998,Temp1999,Temp2001,Temp2002,Temp2003,Temp2007,Temp2008,Temp2009,Temp2010,Temp2011,Temp2012,Temp2013,Temp2014,Temp2015,Temp2016,Temp2017,Temp2018,
                             Temp2019,Temp2020,Temp2021,Temp2022,Temp2023,Temp2024,Temp2026,Temp2027,Temp2028,Temp2032,Temp2033,Temp2034,Temp2035,Temp2036,Temp2036,Temp2038,Temp2039,Temp2040,Temp2041,Temp2042,Temp2043))
Ave_Temp_diff<-rep(Temp_diff_m,46) # we obtained average temperature gradient and used for 46 years (Actually 38 years since years 2000, 2004-2006 are excluded from both sites in this calculation) 

#determine date and value of first date of Bosmina observation each year for pulse of Bosmina
zoop$year<-zoop$Year
surveyyears<-c(1998,1999,2001:2003,2007:2024,2026:2028,2032:2043)
firstdaph<-rep(0,38)
firstdaphjday<-rep(0,38)
firstdaphdens<-rep(0,38)

for (i in 1:38) {
  zoopyear<-subset(zoop,year==surveyyears[i])
  zoopyearno0<-subset(zoopyear,Bosmina_longirostris>0)
  firstdaph[i]<-zoopyearno0$day[which.min(zoopyearno0$DOY)]
  firstdaphjday[i]<-zoopyearno0$DOY[which.min(zoopyearno0$DOY)]
  firstdaphdens[i]<-zoopyearno0$Bosmina_longirostris[which.min(zoopyearno0$DOY)]
}

daphstart<-firstdaph-7 

#create vector with 1s at time of first observation
pulse<-rep(0,16790)
pulse[daphstart]<-1

#need to generate subset of Daphnia observations excluding early year zero observations
rem.days<-0
for (i in 1:38) {
  zoopyear<-subset(zoop,year==surveyyears[i])
  zoopyear0<-subset(zoopyear,day<firstdaph[i])
  rem.days<-c(rem.days,zoopyear0$day)
}
rem.days<-rem.days[2:65]
zoopsubset<-zoop
for (i in 1:64) {
  zoopsubset<-subset(zoopsubset,day!=rem.days[i])
}

# We set Swich for measurement error paramaters for 45m and 110m
S_45<-c(rep(1,7665), rep(0,9125))
S_110<-c(rep(0,7665), rep(1,9125))

# covariate table for pomp 
covartable <- covariate_table(
  t=seq(from=0,to=16789,by=1),
  seas=periodic_bspline_basis(t,nbasis=5,degree=3,period=365),
  byth=BythMA45,
  nocoup= nocoup,
  S_45=S_45,
  S_110=S_110,
  Ave_Temp_diff=Ave_Temp_diff,
  JD= JD,
  pulse=pulse,
  times="t"
)


Bosmina_rmeasure <- Csnippet("
                                 
                                  double tau, bosmina;
                                  tau=sqrt(taudem*taudem*prey*S_45 + tauenv*tauenv*prey*prey*S_45 +taudem2*taudem2*prey*S_110 + tauenv2*tauenv2*prey*prey*S_110);
                                  if ((error_count > 0.0) || (!(R_FINITE(prey)))) {
                                  Bosmina_longirostris = R_NaReal;
                                  } else {
                                  bosmina = rnorm(prey,tau);
                                  if (bosmina<=0) {
                                  Bosmina_longirostris=0;
                                  }
                                  else {
                                  Bosmina_longirostris=bosmina;
                                  }
                                  }
                                  ")

Bosmina_dmeasure <- Csnippet("
                                   double tau, tol = 1.0e-18;
                                  tau=sqrt(taudem*taudem*prey*S_45 + tauenv*tauenv*prey*prey*S_45 +taudem2*taudem2*prey*S_110 + tauenv2*tauenv2*prey*prey*S_110);
                                  double f = 0.0;
                                  if ((error_count > 0.0) || (!(R_FINITE(prey)))) {
                                  lik = tol;
                                  } else {
                                  if (Bosmina_longirostris==0) {
                                  f += pnorm(0,prey,tau,1,1)+tol;
                                  }
                                  else {
                                  f += dnorm(Bosmina_longirostris, prey, tau, 1)+tol;
                                  }
                                  lik = (give_log) ? f : exp(f);
                                  }
                                  ")


Bosmina_untrans <- Csnippet("
                       T_taudem = log(taudem);
                       T_tauenv = log(tauenv);
                       T_taudem2 = log(taudem2);
                       T_tauenv2 = log(tauenv2);
                       T_cc = log(cc);
                       T_sdbeta = log(sdbeta);
                       T_mu = log(mu);
                       T_sigmap = log(sigmap);
                       T_eta = log(eta);
                       T_alpha = log(alpha);
                       T_IP = log(IP);
                       T_IP2 = log(IP2);
                      
                       ")

Bosmina_trans <- Csnippet("
                     taudem = exp(T_taudem);
                     tauenv = exp(T_tauenv);
                     taudem2 = exp(T_taudem2);
                     tauenv2 = exp(T_tauenv2);
                     cc = exp(T_cc);
                     sdbeta = exp(T_sdbeta);
                     mu = exp(T_mu);
                     sigmap = exp(T_sigmap);
                     eta = exp(T_eta);
                     alpha = exp(T_alpha);
                     IP = exp(T_IP);
                     IP2 = exp(T_IP2);
                
                     ")

# pomp object for the package pomp 2
Bosmina_pomp <- pomp(
  data=subset(zoopsubset,select=c("Bosmina_longirostris","day")),
  times="day",
  t0=0,
  params=params.init,
  rprocess = euler(step.fun = Bosmina_rprocess, delta.t=1), #  name changed as euler, euler.sim is removed
  rmeasure= Bosmina_rmeasure,
  dmeasure = Bosmina_dmeasure,
  covar=covartable,  # See covartable above, small changes in covar, and do not use tcovar anymore
  obsnames = c("Bosmina_longirostris"),
  accumvars = c("error_count"), # zeronames ---> accumvars
  statenames = c("prey","noise","error_count"),
  paramnames = c("taudem","tauenv","taudem2","tauenv2","log.betaz1","alpha","FF","IP2","IP","cc","mu","eta","sdbeta","meanp","sigmap","prey.0","noise.0","error_count.0"),
  covarnames = c("seas_1","byth","nocoup","pulse","JD","Ave_Temp_diff","S_45","S_110"),
  partrans=parameter_trans( toEst=Bosmina_untrans, fromEst=Bosmina_trans) 
)

Bosmina_sim<-simulate(Bosmina_pomp)
#plot(Bosmina_sim)
#pf<-pfilter(Bosmina_pomp,params=params.init,Np=2000)
#pf@loglik
#plot(Daphnia_sim, variables=c("Bosmina_longirostris"))


