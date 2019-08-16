# script for analyzing bluehead wrasse otolith data 
# for Bogdan et al. ms (Ethology)

# NOTE: various functions appear at the end of this file; they must be added to the workspace first.


# Looking for patterns of within-group variance in post-settlement growth
require(lme4)
require(ggplot2)
require(MuMIn)
require(gridExtra)
require(nlme)
require(egg)
require(pwr)
require(directlabels)

# Read-in data:
D = read.csv('MASTER_otolith_data.csv')
# Note that data for fish ID 406 were deleted bc of mis-entry

D$month = as.factor(D$month)
D$site = as.factor(D$site)

D$grpsize[D$grpsize==0] = 1
D$ps.rate = D$ps.rate*1e3 # convert to µm
D$grpsize.inv = 1/D$grpsize # competition should be proportional to 1/group size

# Mean post-settlement growth rate within each group
Dm = aggregate(D$ps.rate,by=list(Site=D$site,Month=D$month,
                         Date=D$date,Grpnum=D$grpnum,
                         Grpsize=D$grpsize,Plankton=D$Plankton,Grpsize.inv=D$grpsize.inv),
                        FUN='mean',na.rm=TRUE)
Dm$isGroup = Dm$Grpsize>1

# Mean plankton density experienced by each fish (this is labeled Plankton.hindcast in the datasheet)
Dmp = aggregate(D$Plankton.hindcast,by=list(Site=D$site,Month=D$month,
                                 Date=D$date,Grpnum=D$grpnum,
                                 Grpsize=D$grpsize,Plankton=D$Plankton,Grpsize.inv=D$grpsize.inv),
               FUN='mean',na.rm=TRUE)
Dm$Hindcast.plankton <- Dmp$x


# Var ps growth rate within each group
Dv = aggregate(D$ps.rate,by=list(Site=D$site,Month=D$month,
                                 Date=D$date,Grpnum=D$grpnum,
                                 Grpsize=D$grpsize,Plankton=D$Plankton,Grpsize.inv=D$grpsize.inv),
               FUN='var',na.rm=TRUE)
Dv$isGroup = Dv$Grpsize>1
Dv$std = sqrt(Dv$x)
Dv$Hindcast.plankton <- Dmp$x

# Number of observations within each group (for later filtering)
D$counter= as.numeric(D$ps.rate > 0 )
Dl = aggregate(D$counter,by=list(Site=D$site,Month=D$month,
                                 Date=D$date,Grpnum=D$grpnum,
                                 Grpsize=D$grpsize,Plankton=D$Plankton,Grpsize.inv=D$grpsize.inv),
               FUN='sum',na.rm=TRUE)

# Max ps growth rate
Dmax = aggregate(D$ps.rate,by=list(Site=D$site,Month=D$month,
                                 Date=D$date,Grpnum=D$grpnum,
                                 Grpsize=D$grpsize,Plankton=D$Plankton,Grpsize.inv=D$grpsize.inv),
                 FUN='max',na.rm=TRUE)
Dmax[Dmax==-Inf]=NA # (get rid of spurious Infs created by missing data)
Dmax$isGroup = Dmax$Grpsize>1
Dmax$Hindcast.plankton <- Dmp$x

# Subsets of the data with different numbers of observations:
OK3 = Dl$x >= 3 # only groups with at least 3 observations
OK2 = Dl$x >= 2 # only groups with at least 2 observations
OK4 = Dl$x >= 4

# Subset models for analysis & plotting, based on minimum group memberships
Dm.OK2 = subset(Dm,OK2)
Dm.OK2 = Dm.OK2[!is.na(Dm.OK2$Grpsize),]
Dmax.OK2 = subset(Dmax,OK2)
Dmax.OK2 = Dmax.OK2[!is.na(Dmax.OK2$Grpsize),]
Dv.OK2 = subset(Dv,OK2)
Dv.OK2 = Dv.OK2[!is.na(Dv.OK2$Grpsize),]
Dv.OK2$x <- sqrt(Dv.OK2$x) # variance --> SD
Dv.OK3 = subset(Dv,OK3)
Dv.OK3 = Dv.OK3[!is.na(Dv.OK3$Grpsize),]



######################################################################
# Some preliminary analyses:

# Check for effects of size at emergence & final size on growth rate.
plot(D$AXIS.SAE[!Sonly],D$ps.rate[!Sonly])
plot(D$AXIS.SAE,D$ps.rate,xlab='Otolith axis size-at-emergence (µm)',ylab='Post-settlement growth rate (µm/d)')
sae.m <- lm(ps.rate~AXIS.SAE,data=D)
c.sae.m<-coef(sae.m)
Xdummy = seq(0.19,0.3,length.out=100)
Yline = c.sae.m[1]+c.sae.m[2]*Xdummy
lines(Xdummy,Yline)
plot(D$AXIS,D$ps.rate)
plot(as.factor(D$grpsize>1),D$AXIS)
# There is a weak sig effect of SAE on ps growth, but does not explain much variance.
# plot this:
Xlabel = 'Otolith axis size-at-emergence (µm)'
Ylabel = 'Post-settlement growth rate (µm/d)'
G.obj<-ggplot(D,aes(x=AXIS.SAE,y=ps.rate))+
  geom_smooth(method='lm',color='black',fill="grey80")+
  geom_point(shape=21)+
  ylab(Ylabel)+xlab(Xlabel)+
  theme_bw()+
G.obj


# Raw data pattern. A little misleading bc some of these points are in the same group and some are not, so not independent.
plot(jitter(D$grpsize),D$ps.rate)

plot((D$Plankton.hindcast),D$ps.rate)
# No obvious effect of plankton on growth rate
######################################################################

#############################################################
# Linear models to test for main effects

# Mean growth
m1 = lm(x~ Grpsize.inv+Month+Site+Hindcast.plankton, data=Dm.OK2) # full model
m1 = lm(x~ Grpsize.inv+Month+Hindcast.plankton, data=Dm.OK2)
m1 = lm(x~ Grpsize.inv+Month, data=Dm.OK2) # reduced
summary(m1)


m1a = lm(x~ isGroup + Site + Month+ Hindcast.plankton, data=Dm.OK2) # full model
m1a = lm(x~ isGroup + Month+ Hindcast.plankton, data=Dm.OK2)
m1a = lm(x~ isGroup + Month, data=Dm.OK2) # reduced
summary(m1a)


# Variance in growth
m2 = lm(std ~ Grpsize.inv+Site+Month+Hindcast.plankton, data=Dv.OK3) # full model
summary(m2)

# remove plankton effect (for plotting):
m2b = lm(std ~ Hindcast.plankton, data=Dv.OK3) 
Dv.OK3$Resid = resid(m2b)
m2pR=(lm(Resid~Grpsize.inv+Site+Month,data=Dv.OK3))



# Max growth
m3 = lm(x~ Grpsize.inv + (Site) + (Month)+Hindcast.plankton, data=Dmax.OK2) # full model
m3 = lm(x~ Grpsize.inv + Month + Hindcast.plankton, data=Dmax.OK2)
m3 = lm(x~ Grpsize.inv + Month, data=Dmax.OK2) # reduced
summary(m3) # effect of month


m3.1 = lm(x~ isGroup + Site + Month+Hindcast.plankton, data=Dmax.OK2) # full model
m3.1 = lm(x~ isGroup +  Month+Hindcast.plankton, data=Dmax.OK2)
m3.1 = lm(x~ isGroup + Month, data=Dmax.OK2) # reduced
summary(m3.1)



#----------------
# Plotting

# Add extra column to this data frame so it works with plotting code...
Dv.OK3$x <- Dv.OK3$std 

# Generate plots:
G.m1 <- create.scatterplot(Dm.OK2,m1,ylab='Mean otolith \n growth rate (µm/d)',xlab='',Month=TRUE)
G.m2 <- create.boxplot(Dm.OK2,m1a)
G.m3 <- create.scatterplot(Dv.OK3,m2pR,ylab='SD otolith \n growth rate (µm/d)',xlab='',Month=TRUE,Site=TRUE,Resid=TRUE)
G.m4 <- create.boxplot(Dv.OK3,m2a)
G.m5 <- create.scatterplot(Dmax.OK2,m3,ylab='Max otolith \n growth rate (µm/d)',xlab='',Month=TRUE)
G.m6 <- create.boxplot(Dmax.OK2,m3a)
G.m7 <- create.scatterplot(Dmin.OK2,m4.1,ylab='Min otolith \n growth rate (µm/d)',xlab='Group size',Month=TRUE)
G.m8 <- create.boxplot(Dmin.OK2,m4a)

quartz(width=2.4)
ggarrange(G.m1,G.m3,G.m5, ncol = 1)

#----------------------------------------------------------
# Analysis of plankton patterns over space and time:
P = read.csv('Plankton_v2.csv')
P$log.c = log10(P$Copepods)


ggplot(data=P,aes(x=Site,y=log.c))+
  geom_jitter(aes(x=factor(Site),y=log.c,color=Month),width=0.1)+
  scale_color_manual(values=c('blue','red'))+
  theme_bw()+
  theme(legend.position='none')

# So plankton varied a lot among sites and months, but for a given month the
# variability was similar across sites. See above for test of effect on growth - there is none.
#-----------------------------------------------------------

############################################################
# Power analysis

# Mean
Actual = abs(coef(m1)[2])
SD.effect = sqrt(vcov(m1)[2,2])
Effects = seq(0,1,length.out=100)
g1 = power.analysis(S=SD.effect,Effects=Effects,Actual=Actual)

# SD
Actual = coef(m2)[2]
SD.effect = sqrt(vcov(m2)[2,2])
Effects = seq(0,10,length.out=100)
g2 = power.analysis(S=SD.effect,Effects=Effects,Actual=Actual)


# Max
Actual = abs(coef(m3)[2])
SD.effect = sqrt(vcov(m3)[2,2])
Effects = seq(0,1,length.out=100)
g3 = power.analysis(S=SD.effect,Effects=Effects,Actual=Actual)

# Plot results:
quartz(width=2.4)
ggarrange(g1,g2,g3, ncol = 1)
############################################################


############################################################
# Some helper functions (run these first!):


#-------------------------------------------------------------------
# POWER ANALYSIS
power.analysis <- function(S,Effects,Actual){
  # loop over effect sizes and sample sizes
  # S = observation standard deviation of effect
  # Effects = range of effect sizes to test
  Ns = seq(10,100,by=1)
  Power <- rep(NA,length(Effects)*length(Ns))
  Effects.df = Power
  Ns.df = Power
  for (e in 1:length(Effects)){
    for (n in 1:length(Ns)){
      
      X <- pwr.t.test(n =  Ns[n], d = Effects[e]/S, sig.level = 0.05, power = NULL, type = c("two.sample"))
      Index <- (e-1)*length(Effects) + n
      Power[Index] = X$power
      Effects.df[Index] = Effects[e]
      Ns.df[Index] = Ns[n]
      
    }}
  
  Power.df  = data.frame(Effect=Effects.df,N=Ns.df,Power=Power)
  
  # Plot results:
  g<-ggplot(data=Power.df,mapping=aes(x=N,y=Effect,z=Power))+
    geom_contour(aes(x=N,y=Effect,color=stat(level)),breaks=seq(0.2,0.8,by=0.2))+
    theme_bw()+
    geom_hline(yintercept=Actual,color='black')+
    xlab('Sample size')+
    ylab('Effect size')
  
  g=direct.label(g, list("bottom.pieces", color='black'))
}
#-----------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------
create.boxplot<-function(D,m,ylab=NULL,xlab=NULL){
G.obj<-ggplot(D,aes(x=as.factor(isGroup),y=x))+
  geom_boxplot(aes(x=as.factor(isGroup),y=x))+
  geom_jitter(size=1,shape=21,col='black',height=0,width=0.2)+ 
  ylab(ylab)+xlab("")+
scale_x_discrete(labels=c('Solitary','Group'))+
  theme(axis.title.y=element_text(size=10))+
  theme_bw()
G.obj
}
#-----------------------------------------------------------------------------------

create.scatterplot<-function(D,m,ylab=NULL,xlab=NULL,Month=FALSE,Site=FALSE,Resid=FALSE){
  # Helper function to create a scatterplot. Only works in context of the lm objects used in this code
  if (Resid){D$x = D$Resid} # if plotting residuals instead of raw data
  
  if (Month & !Site){
    Months = c('7','8')
    # intervals:
    D$fitted <- NA
    D$lm.lower <- NA
    D$lm.upper <- NA
    
    for (i in 1:2){
      OKm = D$Month==Months[i]
      n = length(D$x[OKm])
      xm = mean(D$Grpsize.inv[OKm])
      Sxy.mean = sqrt(sum(resid(m)[OKm]^2)/(n-2)) # SD of residuals
      SEline.mean = Sxy.mean * sqrt(1/n + ((D$Grpsize.inv[OKm] - xm)^2)/((n-1)*var(D$Grpsize.inv[OKm]))) # equation for prediction SE of line
      
      if (class(m)[1]=='lmerMod'){
        c = fixef(m)
      }else{
        c = coef(m)}
      if (i == 2){
        isMonth = 1}else{isMonth = 0}
      D$fitted[OKm] <- ((1/D$Grpsize[OKm])*c[2]+c[1]+isMonth*c[3])
      D$lm.lower[OKm] <- (D$fitted[OKm] - 1.96*SEline.mean)-1
      D$lm.upper[OKm] <- (D$fitted[OKm] + 1.96*SEline.mean)-1
      D$fitted[OKm] <- (D$fitted[OKm])-1
    }
    
    D1 = D[OKm,]
    D2 = D[!OKm,]
    
    # plot:
    G.obj<-ggplot(D1,aes(x=Grpsize,y=x))+
      geom_ribbon(data=D1,aes(x=Grpsize,ymin=lm.lower,ymax=lm.upper),alpha=0.5,fill="grey80")+
      geom_ribbon(data=D2,aes(x=Grpsize,ymin=lm.lower,ymax=lm.upper),alpha=0.5,fill="grey40")+
      geom_line(data=D1,aes(x=Grpsize,y=fitted),size=0.5)+
      geom_line(data=D2,aes(x=Grpsize,y=fitted),size=0.5,linetype=2)+
      geom_jitter(data=D1,size=1,shape=21,height=0,width=0.1)+  
      geom_jitter(data=D2,size=1,shape=24,height=0,width=0.1)+
      scale_x_continuous(breaks=c(2,6,10,14))+
      ylab(ylab)+xlab(xlab)+
      theme_bw()+
      theme(axis.title.y=element_text(size=10))+
      guides(color="none")
    G.obj
  }
  else if (Site & !Month){
    Sites = c('1','2','3')
    # intervals:
    D$fitted <- NA
    D$lm.lower <- NA
    D$lm.upper <- NA
    
    for (i in 1:3){
      OKm = D$Site==Sites[i]
      n = length(D$x[OKm])
      xm = mean(D$Grpsize.inv[OKm])
      Sxy.mean = sqrt(sum(resid(m)[OKm]^2)/(n-2)) # SD of residuals
      SEline.mean = Sxy.mean * sqrt(1/n + ((D$Grpsize.inv[OKm] - xm)^2)/((n-1)*var(D$Grpsize.inv[OKm]))) # equation for prediction SE of line
      
      if (class(m)[1]=='lmerMod'){
        c = fixef(m)
      }else{
        c = coef(m)}
      if (i == 2){isSite1 = 1}else{isSite1 = 0}
      if (i == 3){isSite2 = 1}else{isSite2 = 0}
      D$fitted[OKm] <- ((1/D$Grpsize[OKm])*c[2]+c[1]+isSite1*c[3]+isSite2*c[4])
      D$lm.lower[OKm] <- (D$fitted[OKm] - 1.96*SEline.mean)-1
      D$lm.upper[OKm] <- (D$fitted[OKm] + 1.96*SEline.mean)-1
      D$fitted[OKm] <- (D$fitted[OKm])-1
    }
    
    D1 = D[D$Site=='1',]
    D2 = D[D$Site=='2',]
    D3 = D[D$Site=='3',]
    
    # plot:
    G.obj<-ggplot(D1,aes(x=Grpsize,y=x))+
      geom_ribbon(data=D1,aes(x=Grpsize,ymin=lm.lower,ymax=lm.upper),alpha=0.5,fill="grey80")+
      geom_ribbon(data=D2,aes(x=Grpsize,ymin=lm.lower,ymax=lm.upper),alpha=0.5,fill="grey60")+
      geom_ribbon(data=D3,aes(x=Grpsize,ymin=lm.lower,ymax=lm.upper),alpha=0.5,fill="grey40")+
      geom_line(data=D1,aes(x=Grpsize,y=fitted),size=0.5,linetype='solid')+
      geom_line(data=D2,aes(x=Grpsize,y=fitted),size=0.5,linetype='dashed')+
      geom_line(data=D3,aes(x=Grpsize,y=fitted),size=0.5,linetype='dotted')+
      geom_jitter(data=D1,size=1,shape=21,height=0,width=0.1)+  
      geom_jitter(data=D2,size=1,shape=24,height=0,width=0.1)+
      geom_jitter(data=D3,size=1,shape=23,height=0,width=0.1)+
      scale_x_continuous(breaks=c(2,6,10,14))+
      ylab(ylab)+xlab(xlab)+
      theme_bw()+
      theme(axis.title.y=element_text(size=10))+
      guides(color="none")
    G.obj
  } 
  else if (Site & Month){
    Sites = c('1','2','3')
    Months = c('7','8')
    # intervals:
    D$fitted <- NA
    D$lm.lower <- NA
    D$lm.upper <- NA
    
    for (i in 1:3){
      OKm = D$Site==Sites[i]
      n = length(D$x[OKm])
      xm = mean(D$Grpsize.inv[OKm])
      Sxy.mean = sqrt(sum(resid(m)[OKm]^2)/(n-2)) # SD of residuals
      SEline.mean = Sxy.mean * sqrt(1/n + ((D$Grpsize.inv[OKm] - xm)^2)/((n-1)*var(D$Grpsize.inv[OKm]))) # equation for prediction SE of line
      
      if (class(m)[1]=='lmerMod'){
        c = fixef(m)
      }else{
        c = coef(m)}
      if (i == 2){isSite1 = 1}else{isSite1 = 0}
      if (i == 3){isSite2 = 1}else{isSite2 = 0}
      D$fitted[OKm] <- ((1/D$Grpsize[OKm])*c[2]+c[1]+isSite1*c[3]+isSite2*c[4]+0.5*c[4])
      D$lm.lower[OKm] <- (D$fitted[OKm] - 1.96*SEline.mean)-1
      D$lm.upper[OKm] <- (D$fitted[OKm] + 1.96*SEline.mean)-1
      D$fitted[OKm] <- (D$fitted[OKm])-1
    }
    
    D1 = D[D$Site=='1',]
    D2 = D[D$Site=='2',]
    D3 = D[D$Site=='3',]
    D1a = D[D$Site=='1'&D$Month=='7',]
    D2a = D[D$Site=='2'&D$Month=='7',]
    D3a = D[D$Site=='3'&D$Month=='7',]
    D1b = D[D$Site=='1'&D$Month=='8',]
    D2b = D[D$Site=='2'&D$Month=='8',]
    D3b = D[D$Site=='3'&D$Month=='8',]
    
    # plot:
    G.obj<-ggplot(D1,aes(x=Grpsize,y=x))+
      geom_ribbon(data=D1,aes(x=Grpsize,ymin=lm.lower,ymax=lm.upper),alpha=0.5,fill="grey80")+
      geom_ribbon(data=D2,aes(x=Grpsize,ymin=lm.lower,ymax=lm.upper),alpha=0.5,fill="grey60")+
      geom_ribbon(data=D3,aes(x=Grpsize,ymin=lm.lower,ymax=lm.upper),alpha=0.5,fill="grey40")+
      geom_line(data=D1,aes(x=Grpsize,y=fitted),size=0.5,linetype='solid')+
      geom_line(data=D2,aes(x=Grpsize,y=fitted),size=0.5,linetype='dashed')+
      geom_line(data=D3,aes(x=Grpsize,y=fitted),size=0.5,linetype='dotted')+
      geom_jitter(data=D1a,size=1,shape=21,height=0,width=0.1,color='red')+  
      geom_jitter(data=D2a,size=1,shape=24,height=0,width=0.1,color='red')+
      geom_jitter(data=D3a,size=1,shape=23,height=0,width=0.1,color='red')+
      geom_jitter(data=D1b,size=1,shape=21,height=0,width=0.1,color='blue')+  
      geom_jitter(data=D2b,size=1,shape=24,height=0,width=0.1,color='blue')+
      geom_jitter(data=D3b,size=1,shape=23,height=0,width=0.1,color='blue')+
      scale_x_continuous(breaks=c(2,6,10,14))+
      ylab(ylab)+xlab(xlab)+
      theme_bw()+
      theme(axis.title.y=element_text(size=10))+
      guides(color="none")
    G.obj
  }
  else{  # no effect of month or site
    n = length(D$x)
    xm = mean(D$Grpsize.inv)
    Sxy.mean = sqrt(sum(resid(m)^2)/(n-2)) # SD of residuals
    SEline.mean = Sxy.mean * sqrt(1/n + ((D$Grpsize.inv - xm)^2)/((n-1)*var(D$Grpsize.inv))) # equation for prediction SE of line
    # intervals:
    D$fitted <- NA
    if (class(m)[1]=='lmerMod'){
      c = fixef(m)
    }else{
      c = coef(m)}
    D$fitted <- ((1/D$Grpsize)*c[2]+c[1])
    D$lm.lower <- (D$fitted - 1.96*SEline.mean)-1
    D$lm.upper <- (D$fitted + 1.96*SEline.mean)-1
    D$fitted <- (D$fitted)-1
    # plot:
    G.obj<-ggplot(D,aes(x=Grpsize,y=x))+
      geom_ribbon(aes(x=Grpsize,ymin=lm.lower,ymax=lm.upper),fill="grey80")+
      geom_line(aes(x=Grpsize,y=fitted),size=0.5)+
      geom_jitter(size=1,col='black',shape=21,height=0,width=0.1)+ 
      scale_x_continuous(breaks=c(2,6,10,14))+
      ylab(ylab)+xlab(xlab)+
      theme_bw()+
      theme(axis.title.y=element_text(size=10))
    G.obj
  } # end if Month
}
#-----------------------------------------------------------------------------------

