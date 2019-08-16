# script for analyzing bluehead wrasse gut contents data 
# for Bogdan et al. ms (Ethology)

# NOTE: various functions appear at the end of this file; they must be added to the workspace first.

require(ggplot2)
require(MuMIn)
require(gridExtra)
require(egg)
require(pwr)
require(directlabels)


# Read in data
D = read.csv('Guts+Otos.csv')

D$Grp = D$G.or.S # new category - grouped or solitary
D = subset(D,D$otolith==TRUE) # subset out the fish that did not have otos read
D = subset(D,D$PS>0) # subset out fish that are < 1 day old
D$Grp = factor(D$Grp)


D$grpsize[D$grpsize==0] = 1 # adjust the groupsize category so solitaires have a size of 1

# Some summaries of the dataset:
D$totalprey = rowSums(D[,2:14])
D$other = rowSums(D[,9:14]) # low frequency items
D$copepods = rowSums(D[,2:6])
D$Breadth = 1/rowSums((D[,2:14]/D$totalprey)^2)
D$Breadth2 = 1/rowSums((D[,c(7:14,44)]/D$totalprey)^2)

# Note: diet is 65% copepods

# Need a variable for the unique group number (each unique group regardless of date or site)
# This is pretty arbitrary but it gives a unique number to each group:
D$Grpnum2 = as.factor(D$month*0.1 + D$site*0.01 + D$day*0.001 + D$grpnum)
D$cyclo = rowSums(D[,c(2:4,6)])

D$month = factor(D$month)
D$site = factor(D$site)

#---------------------------------------------------------------------
# Summarize overall diet composition:
TotalDiet1 = colSums(D[,c(2:14)]) # all categories
TotalDiet2 = colSums(D[,c(2:8,44)]) # pooling
TotalDiet3 = colSums(D[,c(5,7:8,44,49)]) # pooling more
TotalDiet.G = colSums(D[D$Grp=='G',c(5,7:8,44,49)])
TotalDiet.S = colSums(D[D$Grp=='S',c(5,7:8,44,49)])
TotalDiet.Glg = colSums(D[D$Grp=='G'&D$TL>14.4,c(5,7:8,44,49)])
TotalDiet.Gsm = colSums(D[D$Grp=='G'&D$TL<=14.4,c(5,7:8,44,49)])
TotalDiet.Slg = colSums(D[D$Grp=='S'&D$TL>14.4,c(5,7:8,44,49)])
TotalDiet.Ssm = colSums(D[D$Grp=='S'&D$TL<=14.4,c(5,7:8,44,49)])

# Make some pie charts:
pie(TotalDiet1)
pie(TotalDiet2)
pie(TotalDiet3)
pie(TotalDiet.G)
pie(TotalDiet.S)
pie(TotalDiet.Glg)
pie(TotalDiet.Gsm)
pie(TotalDiet.Slg)
pie(TotalDiet.Ssm)
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Test for consistency in diet:
m = chisq.test(rbind(TotalDiet.G,TotalDiet.S),simulate.p.value=FALSE)
m = chisq.test(rbind(TotalDiet.Gsm,TotalDiet.Glg),simulate.p.value=FALSE)
m = chisq.test(rbind(TotalDiet.Ssm,TotalDiet.Slg),simulate.p.value=FALSE)
m = chisq.test(rbind(TotalDiet.G,TotalDiet.Ssm),simulate.p.value=FALSE)
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Exploratory plots:
plot(D$Grp,D$Harpac)
plot(D$Grp,D$totalprey)
plot(D$Grp,D$Breadth)
plot(D$Grp,jitter(D$copepods))

plot(D$Grp,D$smCyclo)
plot(D$Grp,D$medCyclo)
plot(D$Grp,D$lgCyclo)
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Following otolith analysis logic, group by site & compare mean, sd, etc.

# Number of observations within each group (for later filtering)
D$counter= as.numeric(D$ps.rate > 0 )
Dl = aggregate(D$counter,by=list(Site=D$site,Month=D$month,
                                 Date=D$date,Grpnum=D$grpnum,
                                 Grpsize=D$grpsize,Plankton=D$Plankton),FUN='sum',na.rm=TRUE)
OK2 = Dl$x >= 2 # only groups with at least 2 observations
OK3 = Dl$x >= 3 # only groups with at least 2 observations

# Can plot Dl$x vs. Dm$x, etc. to see that there somewhat of an effect of # observations on the statistics
# but sample size is evenly distributed across groups.

# Mean:
Dm = aggregate(D$copepods,by=list(Site=D$site,Month=D$month,
                                 Date=D$date,Grpnum=D$grpnum,
                                 Grpsize=D$grpsize,Plankton=D$Plankton),
                              FUN='mean',na.rm=TRUE)
Dm$isGroup = Dm$Grpsize>1
Dm$Grpsize.inv = 1/Dm$Grpsize # we model the effects of group size as 1/grpsize because one expects competition to depend on the proportion of a group 1 fish represents.

# Variance
Dv = aggregate(D$copepods,by=list(Site=D$site,Month=D$month,
                                  Date=D$date,Grpnum=D$grpnum,
                                  Grpsize=D$grpsize,Plankton=D$Plankton),
               FUN='var',na.rm=TRUE)
Dv$isGroup = Dv$Grpsize>1
Dv$Grpsize.inv = 1/Dv$Grpsize


# Maximum
Dmax = aggregate(D$copepods,by=list(Site=D$site,Month=D$month,
                                  Date=D$date,Grpnum=D$grpnum,
                                  Grpsize=D$grpsize,Plankton=D$Plankton),
               FUN='max',na.rm=TRUE)
Dmax$isGroup = Dmax$Grpsize>1
Dmax$Grpsize.inv = 1/Dmax$Grpsize


##############################################################
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
Dv.OK3$x <- sqrt(Dv.OK3$x) # variance --> SD
##############################################################

##############################################################
# Linear models that include ambient plankton density as a factor:

# Just plankton:
m0 = lm(log10(x+1)~Plankton,data=Dm[OK2,])

# Mean
m1p = lm(log10(x+1)~Grpsize.inv +Plankton + Site + Month,data=Dm.OK2)# full model
m1p = lm(log10(x+1)~Grpsize.inv + Plankton + Month,data=Dm.OK2)
m1p = lm(log10(x+1)~Grpsize.inv + Month,data=Dm[OK2,]) # reduced model
summary(m1p)

m1p.1 = lm(log10(x+1)~isGroup+Site + Plankton + Month,data=Dm.OK2) # full model
m1p.1 = lm(log10(x+1)~isGroup + Plankton+Month,data=Dm.OK2)
m1p.1 = lm(log10(x+1)~isGroup + Month,data=Dm.OK2) #reduced
summary(m1p.1)


# Standard deviation
m2p=(lm(log10(sqrt(x)+1)~Grpsize.inv+Plankton+Site+Month,data=Dv.OK3)) # full model
m2p=(lm(log10(sqrt(x)+1)~Grpsize.inv+Plankton+Site,data=Dv.OK3))
summary(m2p)
# create residuals vs plankton for plotting:
m2p0=(lm(log10(sqrt(x)+1)~Plankton,data=Dv.OK3))
Dv.OK3$Resid = resid(m2p0)
m2pR=(lm(Resid~Grpsize.inv+Site,data=Dv.OK3))


# Maximum
m3p = lm(log10(x+1)~ Grpsize.inv + Plankton+Site+Month, data=Dmax.OK2) # full
m3p = lm(log10(x+1)~ Grpsize.inv +Site+Month, data=Dmax.OK2)
m3p = lm(log10(x+1)~ Grpsize.inv+Month, data=Dmax.OK2) # reduced
summary(m3p)

##############################################################


# Generate plots:
G.m1 <- create.scatterplot(Dm.OK2,m1p,ylab='Mean no. copepods',xlab='',Month=TRUE,Site=FALSE)
G.m2 <- create.boxplot(Dm.OK2,m1p.1,Month=TRUE,Site=FALSE)
G.m3 <- create.scatterplot(Dv.OK3,m2pR,ylab='Residual SD no. copepods',xlab='',Month=FALSE,Site=TRUE,Resid=TRUE)
G.m4 <- create.boxplot(Dv.OK3,m2p.1R,Month=TRUE,Site=FALSE,Resid=TRUE)
G.m5 <- create.scatterplot(Dmax.OK2,m3p,ylab='Max no. copepods',xlab='',Month=TRUE,Site=FALSE)
G.m6 <- create.boxplot(Dmax.OK2,m3p.1,Month=TRUE,Site=FALSE)


quartz(width=2.4)
ggarrange(G.m1,G.m3,G.m5, ncol = 1)
##############################################################

##############################################################
# Power analysis

# Mean
Actual = coef(m1p)[2]
SD.effect = sqrt(vcov(m1p)[2,2])
Effects = seq(0,1,length.out=100)
g1 = power.analysis(S=SD.effect,Effects=Effects,Actual=Actual)

# SD
Actual = coef(m2p)[2]
SD.effect = sqrt(vcov(m2p)[2,2])
Effects = seq(0,1,length.out=100)
g2 = power.analysis(S=SD.effect,Effects=Effects,Actual=Actual)


# Max
Actual = coef(m3p)[2]
SD.effect = sqrt(vcov(m3p)[2,2])
Effects = seq(0,1,length.out=100)
g3 = power.analysis(S=SD.effect,Effects=Effects,Actual=Actual)

# plot results
quartz(width=2.4)
ggarrange(g1,g2,g3, ncol = 1)

##############################################################
# HELPER FUNCTIONS (run these first!):

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
create.scatterplot<-function(D,m,ylab=NULL,xlab=NULL,Month=FALSE,Site=FALSE,Resid=FALSE){
  # Helper function to create a scatterplot. Only works in context of the lm objects used in this code
  if (Resid){D$x = D$Resid} # if plotting residuals instead of raw data
  
  if (Month){
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
      D$lm.lower[OKm] <- 10^(D$fitted[OKm] - 1.96*SEline.mean)-1
      D$lm.upper[OKm] <- 10^(D$fitted[OKm] + 1.96*SEline.mean)-1
      D$fitted[OKm] <- 10^(D$fitted[OKm])-1
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
  else if (Site){
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
      D$lm.lower[OKm] <- 10^(D$fitted[OKm] - 1.96*SEline.mean)-1
      D$lm.upper[OKm] <- 10^(D$fitted[OKm] + 1.96*SEline.mean)-1
      D$fitted[OKm] <- 10^(D$fitted[OKm])-1
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
    D$lm.lower <- 10^(D$fitted - 1.96*SEline.mean)-1
    D$lm.upper <- 10^(D$fitted + 1.96*SEline.mean)-1
    D$fitted <- 10^(D$fitted)-1
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

#-----------------------------------------------------------------------------------
create.boxplot<-function(D,m,ylab=NULL,xlab=NULL,Month=FALSE,Site=FALSE,Resid=FALSE){
  
  if (Resid){D$x = D$Resid} # if plotting residuals instead of raw data
  
  
  if (Month){
  G.obj<-ggplot(D,aes(x=as.factor(isGroup),y=x))+
    geom_boxplot(aes(x=as.factor(isGroup),y=x),outlier.shape=NA)+
    geom_jitter(size=1,aes(shape=Month),fill=NA,height=0,width=0.2)+     
    scale_shape_manual(values=c(24,21))+
    ylab(ylab)+xlab("")+
    scale_x_discrete(labels=c('Solitary','Group'))+
    theme(axis.title.y=element_text(size=10))+
    theme_bw()+
    guides(shape="none")
  }
  else if (Site){
    G.obj<-ggplot(D,aes(x=as.factor(isGroup),y=x))+
      geom_boxplot(aes(x=as.factor(isGroup),y=x),outlier.shape=NA)+
      geom_jitter(size=1,aes(shape=Site),fill=NA,height=0,width=0.2)+     
      scale_shape_manual(values=c(24,21))+
      ylab(ylab)+xlab("")+
      scale_x_discrete(labels=c('Solitary','Group'))+
      theme(axis.title.y=element_text(size=10))+
      theme_bw()+
      guides(shape="none")
  }  
    
  
  
  G.obj
}

