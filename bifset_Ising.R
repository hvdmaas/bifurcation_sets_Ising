# R Code for 
# Analytical Bifurcation Analysis of Mean-Field Ising Models Reveals Connectivity as a
# Risk Factor for Psychopathology

# Han L. J. van der Maas & Lourens Waldorp
# Uva
# May 2025

###### MF Ising -11
##Figure 2
# beta tau
pdf('figure2.pdf',h=8,w=10)#,paper='a4r')
par(mar=c(4, 3, 2, 2))
layout(matrix(c(1,2,3,4),2,2,byrow=T))

MFtau_critical1 <- function(beta,sigma=1) # equation 4, solution 1
(1/beta)*atanh(sqrt(1-1/(sigma*beta)))-sigma*sqrt(1-1/(sigma*beta))

MFtau_critical2 <- function(beta,sigma=1) # equation 4, solution 2
  (-1/beta)*atanh(sqrt(1-1/(sigma*beta)))+sigma*sqrt(1-1/(sigma*beta))

sigma=1
# panel 1: bifurcation set beta tau
curve(MFtau_critical1,0,8,ylim=c(-1.3,1.3),ylab=expression(tau),xlab=expression(beta),bty='n', 
      main=expression(paste('Bifurcation set for ',beta,',',tau,' (',sigma,' = 1)')),cex.main=1.6,cex.lab=1.8,cex.axis=1.2)
curve(MFtau_critical2,0,8,add=T)

# add two lines to compare with hysteresis plots
beta=4
abline(h= MFtau_critical1(beta),col='blue')
abline(h= MFtau_critical2(beta),col='blue')

# panel 2: bifurcation set beta tau, wider interval of beta

sigma=1
curve(MFtau_critical1,1,200,ylim=c(-1,1),ylab=expression(tau),xlab=expression(beta),bty='n', 
      main=expression(paste('Bifurcation set for ',beta,',',tau,' (',sigma,' = 1)')),cex.main=1.6,cex.lab=1.8,cex.axis=1.2)
curve(MFtau_critical2,1,200,add=T)

MFIsing <- function(mu) {  tanh(beta * (tau + mu * d * sigma)) - mu } # MF equation, equation 2
d=1;sigma=1
taus=seq(-1,1,length=200)
betas=c(1,4,7,100)
dat=matrix(0,length(taus)*length(betas),5)
i=0
for(beta in betas)
  for(tau in taus)
  {
    i=i+1
    sol <- uniroot.all(MFIsing, c(-1.1, 1.1))
    if(length(sol)==1) sol=c(sol, NA, NA)
    dat[i,]=c(beta,tau,sol)
  }


# hysteresis plot panel 3
matplot(dat[,2],dat[,-1:-2],type='p',cex.main=1.6,cex.lab=1.8,cex.axis=1.2,
        pch=1,col=c(1,2,1),bty='n',xlab=expression(tau),ylab=expression(mu),
        cex=.2,main=expression(paste('Hysteresis for ',beta, '= 1,4,7,100')))

# check jump positions
abline(v= MFtau_critical1(1)); abline(v= MFtau_critical2(1))
abline(v= MFtau_critical1(4),col='blue'); abline(v= MFtau_critical2(4),col='blue')
abline(v= MFtau_critical1(7)); abline(v= MFtau_critical2(7))
abline(v= MFtau_critical1(100)); abline(v= MFtau_critical2(100))


# check with real Ising simulation:
# panel 4
resample <- function(x, ...) x[sample.int(length(x), ...)] # as sample(5,1) acts odd
# parmaters
n=200
conn=matrix(1,n,n)
diag(conn)=0
sigma=1
betas=c(4,7,100)
conn=conn*sigma/n
spins=c(-1,1)

iterations=2000
iter=0
for(beta in betas)
{
  iter=iter+1  
  mm=rep(0,iterations)
  tau=rep(-1.5,iterations)
  d_tau=.01;dir=1
  x=sample(spins,n,TRUE)
  for ( i in 1:iterations) # vary tau
  {
    if(i>1)
    {
      dir[i]=dir[i-1]
      if(tau[i-1]  > 2.0) dir[i]=-1 
      if(tau[i-1] < -2.0) dir[i]=1 
      tau[i]=tau[i-1]+dir[i]*d_tau
    }
    for(ii in 1:1000) # do 1000 metropolis steps
    {
      j=sample(n,1)
      spin=x[j]
      spin_n=resample(spins[-which(abs(spins-spin)<.000001)],1)
      m=sum(x*conn[j,])
      dH=-(spin-spin_n)*m - tau[i]*(spin-spin_n);
      p=min(1,exp(beta*dH))  # metropolis
      #p=1/(1+exp(-beta*dH)) # glauber
      if(runif(1) < p) x[j] = spin_n # update state
    }
    mm[i]=sum(x)
  }
  
  tau=round(100000*tau,0)/100000
  t=tapply(mm/n,list(tau,dir),mean,na.rm=T)
  t=cbind(t,3)
  if(iter==1) data=t else data=rbind(data,t)
}

# make plot
datal=cbind(as.numeric(rownames(data)),data)
colnames(datal)[1]='tau'

datal=reshape(as.data.frame(datal),
              varying = c("-1", "1"),  # Columns to pivot
              v.names = c("magn."),                   # Name for values
              timevar = "Dir",                    # Name for the time variable
              times = c("down", "up"),         # New variable levels
              direction = "long")


plot(datal[,1],datal[,4],type='b',bty='n',xlab= expression(tau),lwd=2,ylab=expression(mu),xlim=c(-1,1),col=as.numeric(factor(datal[,2])),pch=1,cex=.2,cex.main=1.6,cex.lab=1.8,cex.axis=1.2,
     main=expression(paste('Hysteresis for ',beta,' = 4,7,100')))

# add predictions
beta=4
abline(v= MFtau_critical1(beta),col='blue'); abline(v= MFtau_critical2(beta),col='blue')
beta=7
abline(v= MFtau_critical1(beta),col='blue'); abline(v= MFtau_critical2(beta),col='blue')
beta=100
abline(v= MFtau_critical1(beta),col='blue'); abline(v= MFtau_critical2(beta),col='blue')



V_Ising <- function(mu) {
  (0.5 * mu^2 - (1 / (beta * sigma)) * log(cosh(beta * (tau + sigma * mu))))
}
d=.06
par(mar=c(0,0,0,0))

beta=1;tau=0
x=.05;y=.75;
par(fig = c(x, x+d, y,y+d), new = TRUE)
curve(V_Ising,-2,2,bty='n',xlab='',ylab='',xaxt = "n",yaxt = "n")

beta=4;tau=0
x=.2;y=.74;
par(fig = c(x, x+d, y,y+d), new = TRUE)
curve(V_Ising,-2,2,bty='n',xlab='',ylab='',xaxt = "n",yaxt = "n")

beta=4;tau=.7
x=.2;y=.85;
par(fig = c(x, x+d, y,y+d), new = TRUE)
curve(V_Ising,-2,2,bty='n',xlab='',ylab='',xaxt = "n",yaxt = "n")

beta=4;tau=-.7
x=.2;y=.63;
par(fig = c(x, x+d, y,y+d), new = TRUE)
curve(V_Ising,-2,2,bty='n',xlab='',ylab='',xaxt = "n",yaxt = "n")

par(fig = c(0, 1, 0, 1), new = FALSE)
par(mar=c(4, 4, 3, 2))



dev.off()

## Figure 3
# sigma tau
pdf('figure3.pdf',h=4,w=10)#,paper='a4r')

layout(matrix(1:2,1,2))
  
MFtau_critical1s <- function(sigma,beta=1)  # same MFtau_critical1 but with different order of arguments
  (1/beta)*atanh(sqrt(1-1/(sigma*beta)))-sigma*sqrt(1-1/(sigma*beta))

MFtau_critical2s <- function(sigma,beta=1) 
  (-1/beta)*atanh(sqrt(1-1/(sigma*beta)))+sigma*sqrt(1-1/(sigma*beta))

# panel 1
beta=1
curve(MFtau_critical1s,1,20,ylim=c(-20,20),ylab=expression(tau),xlab=expression(sigma),bty='n', 
      main=expression(paste('Bifurcation set for ',sigma,',',tau,' (',beta,' = 1)')),cex.main=1.6,cex.lab=1.8,cex.axis=1.2)
curve(MFtau_critical2s,1,20,add=T)

# check hysteresis at sigma = 10
sigma=10
abline(h= MFtau_critical1s(sigma),col='blue')
abline(h= MFtau_critical2s(sigma),col='blue')

# hysteresis plot using the MF equation
g <- function(mu) {  tanh(beta * (tau + mu * d * sigma)) - mu }
d=1;beta=1
taus=seq(-10,10,length=200)
sigmas=seq(1,10,by=3)
dat=matrix(0,length(taus)*length(sigmas),5)
i=0
for(sigma in sigmas)
  for(tau in taus)
  {
    i=i+1
    sol <- uniroot.all(g, c(-10, 10),n=2000)
    if(length(sol)==1) sol=c(sol, NA, NA)
    dat[i,]=c(sigma,tau,sol)
  }
matplot(dat[,2],dat[,-1:-2],type='p',pch=1,col=c(1,2,1),bty='n',cex.main=1.6,cex.lab=1.8,cex.axis=1.2,
        xlab=expression(tau),ylab=expression(mu),
        cex=.4,main = expression(paste("Hysteresis for ", sigma, " = 1, 4, 7, 10")))


abline(v= MFtau_critical1s(4)); abline(v= MFtau_critical2s(4))
abline(v= MFtau_critical1s(7)); abline(v= MFtau_critical2s(7))
abline(v= MFtau_critical1s(10),col='blue'); abline(v= MFtau_critical2s(10),col='blue')
dev.off()
###### MF Ising 01

#figure 4
pdf('figure4.pdf',h=8,w=10)#,paper='a4r')


layout(matrix(1:4,2,2))
sigma=1

# equation 8: both solutions
MF01tau_critical1=function(b,s=sigma) 
{q=sqrt(1-4/(s*b))
(-s/2)*(1+q)+(1/b)*log((1+q)/(1-q))
}

MF01tau_critical2=function(b,s=sigma) 
{q=sqrt(1-4/(s*b))
(-s/2)*(1-q)+(1/b)*log((1-q)/(1+q))
}

#panel 1 beta
curve(MF01tau_critical1,4,15,ylim=c(-1,0),ylab=expression(tau),xlab=expression(beta),bty='n',cex.main=1.6,cex.lab=1.8,cex.axis=1.2,
      main=expression(paste('Bifurcation set for ',beta,',',tau,' (',sigma,' = 1)')))
curve(MF01tau_critical2,4,15,add=T)

# equation 5
MF01 <-function(p) 1/(1+exp(-beta*(d*sigma*p+tau)))-p
d=1

# panel 3: hysteresis
taus=seq(-1.5,.2,length=200)
betas=c(4,7,100)
dat=matrix(0,length(taus)*length(betas),5)
i=0
for(beta in betas)
 for(tau in taus)
{
  i=i+1
 sol <- uniroot.all(MF01, c(-3, 3),n=2000)
  if(length(sol)==1) sol=c(sol, NA, NA)
  dat[i,]=c(tau,beta,sol)
}
matplot(dat[,1],dat[,-1:-2],type='p',
        pch=1,col=c(1,2,1),bty='n',xlab=expression(tau),ylab='p',cex.main=1.6,cex.lab=1.8,cex.axis=1.2,
        cex=.4,main=expression(paste('Hysteresis for ',beta, '= 4,7,10,100')))

# checks
beta=4
abline(v= MF01tau_critical1(b=beta))
abline(v= MF01tau_critical2(b=beta))
beta=7
abline(v= MF01tau_critical1(b=beta))
abline(v= MF01tau_critical2(b=beta))
beta=100
abline(v= MF01tau_critical1(b=beta))
abline(v= MF01tau_critical2(b=beta))

#Panel 2
# bifurcation sigma, tau
MF01tau_critical1=function(s,b=1) 
{q=sqrt(1-4/(s*b))
(-s/2)*(1+q)+(1/b)*log((1+q)/(1-q))
}

MF01tau_critical2=function(s,b=1) 
{q=sqrt(1-4/(s*b))
(-s/2)*(1-q)+(1/b)*log((1-q)/(1+q))
}


curve(MF01tau_critical1,4,15,ylab=expression(tau),xlab=expression(sigma),bty='n', ylim=c(-15,0),
      main=expression(paste('Bifurcation set for ',sigma,',',tau,' (',beta,' = 1)')),cex.main=1.6,cex.lab=1.8,cex.axis=1.2)
curve(MF01tau_critical2,4,15,add=T)

# check hysteresis
d=1
beta=1
taus=seq(-12,0,length=200)
sigmas=c(4,7,10,15)
dat=matrix(0,length(taus)*length(sigmas),5)
i=0
for(sigma in sigmas)
  for(tau in taus)
  {
    i=i+1
    sol <- uniroot.all(MF01, c(-2, 2))
    if(length(sol)==1) sol=c(sol, NA, NA)
    dat[i,]=c(tau,sigma,sol)
  }
matplot(dat[,1],dat[,-1:-2],type='p',
        pch=1,col=c(1,2,1),bty='n',xlab=expression(tau),ylab='p',cex.main=1.6,cex.lab=1.8,cex.axis=1.2,
        cex=.4,main = expression(paste("Hysteresis for ", sigma, " = 4, 7, 10, 15")))

b=1
sigma=4
abline(v= MF01tau_critical1(s=sigma))
abline(v= MF01tau_critical2(s=sigma))
sigma=7
abline(v= MF01tau_critical1(s=sigma))
abline(v= MF01tau_critical2(s=sigma))
sigma=10
abline(v= MF01tau_critical1(s=sigma))
abline(v= MF01tau_critical2(s=sigma))
sigma=15
abline(v= MF01tau_critical1(s=sigma))
abline(v= MF01tau_critical2(s=sigma))

dev.off()

##Figure 5
#zooming in on sigma
pdf('figure5.pdf',h=6,w=10)#,paper='a4r')
layout(matrix(1:3,1,3))

beta=1
d=1
# bifurcation set
MF01tau_critical1=function(s,b=beta) 
{q=sqrt(1-4/(s*b))
(-s/2)*(1+q)+(1/b)*log((1+q)/(1-q))
}

MF01tau_critical2=function(s,b=beta) 
{q=sqrt(1-4/(s*b))
(-s/2)*(1-q)+(1/b)*log((1-q)/(1+q))
}

curve(MF01tau_critical1,4,20,ylim=c(-14,0),ylab=expression(tau),xlab=expression(sigma),bty='n',cex.main=1.6,cex.lab=1.8,cex.axis=1.2,
      main=expression(paste('Bifurcation set for ',sigma,',',tau,' (',beta,' = 1)')),lwd=2)
curve(MF01tau_critical2,4,20,add=T,lwd=2)
abline(h=-1,col='green')
abline(h=-3,col='red')
abline(h=-5,col='blue')
abline(h=-12,col='black')

text(15.5,-4.1,expression(paste('logarithmically to -', infinity)),cex=1.2)
text(16.5,-9.5,expression(paste('linearly to -', infinity)),cex=1.2)

# Panel 2: 1 dimn bif diagram sigma MF01
tau=-3
sigmas=seq(4,18,length=200)
dat=matrix(0,length(sigmas),4)
i=0
for(sigma in sigmas)
{
  i=i+1
  sol <- uniroot.all(MF01, c(-3, 3),n=2000)
  if(length(sol)==1) sol=c(sol, NA, NA)
  dat[i,]=c(sigma,sol)
}
matplot(dat[,1],dat[,-1],pch=1,col='red',bty='n',cex=c(.4,.1,.4),xlab=expression(sigma),ylab='p',ylim=0:1,cex.main=1.6,cex.lab=1.8,cex.axis=1.2,
        main=expression(paste('Bifurcation diagram for ',sigma)))

tau=-5
sigmas=seq(4,18,length=200)
dat=matrix(0,length(sigmas),4)
i=0
for(sigma in sigmas)
{
  i=i+1
  sol <- uniroot.all(MF01, c(-3, 3),n=2000)
  if(length(sol)==1) sol=c(sol, NA, NA)
  dat[i,]=c(sigma,sol)
}
matplot(dat[,1],dat[,-1],pch=1,col='blue',bty='n',cex=c(.4,.1,.4),xlab=expression(sigma),ylab=expression(mu),add=T)

tau=-1
sigmas=seq(4,18,length=200)
dat=matrix(0,length(sigmas),4)
i=0
for(sigma in sigmas)
{
  i=i+1
  sol <- uniroot.all(MF01, c(-3, 3),n=2000)
  if(length(sol)==1) sol=c(sol, NA, NA)
  dat[i,]=c(sigma,sol)
}
matplot(dat[,1],dat[,-1],pch=1,col='green',bty='n',cex=c(.4,.1,.4),xlab=expression(sigma),ylab='p',add=T)

# panel 3: potential function
V= function (m) (1/2)*m^2-(1/(sigma*beta))*log(1+exp(beta*(tau+sigma*m)))
MF01tau_critical1(15.684)

sigmas=c(15.684,seq(5,60,5))
i=0
beta=1
tau=-12
for(sigma in sigmas)
{
  i=i+1
  curve(V,0,1,add=i>1,ylim=c(-.3,.5),col=1+(i==1),lwd=1+(i==1),bty='n',ylab='V',xlab='p',
        main=expression(paste('Energy landscape ',tau,' = -12')),cex.main=1.6,cex.lab=1.8,cex.axis=1.2)#,xaxt = "n", yaxt = "n")
}
text(.6,.25,expression(paste(sigma,' = 5')),cex=1.2)
text(.44,-.2,expression(paste(sigma,' = 60')),cex=1.2)

points(.9,.43,pch=19,cex=3)
arrows(.9,.43,.8,.352,length = 0.1,col='blue',lwd=1.5)

points(.28,-.027,pch=19,cex=3)
arrows(.28,-.027,.38,-.09,length = 0.1,col='blue',lwd=1.5)


dev.off()

## Figure 6
pdf('figure6.pdf',h=5,w=10)#,paper='a4r')

layout(matrix(1:5,1,5,byrow=T))


MF01 <-function(p) 1/(1+exp(-beta*(sigma*p+tau)))-p
n=100 # network size
beta <- 1    # inverse temperature
tau <- -3   # external field
samples=1000

sigmas=seq(5,8,1)

curve(MF01tau_critical2,4,9,ylim=c(-4,-1.5),ylab=expression(tau),xlab=expression(sigma),bty='n',cex.main=1.6,cex.lab=1.8,cex.axis=1.2,
      main=expression(paste('Bifurcation set for ',sigma,',',tau)),lwd=2)
curve(MF01tau_critical1,4,6.5,add=T,lwd=2)
abline(h=-3,col='red')
points(x=sigmas,y=rep(-3,4))


for (sigma in sigmas)
{
network=matrix(sigma/(n-1),nrow=n,ncol=n); diag(network)=0
data=IsingSampler(samples, network, nIter = 100, thresholds=rep(tau,n), beta = beta, responses = c(0, 1))
h=hist(apply(data,1,sum),breaks =seq(0,n,2),main=bquote(sigma == .(sigma)),xlab=expression(sum(x)),freq=F,col='white',cex.main=1.8,cex.lab=1.8,cex.axis=1.2,ylab='')
equilibria <- uniroot.all(MF01, c(-1.1, 1.1),n=500)
equilibria
y_max <- max(h$counts)/samples
segments(x0 = equilibria*n, y0 = 0, x1 = equilibria*n   , y1 = -.1, col = c(2,1,2), lwd = 5)
}
dev.off()

#### extra
library(plotly)
MF01 <-function(p) 1/(1+exp(-beta*(sigma*p+tau)))-p
taus=seq(-25,0,length=140)
sigmas=seq(4,25,length=140)
dat=matrix(0,length(sigmas)*length(taus),5)
dat=matrix(NA,0,3)
i=0
for(tau in taus)
for(sigma in sigmas)
{
  i=i+1
  mu <- uniroot.all(MF01, c(-3, 3),n=2000)
  stable=3
  if(length(mu)==3) stable=c(3,1,3)
  dat=rbind(dat,data.frame(tau,sigma,mu,stable))
}

plot_ly(
  dat, x = ~-1*tau, y = ~sigma, z = ~mu,
  type = 'scatter3d', mode = 'markers',
  marker = list(size = 2, opacity = 0.6, color = ~stable,
                colorscale = list(c(0, 'blue'), c(1, 'green'), c(2, 'red')))) %>%
  layout(
    scene = list(
      xaxis = list(title = list(text = "τ", font = list(size = 18))),
      yaxis = list(title = list(text = "σ", font = list(size = 18))),
      zaxis = list(title = list(text = "μ", font = list(size = 18)))
    )
  )

