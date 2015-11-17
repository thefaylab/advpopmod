
setwd("~/classes/advpopmod/labs/lab8/")
source('lab8_functions.r')

#model <- "prodmodel"
model <- "tuna"
mod1 <- get.pars(model)

# plot the fit of the south atlantic albacore model (Polacheck et al 1993)
plot(0.5:23.5,mod1$bio[,2],type='l',lwd=3,xlab="Year",ylab="Biomass",ylim=c(0,300),
     xlim=c(0,24))
CatObs <- scan('tuna.dat',n=23,skip=23)
barplot(CatObs,add=TRUE,space=0)

#generate a new data set based on the fitted model
sigma <- mod1$par[6]
q <- mod1$par[5]
#assume that CPUE index is in the middle of the year (same as our assessment model, see prodmodel.tpl)
true.bio.pred <- 0.5*(mod1$bio[-1,2]+mod1$bio[-nrow(mod1$bio),2])
#rnorm() generates normally distributed random variables
new.dat <- q*true.bio.pred*exp(rnorm(length(pred.means),0,sigma)-0.5*sigma^2)
#add the new data to the previous plot
points(new.dat/q,pch=16,cex=0.8,col=rgb(0,0,1,1))


#write the new data file
BioObs <- cbind(1:23,format(new.dat,digits=5))
write.datfile(Nsp=1,BioObs,CatObs)
#have a look at new.dat

# Now fit the model to the new data set
# If using Windows, change which lines are commented out
#### WINDOWS
#command <- paste("C: & cd ",getwd()," & prodmodel -ind new.dat",sep="")
#shell(command,wait=TRUE,invisible=TRUE)
#### UNIX (Mac/linux)
command <- paste("cd ",getwd()," & ./prodmodel -ind new.dat",sep="")
system(command,wait=TRUE,show.output.on.console=FALSE)

#now read in the model estimates from the fit to the simulated data
model <- "prodmodel"
results <- get.pars(model)
str(results)
#the parameter estimates are in results$par, biomass in results$bio


# Now let's do a 100 simulations to look at estimation performance

#First set up some storage variables to put the results in.
est.pars <- NULL
gdt <- NULL
nll <- NULL
est.bio <- NULL

set.seed(42)
Nsim = 100
#Perform the following instructions Nsim (100) times.
for (isim in 1:Nsim)
{
  sigma <- mod1$par[6]
  q <- mod1$par[5]
  true.bio.pred <- 0.5*(mod1$bio[-1,2]+mod1$bio[-nrow(mod1$bio),2]) 
  new.dat <- q*true.bio.pred*exp(rnorm(length(pred.means),0,sigma)-0.5*sigma^2)

  BioObs <- cbind(1:23,format(new.dat,digits=5))
  write.datfile(Nsp=1,BioObs,CatObs)

  # Now fit the model to the new data set.
  # If using Windows, change which lines are commented out
  #### WINDOWS
  #command <- paste("C: & cd ",getwd()," & prodmodel -ind new.dat",sep="")
  #shell(command,wait=TRUE,invisible=TRUE)
  #### UNIX (Mac/linux)
  command <- paste("cd ",getwd()," & ./prodmodel -ind new.dat",sep="")
  system(command,wait=TRUE,show.output.on.console=FALSE)
  
  #now read in the results and save them.
  results <- get.pars("prodmodel")
  #NB, we are appending to the storage variables each time (hence rbind() )
  est.pars <- rbind(est.pars,results$par)
  gdt <- c(gdt,results$grad)
  nll <- c(nll,results$nll)
  est.bio <- cbind(est.bio,results$bio[,2])
}

#We might want to check if there were any convergence problems
pick <- which(gdt<=0.1)
print(length(pick))
#Perhaps only proceed with those simulations that 'worked'

#plot histograms of the parameter estimates (plus MSY)
par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(2,3,2,2))
hist(est.pars[pick,1],main="r",xlab="",xlim=c(0,0.6))
abline(v=mod1$par[1],lty=2)
hist(est.pars[pick,2],main="K",xlab="",xlim=c(0,600))
abline(v=mod1$par[2],lty=2)
hist(est.pars[pick,4],main="B1/K",xlab="",xlim=c(0.75,2))
abline(v=mod1$par[4],lty=2)
hist(est.pars[pick,1]*est.pars[pick,2]/4,main="MSY",xlab="",xlim=c(0,25))
abline(v=mod1$par[1]*mod1$par[2]/4,lty=2)

#calculate % Relative Error of parameter esimtates
par(mfrow=c(1,1))
ree <- est.pars[pick,]
ree <- cbind(ree,ree[,1]*ree[,2]/4)
true <- c(mod1$par,mod1$par[1]*mod1$par[2]/4)
for (irow in 1:nrow(ree))
 ree[irow,] <-(ree[irow,]-true)/ree[irow,]
boxplot(100*ree[,-3],names=c("r","K","B1/K","q","sigma","MSY"),ylim=c(-100,100),
        ylab="% REE")  
abline(h=0)
#look at correlations between r & K
pairs(ree[,1:2],labels=c("r","K"))




