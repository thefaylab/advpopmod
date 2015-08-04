OutDir <- "F:\\"
InpDir <- "C:\\courses\\fish 559_14\\ADMB Workshop\\"

lectH<-function()
{
 par(mfrow=c(2,3))
 Xinit <- Case1()
# Xinit <- c(log(40),log(60),-1)

 # Poor choice
 sd <- c(0.1,0.1,0.1)
 cor <- matrix(c(1,0,0,0,1,0,0,0,1),ncol=3,nrow=3)
 
 # Good choice
 sd <- c(0.015702,0.0093949,0.17678)
 cor <- matrix(c(1,0.3769,0,0.3769,1,0,0,0,1),ncol=3,nrow=3)
 print(cor)
 DoMCMC(Xinit,3,sd,cor,Nsim=10000,Nburn=1,Nthin=5)
 
 Case3(paste(OutDir,"Res.Out",sep=""))
}
Case3<-function(FileN)
{
 TheData<-scan(FileN,what=list(A50=0,A95=0,Sigma=0,Post=0))
 par(mfrow=c(2,2))
 xx <- seq(1,length(TheData$A50))
 plot(xx,TheData$A50,xlab="Cycle number",ylab="Parameter 1",pch=16)
 plot(xx,TheData$A95,xlab="Cycle number",ylab="Parameter 2",pch=16)
 plot(xx,TheData$Sigma,xlab="Cycle number",ylab="Parameter 3",pch=16)
 plot(xx,-1*TheData$Post,xlab="Cycle number",ylab="Post Density",pch=16)
 dataout<-matrix(c(TheData$A50,TheData$A95,TheData$Sigma,TheData$Post),nrow=length(TheData$A50))
 pairs(dataout,labels=c("A50","A95","Sigma","Post Density"),pch=16)
}
Case1<-function()
{
 FileN <- paste(InpDir,"lectH.txt",sep="")    
 TheData<-scan(file=FileN,what=list(Length=0,Prob=0))
 co <- NULL
 co$Length <- TheData$Length
 co$Prob1 <- TheData$Prob
 co$Prob2 <- log(TheData$Prob/(1-TheData$Prob))
 assign("co",co,pos=1)

 xinit <- c(log(40),log(80),log(5))
 result<- optim(xinit,f1)   
 print(result)
 xfinal <- result$par
 Pred <- MakePred(TheData$Length,exp(xfinal[1]),exp(xfinal[2]))
 plot(TheData$Length,TheData$Prob,type='p',pch=16,lty=1,xlab="Length",ylab="Probability")
 lines(TheData$Length,Pred,lty=1,lwd=2)
 print(Pred)
 return(xfinal)
}

# =================================================================================================================

MakePred<-function(Length,A50,A95)
{
 pred <- 1/(1+exp(-1*log(19)*(Length-A50)/(A95-A50)))
 return(pred)   
}

# =================================================================================================================

f1 <- function(x)
{
 A50 <- exp(x[1])   
 A95 <- exp(x[2])   
 Sigma <- exp(x[3])
 pred <- MakePred(co$Length,A50,A95)
 pred <- log(pred/(1-pred))
 SS <- (pred-co$Prob2)^2     
 SS <- sum(SS)
 LogLike <- length(pred)*log(Sigma) + sum(SS)/(2.0*Sigma^2)
 return(LogLike)
}

# =================================================================================================================

DoMCMC<-function(Xinit,Ndim,sd,cor,Nsim=1000,Nburn=0,Nthin=1)
{
    
 library(mvtnorm)
 FileN <- paste(InpDir,"lectH.txt",sep="")    
 TheData<-scan(file=FileN,what=list(Length=0,Prob=0))
 co <- NULL
 co$Length <- TheData$Length
 co$Prob1 <- TheData$Prob
 co$Prob2 <- log(TheData$Prob/(1-TheData$Prob))
# assign("co",co,pos=1)
 environment(f1) = environment()

 covar <- matrix(0,ncol=Ndim,nrow=Ndim)
 print(covar)
 for (II in 1: Ndim)
  for (JJ in 1:Ndim)
   covar[II,JJ] <- cor[II,JJ]*sd[II]*sd[JJ] 
 print(covar)

 Xcurr <- Xinit
 Fcurr <- -1*f1(Xcurr)
 Outs <- matrix(0,nrow=(Nsim-Nburn),ncol=(Ndim+1))
 Ipnt <- 0; Icnt <- 0
 for (Isim in 1:Nsim)
  {
    Xnext <- rmvnorm(1, mean=Xcurr, sigma=covar)
    Fnext <- -1*f1(Xnext)
    Rand1 <- log(runif(1,0,1))
    if (Fnext > Fcurr+Rand1)
     {Fcurr <- Fnext; Xcurr <- Xnext }   
    if (Isim %% Nthin == 0)
     {
      Ipnt <- Ipnt + 1
      if (Ipnt > Nburn) { Icnt <- Icnt + 1; Outs[Icnt,] <- c(Xcurr,Fcurr); }    
     }
  } 
 xx <- seq(1,Icnt)
 for (II in 1:(Ndim+1))
  {
   yy <- Outs[,II][1:Icnt]
   if (II <= Ndim)
    lab1 <- paste("Parameter ",II)
   else
    lab1 <- "Posterior density"
   if (II <= Ndim) yy <- exp(yy)
   plot(xx,yy,xlab="Cycle number",ylab=lab1,type='b',pch=16)
  }
 Outs2 <- matrix(c(Outs[,1][1:Icnt],Outs[,2][1:Icnt],Outs[,3][1:Icnt],Outs[,4][1:Icnt]),nrow=Icnt)
 print(Outs2)
 pairs(Outs2,labels=c("A50","A95","Sigma","Post Density"),pch=16)
 
 write(t(Outs2),paste(OutDir,"Res.Out",sep=""),ncolumns=4)
 
}
lectH()
