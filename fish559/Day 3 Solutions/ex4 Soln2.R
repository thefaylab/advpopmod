FileN <- "C:\\courses\\FISH 559_14\\ClassEx\\Ex4\\EX4b.INP"

setwd("C:\\courses\\FISH 559_14\\ClassEx\\Ex4\\")


MinF <- 0
MaxF <- 4
for (II in 1:15)
 {
  Fbar <- (MinF+MaxF)/2.0
  
  write("# number of years",file=FileN) 
  write(25,file=,FileN,append=T)
  write("# fishing mortality",file=FileN,append=T) 
  write(Fbar,file=,FileN,append=T)
  shell("ex4b -mceval > res.out")

  resu <- read.table("res.out")
  Nyears <- length(resu[1,])
  Ndraw <- length(resu[,1])
  Nprob <- 0
  for (Isim in 1:Ndraw)
   if (Isim > 100 & resu[Isim,46] > 1000) Nprob <- Nprob + 1
  Nprob <- Nprob / (Ndraw-100)
  cat(Fbar,Nprob,"\n")
  if (Nprob > 0.5)
   MinF <- Fbar
  else
   MaxF <- Fbar
    
  
  
 }  


Stats <- matrix(0,nrow=Nyears,ncol=5)
for (Iyear in 1:Nyears)
{
  Stats[Iyear,] <- quantile(resu[101:Ndraw,Iyear],prob=c(0.05,0.25,0.5,0.75,0.95)) 
}  
print(Stats)

par(mfrow=c(2,2))
ymax <- max(Stats)*1.05
plot(1:Nyears,Stats[,3],xlab="Year",ylab="Population size",ylim=c(0,ymax),type='n',yaxs="i",xaxs="i")
xx <- c(c(1:Nyears),Nyears:1)
yy <- c(Stats[,1],rev(Stats[,5]))
polygon(xx,yy,col="gray50")
yy <- c(Stats[,2],rev(Stats[,4]))
polygon(xx,yy,col="gray95")
lines(1:Nyears,Stats[,3],lty=1,lwd=2)
abline(h=1000,lty=2)
