setwd("C:\\courses\\FISH 559_14\\ClassEx\\Ex4\\")

resu <- read.table("resu.out")
print(head(resu))
Nyears <- length(resu[1,])
Ndraw <- length(resu[,1])
Stats <- matrix(0,nrow=Nyears,ncol=5)

for (Iyear in 1:Nyears)
 {
  Stats[Iyear,] <- quantile(resu[11:Ndraw,Iyear],prob=c(0.05,0.25,0.5,0.75,0.95)) 
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




