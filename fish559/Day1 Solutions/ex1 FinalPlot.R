setwd("C:\\courses\\FISH 559_14\\ADMB Workshop\\Test\\")

TheData <- read.table("ex1a.dat",skip=3)
print(TheData)

ThePred1 <- read.table("ex1a.rep",skip=1)
ThePred2 <- read.table("ex1b.rep",skip=1)

par(mfrow=c(2,2))
plot(TheData[,1],TheData[,2],xlab="Age",ylab="Length",pch=16)
lines(ThePred1[,1],ThePred1[,2],lty=1)
lines(ThePred2[,1],ThePred2[,2],lty=2)
legend("topleft",legend=c("Logistic","Von Bertlanffy"),lty=1:2)

