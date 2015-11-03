get.pars <- function(model="prodmodel")
{
  repfile <- paste(model,".rep",sep="")
  parfile <- paste(model,".par",sep="")
  results <- NULL
  # negative log-likelihood
  results$nll <- scan(repfile,n=1)
  # gradient
  pardat <- read.table(parfile,nrows=1,col.names=1:16,header=FALSE,fill=TRUE,
                       comment.char="")
  results$grad <- pardat[1,16]
  #parameter estimates
  results$par <- scan(repfile,n=6,skip=1)
  #biomass estimates
  results$bio <- read.table(repfile,skip=7,header=FALSE)
  return(results)
}


write.datfile <- function(Nsp=1,BioObs=NULL,CatObs=1:10)
{
  fyear <- 1
  lyear <- length(CatObs)
  outfile <- "new.dat"
  write("#Nsp",outfile)
  write(1,outfile,append=TRUE)
  write("# r phase",outfile,append=TRUE)
  write(1,outfile,append=TRUE)
  write("# rinit",outfile,append=TRUE)
  write(0.4,outfile,append=TRUE)
  write("# k phase",outfile,append=TRUE)
  write(1,outfile,append=TRUE)
  write("# Kinit",outfile,append=TRUE)
  write(300,outfile,append=TRUE)
  write("# z phase",outfile,append=TRUE)
  write(-3,outfile,append=TRUE)
  write("# Z init",outfile,append=TRUE)
  write(2,outfile,append=TRUE)
  write("# theta phase",outfile,append=TRUE)
  write(2,outfile,append=TRUE)
  write("# Theta init",outfile,append=TRUE)
  write(0.9,outfile,append=TRUE)
  write("# fyear",outfile,append=TRUE)
  write(fyear,outfile,append=TRUE)
  write("# lyear",outfile,append=TRUE)
  write(lyear,outfile,append=TRUE)
  write("# catches",outfile,append=TRUE)
  write(CatObs,outfile,append=TRUE) #,row.names=FALSE,col.names=FALSE)
  write("# nbio",outfile,append=TRUE)
  write(nrow(BioObs),outfile,append=TRUE)
  write("# obs bio",outfile,append=TRUE)
  write.table(BioObs,outfile,append=TRUE,col.names=FALSE,
              row.names=FALSE,quote=FALSE)
}


