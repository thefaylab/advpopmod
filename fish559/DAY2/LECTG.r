setwd("C:\\courses\\FISH 559_14\\ADMB Workshop\\Examples\\")
dll <- "simpdll2.dll"
x <- runif(100,0,1)
y <- 3+5*x+rnorm(100,0,0.1)

# Need to create the variables
f <- 0
a <- 0
b <- 0


nsim <- 100
Results <- rep(NA,nsim)
for (II in 1:nsim)
{
 y <- 3+5*x+rnorm(100,0,1)
 
 # Load the dll
 dyn.load(dll)
 
 # call the dll
 xx <- .C("simpdll2", as.integer(length(x)), as.double(x), as.double(y),as.double(a),as.double(b),as.double(f), "-nox")
 
 # Store the results
 Results[II] <- xx[[4]]
 print(xx)
 AAA
 
 # upload the dll
 dyn.unload(dll)
} 
hist(Results)