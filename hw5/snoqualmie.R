#
# 36 years data with precipitation measured in inches
#
y <- scan("snoqualmie.txt")
nodays <- rep(c(366,365,365,365),9)
sum(nodays) # check
#
# Look at June only data
z <- matrix(0,nrow=36,ncol=30)
daysum <- 152
for (i in 1:36){
    if (i>1) daysum <- daysum + nodays[i]
#    cat("daysum: ",daysum,"\n")
    z[i,] <- y[daysum+1:30]
}
n00 <- n01 <- n10 <- n11 <- 0
for (i in 1:36){
   for (j in 2:30){
      if ( z[i,j-1]==0 && z[i,j]==0) n00 <- n00+1
      if ( z[i,j-1]==0 && z[i,j]>0) n01 <- n01+1
      if ( z[i,j-1]>0 && z[i,j]==0) n10 <- n10+1
      if ( z[i,j-1]>0 && z[i,j]>0) n11 <- n11+1
   }
}
n <- matrix(c(n00, n01,n10,n11),nrow=2,ncol=2,byrow=T)
n


