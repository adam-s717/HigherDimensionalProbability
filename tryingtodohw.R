

normal <- as.data.frame(runif(10000))
u <- as.data.frame(runif(10000))
u[,c(1:5)] <- as.data.frame(runif(10000))
result <- as.data.frame(runif(10000))

for (i in 2:400){
  u <- as.data.frame(runif(10000))
  u[1:i] <- as.data.frame(runif(10000))
  normal <- as.data.frame(runif(10000))
  u[i] <-as.data.frame(runif(10000))
  x <- as.data.frame((1-2*u))
  xsq <- as.data.frame(x^2)
  normal <- as.data.frame(rowSums(xsq))
  result[i] <- as.data.frame(ifelse(normal$`rowSums(xsq)`<=1, 0, 1))
  remove(u,x,xsq,normal)
}
count <- as.data.frame(colSums(result))
count$numberrejected <- count$`colSums(result)`
count$`colSums(result)`<- NULL
plot(count$numberrejected)

a <-as.data.frame(runif(100))
for (i in 2:400){
  a[i] <-as.data.frame(runif(100))
}

d<- as.data.frame(runif(1))

for(j in 1:100){
  for (i in 2:400){
    d[j,i] <- as.data.frame(a[j,1]-a[j,i])
  }
}

dsq <- as.data.frame(d^2)
dsq$`runif(1)` <- NULL

d1_else <- as.data.frame(colSums(dsq))
d1_else <- sqrt(d1_else)
hist(d1_else$`colSums(dsq)`)
