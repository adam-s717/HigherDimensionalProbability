library(R.matlab)
library(ggplot2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#      Importing The Data      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

data1 <- readMat('/Users/Adam/Downloads/genomics/data1.mat')
data1 <- as.data.frame(data1)
data1 <- data.matrix(data1)

data2 <- readMat('/Users/Adam/Downloads/genomics/data2.mat')
data2 <- as.data.frame(data2)
data2 <- data.matrix(data2)

data3 <- readMat('/Users/Adam/Downloads/genomics/data3.mat')
data3 <- as.data.frame(data3)
data3 <- data.matrix(data3)

data4 <- readMat('/Users/Adam/Downloads/genomics/data4.mat')
data4 <- as.data.frame(data4)
data4 <- data.matrix(data4)

data5 <- readMat('/Users/Adam/Downloads/genomics/data5.mat')
data5 <- as.data.frame(data5)
data5 <- data.matrix(data5)

data6 <- readMat('/Users/Adam/Downloads/genomics/data6.mat')
data6 <- as.data.frame(data6)
data6 <- data.matrix(data6)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#      Creating Our Algos      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


JLT <- function(X,k){
  A = matrix(rnorm(n=k*nrow(X), mean = 0, sd = 1), ncol = nrow(X))
  Y = A %*% X
  Y = (1/sqrt(k))*Y
  return(Y)
}

FJLT <- function(X,k){
  A = matrix(sample(x = c(1,-1,0), prob = c((1/6),(1/6),(2/3)), size = k*nrow(X), replace = TRUE)
             , ncol = nrow(X))
  Y = A %*% X
  Y = sqrt(3/k)*Y
  return(Y)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#      Running the FJLT Algo    #~~~~~~~~~~~~#
#~~~~~~~~~~~~#         on data set 1       #~~~~~~~~~~~~#
#~~~~~~~~~~~~#  For k = 25, 100, 225, 400   #~~~~~~~~~~~~#
#~~~~~~~~~~~~# and recording the  run time  #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

time25 <- system.time(Y_s1_25 <- FJLT(data1, 25)) 
time100 <- system.time(Y_s1_100 <- FJLT(data1, 100))
time225 <- system.time(Y_s1_225 <- FJLT(data1, 225))
time400 <- system.time(Y_s1_400 <- FJLT(data1, 400))

runtime  = matrix(nrow = 4, ncol = 1)
runtime <- as.data.frame(runtime)
runtime <- data.matrix(runtime)

runtime[1,1] <- time25[3]
runtime[2,1] <- time100[3]
runtime[3,1] <- time225[3]
runtime[4,1] <- time400[3]

#system.time(timetest <- JLT(data1, 25))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS1, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C25 = matrix(nrow = ncol(Y_s1_25), ncol = ncol(Y_s1_25))
C25 <- as.data.frame(C25)
C25  <- data.matrix(C25)

j=1
while(j <= ncol(Y_s1_25)){
  l <- 1
  while (l <= ncol(Y_s1_25)){
    i = 1:nrow(Y_s1_25)
    C25[j,l] <- sum((Y_s1_25[i,j]-Y_s1_25[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

##################
#   Creating D   #
##################
D = matrix(nrow = ncol(data1), ncol = ncol(data1))
D <- as.data.frame(D)
D <- data.matrix(D)


for(j in 1:ncol(data1)){
  l <- 1
  while (l <= ncol(data1)){
    i = 1:nrow(data1)
    D[j,l] <- sum((data1[i,j]-data1[i,l])^2)
  
    l=l+1
  }
}


dat25 <- data.frame(d=as.vector(D[upper.tri(D)]),
                  c25=as.vector(C25[upper.tri(C25)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS1, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat25$r <- dat25$c/dat25$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS1, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 48)
j <- j[43:47]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C25,index.return = TRUE)$ix, 48)
l <- l[43:47]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C25 == C25[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C25 == C25[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C25 == C25[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C25 == C25[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C25 == C25[l[5]], arr.ind = TRUE)[1,2]

nn1_25 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
   nn1_25 = nn1_25+1
  }else{nn1_25 = nn1_25+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn1_25 = nn1_25/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS1, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C100 = matrix(nrow = ncol(Y_s1_100), ncol = ncol(Y_s1_100))
C100 <- as.data.frame(C100)
C100 <- data.matrix(C100)

j=1
while(j <= ncol(Y_s1_100)){
  l <- 1
  while (l <= ncol(Y_s1_100)){
    i = 1:nrow(Y_s1_100)
    C100[j,l] <- sum((Y_s1_100[i,j]-Y_s1_100[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

##################
#   Creating D   #
##################
D = matrix(nrow = ncol(data1), ncol = ncol(data1))
D <- as.data.frame(D)
D <- data.matrix(D)


for(j in 1:ncol(data1)){
  l <- 1
  while (l <= ncol(data1)){
    i = 1:nrow(data1)
    D[j,l] <- sum((data1[i,j]-data1[i,l])^2)
    
    l=l+1
  }
}


dat100 <- data.frame(d=as.vector(D[upper.tri(D)]),
                  c100=as.vector(C100[upper.tri(C100)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS1, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat100$r <- dat100$c/dat100$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS1, k=100     #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 48)
j <- j[43:47]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C100,index.return = TRUE)$ix, 48)
l <- l[43:47]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C100 == C100[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C100 == C100[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C100 == C100[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C100 == C100[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C100 == C100[l[5]], arr.ind = TRUE)[1,2]

nn1_100 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn1_100 = nn1_100+1
  }else{nn1_100 = nn1_100+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn1_100 = nn1_100/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS1, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C225 = matrix(nrow = ncol(Y_s1_225), ncol = ncol(Y_s1_225))
C225 <- as.data.frame(C225)
C225 <- data.matrix(C225)

j=1
while(j <= ncol(Y_s1_225)){
  l <- 1
  while (l <= ncol(Y_s1_225)){
    i = 1:nrow(Y_s1_225)
    C225[j,l] <- sum((Y_s1_225[i,j]-Y_s1_225[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

##################
#   Creating D   #
##################
D = matrix(nrow = ncol(data1), ncol = ncol(data1))
D <- as.data.frame(D)
D <- data.matrix(D)


for(j in 1:ncol(data1)){
  l <- 1
  while (l <= ncol(data1)){
    i = 1:nrow(data1)
    D[j,l] <- sum((data1[i,j]-data1[i,l])^2)
    
    l=l+1
  }
}


dat225 <- data.frame(d=as.vector(D[upper.tri(D)]),
                  c225=as.vector(C225[upper.tri(C225)]))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS1, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat225$r <- dat225$c/dat225$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS1, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 48)
j <- j[43:47]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C225,index.return = TRUE)$ix, 48)
l <- l[43:47]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C225 == C225[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C225 == C225[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C225 == C225[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C225 == C225[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C225 == C225[l[5]], arr.ind = TRUE)[1,2]

nn1_225 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn1_225 = nn1_225+1
  }else{nn1_225 = nn1_225+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn1_225 = nn1_225/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS1, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C400 = matrix(nrow = ncol(Y_s1_400), ncol = ncol(Y_s1_400))
C400 <- as.data.frame(C400)
C400 <- data.matrix(C400)

j=1
while(j <= ncol(Y_s1_400)){
  l <- 1
  while (l <= ncol(Y_s1_400)){
    i = 1:nrow(Y_s1_400)
    C400[j,l] <- sum((Y_s1_400[i,j]-Y_s1_400[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

##################
#   Creating D   #
##################
D = matrix(nrow = ncol(data1), ncol = ncol(data1))
D <- as.data.frame(D)
D <- data.matrix(D)


for(j in 1:ncol(data1)){
  l <- 1
  while (l <= ncol(data1)){
    i = 1:nrow(data1)
    D[j,l] <- sum((data1[i,j]-data1[i,l])^2)
    
    l=l+1
  }
}


dat400 <- data.frame(d=as.vector(D[upper.tri(D)]),
                  c400=as.vector(C400[upper.tri(C400)]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS1, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat400$r <- dat400$c/dat400$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS1, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 48)
j <- j[43:47]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C400,index.return = TRUE)$ix, 48)
l <- l[43:47]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C400 == C400[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C400 == C400[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C400 == C400[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C400 == C400[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C400 == C400[l[5]], arr.ind = TRUE)[1,2]

nn1_400 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn1_400 = nn1_400+1
  }else{nn1_400 = nn1_400+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn1_400 = nn1_400/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 1 for DS1 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
#k=25
plot(dat25$d, dat25$c25, main = "DS1, K=25")
#k=100
plot(dat100$d, dat100$c100, main = "DS1, K=100")
#k=225
plot(dat225$d, dat225$c225, main = "DS1, K=225")
#k=400
plot(dat400$d, dat400$c400, main = "DS1, K=400")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 2 for DS1 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
require(gridExtra)
#k=25
p1 <- ggplot(dat25, aes(x=r)) + geom_histogram()+ggtitle("DS1 K=25")
#k=100 
p2 <- ggplot(dat100, aes(x=r)) + geom_histogram()+ggtitle("DS1 K=100")
#k=225 
p3 <- ggplot(dat225, aes(x=r)) + geom_histogram() +ggtitle("DS1 K=225")
#k=400
p4 <- ggplot(dat400, aes(x=r)) + geom_histogram() +ggtitle("DS1 K=400")

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#      Running the FJLT Algo    #~~~~~~~~~~~~#
#~~~~~~~~~~~~#         on data set 2       #~~~~~~~~~~~~#
#~~~~~~~~~~~~#  For k = 25, 100, 225, 400   #~~~~~~~~~~~~#
#~~~~~~~~~~~~# and recording the  run time  #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

time25 <- system.time(Y_s2_25 <- FJLT(data2, 25)) 
time100 <- system.time(Y_s2_100 <- FJLT(data2, 100))
time225 <- system.time(Y_s2_225 <- FJLT(data2, 225))
time400 <- system.time(Y_s2_400 <- FJLT(data2, 400))

runtime2  = matrix(nrow = 4, ncol = 1)
runtime2 <- as.data.frame(runtime)
runtime2 <- data.matrix(runtime)

runtime2[1,1] <- time25[3]
runtime2[2,1] <- time100[3]
runtime2[3,1] <- time225[3]
runtime2[4,1] <- time400[3]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS2, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C25 = matrix(nrow = ncol(Y_s2_25), ncol = ncol(Y_s2_25))
C25 <- as.data.frame(C25)
C25  <- data.matrix(C25)

j=1
while(j <= ncol(Y_s2_25)){
  l <- 1
  while (l <= ncol(Y_s2_25)){
    i = 1:nrow(Y_s2_25)
    C25[j,l] <- sum((Y_s2_25[i,j]-Y_s2_25[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

##################
#   Creating D   #
##################
D = matrix(nrow = ncol(data2), ncol = ncol(data2))
D <- as.data.frame(D)
D <- data.matrix(D)


for(j in 1:ncol(data2)){
  l <- 1
  while (l <= ncol(data2)){
    i = 1:nrow(data2)
    D[j,l] <- sum((data2[i,j]-data2[i,l])^2)
    
    l=l+1
  }
}


dat25 <- data.frame(d=as.vector(D[upper.tri(D)]),
                    c25=as.vector(C25[upper.tri(C25)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS2, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat25$r <- dat25$c/dat25$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS1, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 68)
j <- j[63:67]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C25,index.return = TRUE)$ix, 68)
l <- l[63:67]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C25 == C25[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C25 == C25[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C25 == C25[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C25 == C25[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C25 == C25[l[5]], arr.ind = TRUE)[1,2]

nn2_25 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn1_25 = nn2_25+1
  }else{nn2_25 = nn2_25+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn2_25 = nn2_25/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS2, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C100 = matrix(nrow = ncol(Y_s2_100), ncol = ncol(Y_s2_100))
C100 <- as.data.frame(C100)
C100 <- data.matrix(C100)

j=1
while(j <= ncol(Y_s2_100)){
  l <- 1
  while (l <= ncol(Y_s2_100)){
    i = 1:nrow(Y_s2_100)
    C100[j,l] <- sum((Y_s2_100[i,j]-Y_s2_100[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat100 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c100=as.vector(C100[upper.tri(C100)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS2, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat100$r <- dat100$c/dat100$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS2, k=100     #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 68)
j <- j[63:67]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C100,index.return = TRUE)$ix, 68)
l <- l[63:67]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C100 == C100[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C100 == C100[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C100 == C100[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C100 == C100[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C100 == C100[l[5]], arr.ind = TRUE)[1,2]

nn2_100 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn2_100 = nn2_100+1
  }else{nn2_100 = nn2_100+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn2_100 = nn2_100/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS2, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C225 = matrix(nrow = ncol(Y_s2_225), ncol = ncol(Y_s2_225))
C225 <- as.data.frame(C225)
C225 <- data.matrix(C225)

j=1
while(j <= ncol(Y_s2_225)){
  l <- 1
  while (l <= ncol(Y_s2_225)){
    i = 1:nrow(Y_s2_225)
    C225[j,l] <- sum((Y_s2_225[i,j]-Y_s2_225[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat225 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c225=as.vector(C225[upper.tri(C225)]))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS2, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat225$r <- dat225$c/dat225$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS1, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 68)
j <- j[63:67]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C225,index.return = TRUE)$ix, 68)
l <- l[63:67]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C225 == C225[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C225 == C225[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C225 == C225[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C225 == C225[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C225 == C225[l[5]], arr.ind = TRUE)[1,2]

nn2_225 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn2_225 = nn2_225+1
  }else{nn2_225 = nn2_225+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn2_225 = nn2_225/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS2, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C400 = matrix(nrow = ncol(Y_s2_400), ncol = ncol(Y_s2_400))
C400 <- as.data.frame(C400)
C400 <- data.matrix(C400)

j=1
while(j <= ncol(Y_s2_400)){
  l <- 1
  while (l <= ncol(Y_s2_400)){
    i = 1:nrow(Y_s2_400)
    C400[j,l] <- sum((Y_s2_400[i,j]-Y_s2_400[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat400 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c400=as.vector(C400[upper.tri(C400)]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS2, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat400$r <- dat400$c/dat400$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS2, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 68)
j <- j[63:67]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C400,index.return = TRUE)$ix, 68)
l <- l[63:67]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C400 == C400[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C400 == C400[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C400 == C400[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C400 == C400[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C400 == C400[l[5]], arr.ind = TRUE)[1,2]

nn2_400 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn2_400 = nn2_400+1
  }else{nn2_400 = nn2_400+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn2_400 = nn2_400/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 1 for DS1 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
#k=25
plot(dat25$d, dat25$c25, main = "DS2, K=25")
#k=100
plot(dat100$d, dat100$c100, main = "DS2, K=100")
#k=225
plot(dat225$d, dat225$c225, main = "DS2, K=225")
#k=400
plot(dat400$d, dat400$c400, main = "DS2, K=400")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 2 for DS1 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
require(gridExtra)
#k=25
p1 <- ggplot(dat25, aes(x=r)) + geom_histogram()+ggtitle("DS2 K=25")
#k=100 
p2 <- ggplot(dat100, aes(x=r)) + geom_histogram()+ggtitle("DS2 K=100")
#k=225 
p3 <- ggplot(dat225, aes(x=r)) + geom_histogram() +ggtitle("DS2 K=225")
#k=400
p4 <- ggplot(dat400, aes(x=r)) + geom_histogram() +ggtitle("DS2 K=400")

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#      Running the FJLT Algo    #~~~~~~~~~~~~#
#~~~~~~~~~~~~#         on data set 3       #~~~~~~~~~~~~#
#~~~~~~~~~~~~#  For k = 25, 100, 225, 400   #~~~~~~~~~~~~#
#~~~~~~~~~~~~# and recording the  run time  #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

time25 <- system.time(Y_s3_25 <- FJLT(data3, 25)) 
time100 <- system.time(Y_s3_100 <- FJLT(data3, 100))
time225 <- system.time(Y_s3_225 <- FJLT(data3, 225))
time400 <- system.time(Y_s3_400 <- FJLT(data3, 400))

runtime3  = matrix(nrow = 4, ncol = 1)
runtime3 <- as.data.frame(runtime3)
runtime3 <- data.matrix(runtime3)

runtime3[1,1] <- time25[3]
runtime3[2,1] <- time100[3]
runtime3[3,1] <- time225[3]
runtime3[4,1] <- time400[3]








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS3, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C25 = matrix(nrow = ncol(Y_s3_25), ncol = ncol(Y_s3_25))
C25 <- as.data.frame(C25)
C25  <- data.matrix(C25)

j=1
while(j <= ncol(Y_s3_25)){
  l <- 1
  while (l <= ncol(Y_s3_25)){
    i = 1:nrow(Y_s3_25)
    C25[j,l] <- sum((Y_s3_25[i,j]-Y_s3_25[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

##################
#   Creating D   #
##################
D = matrix(nrow = ncol(data3), ncol = ncol(data3))
D <- as.data.frame(D)
D <- data.matrix(D)


for(j in 1:ncol(data3)){
  l <- 1
  while (l <= ncol(data3)){
    i = 1:nrow(data3)
    D[j,l] <- sum((data3[i,j]-data3[i,l])^2)
    
    l=l+1
  }
}


dat25 <- data.frame(d=as.vector(D[upper.tri(D)]),
                    c25=as.vector(C25[upper.tri(C25)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS3, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat25$r <- dat25$c/dat25$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS3, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 78)
j <- j[73:77]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C25,index.return = TRUE)$ix, 78)
l <- l[73:77]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C25 == C25[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C25 == C25[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C25 == C25[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C25 == C25[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C25 == C25[l[5]], arr.ind = TRUE)[1,2]

nn3_25 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn3_25 = nn3_25+1
  }else{nn3_25 = nn3_25+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn3_25 = nn3_25/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS2, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C100 = matrix(nrow = ncol(Y_s3_100), ncol = ncol(Y_s3_100))
C100 <- as.data.frame(C100)
C100 <- data.matrix(C100)

j=1
while(j <= ncol(Y_s3_100)){
  l <- 1
  while (l <= ncol(Y_s3_100)){
    i = 1:nrow(Y_s3_100)
    C100[j,l] <- sum((Y_s3_100[i,j]-Y_s3_100[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat100 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c100=as.vector(C100[upper.tri(C100)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS3, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat100$r <- dat100$c/dat100$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS3, k=100     #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 78)
j <- j[73:77]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C100,index.return = TRUE)$ix, 78)
l <- l[73:77]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C100 == C100[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C100 == C100[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C100 == C100[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C100 == C100[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C100 == C100[l[5]], arr.ind = TRUE)[1,2]

nn3_100 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn3_100 = nn3_100+1
  }else{nn3_100 = nn3_100+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn3_100 = nn3_100/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS3, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C225 = matrix(nrow = ncol(Y_s3_225), ncol = ncol(Y_s3_225))
C225 <- as.data.frame(C225)
C225 <- data.matrix(C225)

j=1
while(j <= ncol(Y_s3_225)){
  l <- 1
  while (l <= ncol(Y_s3_225)){
    i = 1:nrow(Y_s3_225)
    C225[j,l] <- sum((Y_s3_225[i,j]-Y_s3_225[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat225 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c225=as.vector(C225[upper.tri(C225)]))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS3, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat225$r <- dat225$c/dat225$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS3, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 78)
j <- j[73:77]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C225,index.return = TRUE)$ix, 78)
l <- l[73:77]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C225 == C225[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C225 == C225[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C225 == C225[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C225 == C225[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C225 == C225[l[5]], arr.ind = TRUE)[1,2]

nn3_225 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn3_225 = nn3_225+1
  }else{nn3_225 = nn3_225+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn3_225 = nn3_225/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS3, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C400 = matrix(nrow = ncol(Y_s3_400), ncol = ncol(Y_s3_400))
C400 <- as.data.frame(C400)
C400 <- data.matrix(C400)

j=1
while(j <= ncol(Y_s3_400)){
  l <- 1
  while (l <= ncol(Y_s3_400)){
    i = 1:nrow(Y_s3_400)
    C400[j,l] <- sum((Y_s3_400[i,j]-Y_s3_400[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat400 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c400=as.vector(C400[upper.tri(C400)]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS3, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat400$r <- dat400$c/dat400$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS3, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 78)
j <- j[73:77]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C400,index.return = TRUE)$ix, 78)
l <- l[73:77]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C400 == C400[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C400 == C400[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C400 == C400[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C400 == C400[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C400 == C400[l[5]], arr.ind = TRUE)[1,2]

nn3_400 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn3_400 = nn3_400+1
  }else{nn3_400 = nn3_400+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn3_400 = nn3_400/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 1 for DS3 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
#k=25
plot(dat25$d, dat25$c25, main = "DS3, K=25")
#k=100
plot(dat100$d, dat100$c100, main = "DS3, K=100")
#k=225
plot(dat225$d, dat225$c225, main = "DS3, K=225")
#k=400
plot(dat400$d, dat400$c400, main = "DS3, K=400")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 2 for DS3 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
require(gridExtra)
#k=25
p1 <- ggplot(dat25, aes(x=r)) + geom_histogram()+ggtitle("DS3 K=25")
#k=100 
p2 <- ggplot(dat100, aes(x=r)) + geom_histogram()+ggtitle("DS3 K=100")
#k=225 
p3 <- ggplot(dat225, aes(x=r)) + geom_histogram() +ggtitle("DS3 K=225")
#k=400
p4 <- ggplot(dat400, aes(x=r)) + geom_histogram() +ggtitle("DS3 K=400")

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)







#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#      Running the FJLT Algo    #~~~~~~~~~~~~#
#~~~~~~~~~~~~#         on data set 4       #~~~~~~~~~~~~#
#~~~~~~~~~~~~#  For k = 25, 100, 225, 400   #~~~~~~~~~~~~#
#~~~~~~~~~~~~# and recording the  run time  #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

time25 <- system.time(Y_s4_25 <- FJLT(data4, 25)) 
time100 <- system.time(Y_s4_100 <- FJLT(data4, 100))
time225 <- system.time(Y_s4_225 <- FJLT(data4, 225))
time400 <- system.time(Y_s4_400 <- FJLT(data4, 400))

runtime4  = matrix(nrow = 4, ncol = 1)
runtime4 <- as.data.frame(runtime4)
runtime4 <- data.matrix(runtime4)

runtime4[1,1] <- time25[3]
runtime4[2,1] <- time100[3]
runtime4[3,1] <- time225[3]
runtime4[4,1] <- time400[3]








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS4, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C25 = matrix(nrow = ncol(Y_s4_25), ncol = ncol(Y_s4_25))
C25 <- as.data.frame(C25)
C25  <- data.matrix(C25)

j=1
while(j <= ncol(Y_s4_25)){
  l <- 1
  while (l <= ncol(Y_s4_25)){
    i = 1:nrow(Y_s4_25)
    C25[j,l] <- sum((Y_s4_25[i,j]-Y_s4_25[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

##################
#   Creating D   #
##################
D = matrix(nrow = ncol(data4), ncol = ncol(data4))
D <- as.data.frame(D)
D <- data.matrix(D)


for(j in 1:ncol(data4)){
  l <- 1
  while (l <= ncol(data4)){
    i = 1:nrow(data4)
    D[j,l] <- sum((data4[i,j]-data4[i,l])^2)
    
    l=l+1
  }
}


dat25 <- data.frame(d=as.vector(D[upper.tri(D)]),
                    c25=as.vector(C25[upper.tri(C25)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS4, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat25$r <- dat25$c/dat25$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS4, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 68)
j <- j[63:67]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C25,index.return = TRUE)$ix, 68)
l <- l[63:67]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C25 == C25[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C25 == C25[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C25 == C25[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C25 == C25[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C25 == C25[l[5]], arr.ind = TRUE)[1,2]

nn4_25 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn4_25 = nn4_25+1
  }else{nn4_25 = nn4_25+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn4_25 = nn4_25/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS4, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C100 = matrix(nrow = ncol(Y_s4_100), ncol = ncol(Y_s4_100))
C100 <- as.data.frame(C100)
C100 <- data.matrix(C100)

j=1
while(j <= ncol(Y_s4_100)){
  l <- 1
  while (l <= ncol(Y_s4_100)){
    i = 1:nrow(Y_s4_100)
    C100[j,l] <- sum((Y_s4_100[i,j]-Y_s4_100[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat100 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c100=as.vector(C100[upper.tri(C100)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS4, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat100$r <- dat100$c/dat100$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS3, k=100     #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 68)
j <- j[63:67]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C100,index.return = TRUE)$ix, 68)
l <- l[63:67]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C100 == C100[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C100 == C100[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C100 == C100[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C100 == C100[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C100 == C100[l[5]], arr.ind = TRUE)[1,2]

nn4_100 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn4_100 = nn4_100+1
  }else{nn4_100 = nn4_100+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn4_100 = nn4_100/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS4, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C225 = matrix(nrow = ncol(Y_s4_225), ncol = ncol(Y_s4_225))
C225 <- as.data.frame(C225)
C225 <- data.matrix(C225)

j=1
while(j <= ncol(Y_s4_225)){
  l <- 1
  while (l <= ncol(Y_s4_225)){
    i = 1:nrow(Y_s4_225)
    C225[j,l] <- sum((Y_s4_225[i,j]-Y_s4_225[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat225 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c225=as.vector(C225[upper.tri(C225)]))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS4, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat225$r <- dat225$c/dat225$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS4, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 68)
j <- j[63:67]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C225,index.return = TRUE)$ix, 68)
l <- l[63:67]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C225 == C225[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C225 == C225[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C225 == C225[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C225 == C225[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C225 == C225[l[5]], arr.ind = TRUE)[1,2]

nn4_225 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn4_225 = nn4_225+1
  }else{nn4_225 = nn4_225+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn4_225 = nn4_225/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS4, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C400 = matrix(nrow = ncol(Y_s4_400), ncol = ncol(Y_s4_400))
C400 <- as.data.frame(C400)
C400 <- data.matrix(C400)

j=1
while(j <= ncol(Y_s4_400)){
  l <- 1
  while (l <= ncol(Y_s4_400)){
    i = 1:nrow(Y_s4_400)
    C400[j,l] <- sum((Y_s4_400[i,j]-Y_s4_400[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat400 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c400=as.vector(C400[upper.tri(C400)]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS4, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat400$r <- dat400$c/dat400$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS4, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 68)
j <- j[63:67]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C400,index.return = TRUE)$ix, 68)
l <- l[63:67]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C400 == C400[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C400 == C400[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C400 == C400[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C400 == C400[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C400 == C400[l[5]], arr.ind = TRUE)[1,2]

nn4_400 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn4_400 = nn4_400+1
  }else{nn4_400 = nn4_400+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn4_400 = nn4_400/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 1 for DS4 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
#k=25
plot(dat25$d, dat25$c25, main = "DS4, K=25")
#k=100
plot(dat100$d, dat100$c100, main = "DS4, K=100")
#k=225
plot(dat225$d, dat225$c225, main = "DS4, K=225")
#k=400
plot(dat400$d, dat400$c400, main = "DS4, K=400")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 2 for DS4 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
require(gridExtra)
#k=25
p1 <- ggplot(dat25, aes(x=r)) + geom_histogram()+ggtitle("DS4 K=25")
#k=100 
p2 <- ggplot(dat100, aes(x=r)) + geom_histogram()+ggtitle("DS4 K=100")
#k=225 
p3 <- ggplot(dat225, aes(x=r)) + geom_histogram() +ggtitle("DS4 K=225")
#k=400
p4 <- ggplot(dat400, aes(x=r)) + geom_histogram() +ggtitle("DS4 K=400")

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#      Running the FJLT Algo    #~~~~~~~~~~~~#
#~~~~~~~~~~~~#         on data set 5       #~~~~~~~~~~~~#
#~~~~~~~~~~~~#  For k = 25, 100, 225, 400   #~~~~~~~~~~~~#
#~~~~~~~~~~~~# and recording the  run time  #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

time25 <- system.time(Y_s5_25 <- FJLT(data5, 25)) 
time100 <- system.time(Y_s5_100 <- FJLT(data5, 100))
time225 <- system.time(Y_s5_225 <- FJLT(data5, 225))
time400 <- system.time(Y_s5_400 <- FJLT(data5, 400))

runtime5  = matrix(nrow = 4, ncol = 1)
runtime5 <- as.data.frame(runtime5)
runtime5 <- data.matrix(runtime5)

runtime5[1,1] <- time25[3]
runtime5[2,1] <- time100[3]
runtime5[3,1] <- time225[3]
runtime5[4,1] <- time400[3]








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS5, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C25 = matrix(nrow = ncol(Y_s5_25), ncol = ncol(Y_s5_25))
C25 <- as.data.frame(C25)
C25  <- data.matrix(C25)

j=1
while(j <= ncol(Y_s5_25)){
  l <- 1
  while (l <= ncol(Y_s5_25)){
    i = 1:nrow(Y_s5_25)
    C25[j,l] <- sum((Y_s5_25[i,j]-Y_s5_25[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

##################
#   Creating D   #
##################
D = matrix(nrow = ncol(data5), ncol = ncol(data5))
D <- as.data.frame(D)
D <- data.matrix(D)


for(j in 1:ncol(data5)){
  l <- 1
  while (l <= ncol(data5)){
    i = 1:nrow(data5)
    D[j,l] <- sum((data5[i,j]-data5[i,l])^2)
    
    l=l+1
  }
}


dat25 <- data.frame(d=as.vector(D[upper.tri(D)]),
                    c25=as.vector(C25[upper.tri(C25)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS5, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat25$r <- dat25$c/dat25$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS5, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 108)
j <- j[103:107]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C25,index.return = TRUE)$ix, 108)
l <- l[103:107]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C25 == C25[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C25 == C25[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C25 == C25[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C25 == C25[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C25 == C25[l[5]], arr.ind = TRUE)[1,2]

nn5_25 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn5_25 = nn5_25+1
  }else{nn5_25 = nn5_25+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn5_25 = nn5_25/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS5, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C100 = matrix(nrow = ncol(Y_s5_100), ncol = ncol(Y_s5_100))
C100 <- as.data.frame(C100)
C100 <- data.matrix(C100)

j=1
while(j <= ncol(Y_s5_100)){
  l <- 1
  while (l <= ncol(Y_s5_100)){
    i = 1:nrow(Y_s5_100)
    C100[j,l] <- sum((Y_s5_100[i,j]-Y_s5_100[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat100 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c100=as.vector(C100[upper.tri(C100)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS5, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat100$r <- dat100$c/dat100$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS5, k=100     #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 108)
j <- j[103:107]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C100,index.return = TRUE)$ix, 108)
l <- l[103:107]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C100 == C100[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C100 == C100[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C100 == C100[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C100 == C100[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C100 == C100[l[5]], arr.ind = TRUE)[1,2]

nn5_100 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn5_100 = nn5_100+1
  }else{nn5_100 = nn5_100+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn5_100 = nn5_100/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS5, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C225 = matrix(nrow = ncol(Y_s5_225), ncol = ncol(Y_s5_225))
C225 <- as.data.frame(C225)
C225 <- data.matrix(C225)

j=1
while(j <= ncol(Y_s5_225)){
  l <- 1
  while (l <= ncol(Y_s5_225)){
    i = 1:nrow(Y_s5_225)
    C225[j,l] <- sum((Y_s5_225[i,j]-Y_s5_225[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat225 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c225=as.vector(C225[upper.tri(C225)]))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS5, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat225$r <- dat225$c/dat225$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS5, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 108)
j <- j[103:107]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C225,index.return = TRUE)$ix, 108)
l <- l[103:107]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C225 == C225[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C225 == C225[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C225 == C225[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C225 == C225[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C225 == C225[l[5]], arr.ind = TRUE)[1,2]

nn5_225 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn5_225 = nn5_225+1
  }else{nn5_225 = nn5_225+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn5_225 = nn5_225/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS5, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C400 = matrix(nrow = ncol(Y_s5_400), ncol = ncol(Y_s5_400))
C400 <- as.data.frame(C400)
C400 <- data.matrix(C400)

j=1
while(j <= ncol(Y_s5_400)){
  l <- 1
  while (l <= ncol(Y_s5_400)){
    i = 1:nrow(Y_s5_400)
    C400[j,l] <- sum((Y_s5_400[i,j]-Y_s5_400[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat400 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c400=as.vector(C400[upper.tri(C400)]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS5, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat400$r <- dat400$c/dat400$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS5, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 108)
j <- j[103:107]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C400,index.return = TRUE)$ix, 108)
l <- l[103:107]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C400 == C400[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C400 == C400[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C400 == C400[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C400 == C400[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C400 == C400[l[5]], arr.ind = TRUE)[1,2]

nn5_400 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn5_400 = nn5_400+1
  }else{nn5_400 = nn5_400+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn5_400 = nn5_400/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 1 for DS5 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
#k=25
plot(dat25$d, dat25$c25, main = "DS5, K=25")
#k=100
plot(dat100$d, dat100$c100, main = "DS5, K=100")
#k=225
plot(dat225$d, dat225$c225, main = "DS5, K=225")
#k=400
plot(dat400$d, dat400$c400, main = "DS5, K=400")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 2 for DS5 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
require(gridExtra)
#k=25
p1 <- ggplot(dat25, aes(x=r)) + geom_histogram()+ggtitle("DS5 K=25")
#k=100 
p2 <- ggplot(dat100, aes(x=r)) + geom_histogram()+ggtitle("DS5 K=100")
#k=225 
p3 <- ggplot(dat225, aes(x=r)) + geom_histogram() +ggtitle("DS5 K=225")
#k=400
p4 <- ggplot(dat400, aes(x=r)) + geom_histogram() +ggtitle("DS5 K=400")

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#      Running the FJLT Algo    #~~~~~~~~~~~~#
#~~~~~~~~~~~~#         on data set 6       #~~~~~~~~~~~~#
#~~~~~~~~~~~~#  For k = 25, 100, 225, 400   #~~~~~~~~~~~~#
#~~~~~~~~~~~~# and recording the  run time  #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

time25 <- system.time(Y_s6_25 <- FJLT(data6, 25)) 
time100 <- system.time(Y_s6_100 <- FJLT(data6, 100))
time225 <- system.time(Y_s6_225 <- FJLT(data6, 225))
time400 <- system.time(Y_s6_400 <- FJLT(data6, 400))

runtime6  = matrix(nrow = 4, ncol = 1)
runtime6 <- as.data.frame(runtime6)
runtime6 <- data.matrix(runtime6)

runtime6[1,1] <- time25[3]
runtime6[2,1] <- time100[3]
runtime6[3,1] <- time225[3]
runtime6[4,1] <- time400[3]








#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS6, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C25 = matrix(nrow = ncol(Y_s6_25), ncol = ncol(Y_s6_25))
C25 <- as.data.frame(C25)
C25  <- data.matrix(C25)

j=1
while(j <= ncol(Y_s6_25)){
  l <- 1
  while (l <= ncol(Y_s6_25)){
    i = 1:nrow(Y_s6_25)
    C25[j,l] <- sum((Y_s6_25[i,j]-Y_s6_25[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

##################
#   Creating D   #
##################
D = matrix(nrow = ncol(data6), ncol = ncol(data6))
D <- as.data.frame(D)
D <- data.matrix(D)


for(j in 1:ncol(data6)){
  l <- 1
  while (l <= ncol(data6)){
    i = 1:nrow(data6)
    D[j,l] <- sum((data6[i,j]-data6[i,l])^2)
    
    l=l+1
  }
}


dat25 <- data.frame(d=as.vector(D[upper.tri(D)]),
                    c25=as.vector(C25[upper.tri(C25)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS6, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat25$r <- dat25$c/dat25$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS6, k=25      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 69)
j <- j[64:68]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C25,index.return = TRUE)$ix, 69)
l <- l[64:68]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C25 == C25[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C25 == C25[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C25 == C25[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C25 == C25[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C25 == C25[l[5]], arr.ind = TRUE)[1,2]

nn6_25 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn6_25 = nn6_25+1
  }else{nn6_25 = nn6_25+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn6_25 = nn6_25/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS6, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C100 = matrix(nrow = ncol(Y_s6_100), ncol = ncol(Y_s6_100))
C100 <- as.data.frame(C100)
C100 <- data.matrix(C100)

j=1
while(j <= ncol(Y_s6_100)){
  l <- 1
  while (l <= ncol(Y_s6_100)){
    i = 1:nrow(Y_s6_100)
    C100[j,l] <- sum((Y_s6_100[i,j]-Y_s6_100[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat100 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c100=as.vector(C100[upper.tri(C100)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS6, k=100      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat100$r <- dat100$c/dat100$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS6, k=100     #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 69)
j <- j[64:68]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C100,index.return = TRUE)$ix, 69)
l <- l[64:68]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C100 == C100[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C100 == C100[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C100 == C100[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C100 == C100[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C100 == C100[l[5]], arr.ind = TRUE)[1,2]

nn6_100 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn6_100 = nn6_100+1
  }else{nn6_100 = nn6_100+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn6_100 = nn6_100/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS6, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C225 = matrix(nrow = ncol(Y_s6_225), ncol = ncol(Y_s6_225))
C225 <- as.data.frame(C225)
C225 <- data.matrix(C225)

j=1
while(j <= ncol(Y_s6_225)){
  l <- 1
  while (l <= ncol(Y_s6_225)){
    i = 1:nrow(Y_s6_225)
    C225[j,l] <- sum((Y_s6_225[i,j]-Y_s6_225[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat225 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c225=as.vector(C225[upper.tri(C225)]))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS6, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat225$r <- dat225$c/dat225$d


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS6, k=225      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 69)
j <- j[64:68]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C225,index.return = TRUE)$ix, 69)
l <- l[64:68]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C225 == C225[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C225 == C225[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C225 == C225[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C225 == C225[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C225 == C225[l[5]], arr.ind = TRUE)[1,2]

nn6_225 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn6_225 = nn6_225+1
  }else{nn6_225 = nn6_225+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn6_225 = nn6_225/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 1 for DS6, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##################
#   Creating C   #
##################

C400 = matrix(nrow = ncol(Y_s6_400), ncol = ncol(Y_s6_400))
C400 <- as.data.frame(C400)
C400 <- data.matrix(C400)

j=1
while(j <= ncol(Y_s6_400)){
  l <- 1
  while (l <= ncol(Y_s6_400)){
    i = 1:nrow(Y_s6_400)
    C400[j,l] <- sum((Y_s6_400[i,j]-Y_s6_400[i,l])^2)
    
    l=l+1
  }
  j=j+1
}

dat400 <- data.frame(d=as.vector(D[upper.tri(D)]),
                     c400=as.vector(C400[upper.tri(C400)]))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 2 for DS6, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

dat400$r <- dat400$c/dat400$d

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#    Step 3 for DS6, k=400      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

j <- head(sort(D,index.return = TRUE)$ix, 69)
j <- j[64:68]
j_1 <- as.array(0,0,0,0,0)
j_1[1] <- which(D == D[j[1]], arr.ind = TRUE)[1,2]  
j_1[2] <- which(D == D[j[2]], arr.ind = TRUE)[1,2]
j_1[3] <- which(D == D[j[3]], arr.ind = TRUE)[1,2]
j_1[4] <- which(D == D[j[4]], arr.ind = TRUE)[1,2]
j_1[5] <- which(D == D[j[5]], arr.ind = TRUE)[1,2]
l <- head(sort(C400,index.return = TRUE)$ix, 69)
l <- l[64:68]
l_1 <- as.array(0,0,0,0,0)
l_1[1] <- which(C400 == C400[l[1]], arr.ind = TRUE)[1,2]  
l_1[2] <- which(C400 == C400[l[2]], arr.ind = TRUE)[1,2]
l_1[3] <- which(C400 == C400[l[3]], arr.ind = TRUE)[1,2]
l_1[4] <- which(C400 == C400[l[4]], arr.ind = TRUE)[1,2]
l_1[5] <- which(C400 == C400[l[5]], arr.ind = TRUE)[1,2]

nn6_400 = 0
for(i in (1:length(l_1))){
  if(l_1[i] %in% j_1){
    nn6_400 = nn6_400+1
  }else{nn6_400 = nn6_400+0}
}
#the fraction of nearest neighbors that are correctly identified after dimension reduction
nn6_400 = nn6_400/5

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 1 for DS6 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
par(mar=c(1,1,1,1))
#k=25
plot(dat25$d, dat25$c25, main = "DS6, K=25")
#k=100
plot(dat100$d, dat100$c100, main = "DS6, K=100")
#k=225
plot(dat225$d, dat225$c225, main = "DS6, K=225")
#k=400
plot(dat400$d, dat400$c400, main = "DS6, K=400")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#Plots of Step 2 for DS4 all Ks#~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
par(mfrow=c(2,2))
require(gridExtra)
#k=25
p1 <- ggplot(dat25, aes(x=r)) + geom_histogram()+ggtitle("DS6 K=25")
#k=100 
p2 <- ggplot(dat100, aes(x=r)) + geom_histogram()+ggtitle("DS6 K=100")
#k=225 
p3 <- ggplot(dat225, aes(x=r)) + geom_histogram() +ggtitle("DS6 K=225")
#k=400
p4 <- ggplot(dat400, aes(x=r)) + geom_histogram() +ggtitle("DS6 K=400")

grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#   Plots of Step 3 for  all   #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

k <- c(25,100,225,400)

nn1 <- c(nn1_25, nn1_100, nn1_225, nn1_400) 
nn2 <- c(nn2_25, nn2_100, nn2_225, nn2_400) 
nn3 <- c(nn3_25, nn3_100, nn3_225, nn3_400) 
nn4 <- c(nn4_25, nn4_100, nn4_225, nn4_400) 
nn5 <- c(nn5_25, nn5_100, nn5_225, nn5_400) 
nn6 <- c(nn6_25, nn6_100, nn6_225, nn6_400) 

s3 <- as.data.frame(k)
s3$nn1 <- nn1
s3$nn2 <- nn2
s3$nn3 <- nn3
s3$nn4 <- nn4
s3$nn5 <- nn5
s3$nn6 <- nn6

matplot(x= s3[,1], y= as.matrix(s3[-1]), type='l', pch=1, 
        col= 2:7, xlab='k', ylab = 'nn')
legend("right", inset=.1, legend=c("DS1", "DS2", "DS3", "DS4", "DS5", "DS6"), 
       pch=2, col= 2:7, cex = 0.5, horiz=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~#     Plots of #2for  all      #~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

t <- as.data.frame(k)
t$DS1 <- runtime
t$DS2 <- runtime2
t$DS3 <- runtime3
t$DS4 <- runtime4
t$DS5 <- runtime5
t$DS6 <- runtime6

matplot(x= t[,1], y= as.matrix(t[-1]), type='l', pch=1, 
        col= 2:7, xlab='k', ylab = 'CPU Time')
legend("topleft", inset=.05, legend=c("DS1", "DS2", "DS3", "DS4", "DS5", "DS6"), 
       pch=1, col= 2:7, cex = 0.5, horiz=TRUE)
