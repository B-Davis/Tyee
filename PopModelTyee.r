source("code/functions.r")
##################
### Parameters ###
##################

n <- 1000 # rows per age class - different threads of stochasticity
z <- 100 # number of matrices - years into future
Sex <- .5 # proportion of females - natural population
YYS <- .18 # YY survival from age 0 to spawn
YYN <- 3000
# YYS0 <- 0
#
# m <- matrix(c(0,.425,0,1.5,0,.12,2.3,0,.015),ncol = 3) # Male survival matrix
# mM <- matrix(c(0,.425,0,0,0,.75,0,0,.035),ncol = 3) # Male survival matrix
# mF <- matrix(c(0,.425,0,3,0,.75,4.6,0,.035),ncol = 3) # Female
# mYY <- matrix(c(0,.18,0,0,0,.25,0,0,.09),ncol = 3)
# st <- c(3000,2000,1000) # Starting abundance - 0,1,ad
# stF <- matrix(rep(c(3000,2000,1000)/2,n,each = n),n,3)
# stM <- matrix(rep(c(3000,2000,1000)/2,n,each = n),n,3)
# stYY <- cbind(rbinom(n = n,size = YYN,prob = YYS),0,0) # Starting abundance - 0,1,ad
# # stYY <- cbind(rep(YYN,n),0,0) # Starting abundance - 0,1,ad
# ST <- cbind(stF,stM,stYY,NA,NA,NA,NA,NA)
S <- c(.44,.18,.125,.44,.18,.125,.18,.05,.01) # Survival: 0-1,1-ad,ad-ad for females, males and YYs
SYY <- c(.18,.18,.18)
st <- c(4020,2744,376) # Starting abundance - 0,1,ad
Fc <- c(0,2.5,7.2) # Fecundity
###################
### Suppression ###
###################
# sup <- c(.2,.2,.2)
sup <- c(.5,.5,.5) # vector == proprtion of fish suppressed, age 0 ,1, and 2
# sup <- c(.3,.5,.4) # vector == proprtion of fish suppressed, age 0 ,1, and 2
sdd <- .024 # standard deviation of suppression rates
# sapply(1 - sup,function(x) rbeta(n,beta.mom(x,sdd)[1],beta.mom(x,sdd)[2]))
CE <- function(x,top = .03,b) (1-exp(-b*x))*top # capture efficiency as a f() of abundance
SuppressionRateDecay <- .1 # try b = 0.1, .03,.008, .001, and .0005
curve(CE(x,top = .3,SuppressionRateDecay),0,1000,lwd = 3)





#
# f <- function(st,n,S,Fc){
# a <- rpois(n,Fc %*% st)
# b <- matrix(rbinom(n*3,st,prob = S),n,3,byrow = T)
# cbind(a,b[,1],apply(b[,2:3],1,sum))
# }
# f(st,n,S,Fc)

# st <- out[,1:9,1]
# YYS <- SYY

f <- function(st,n,S,Fc,YYS,sup,sdd,YYN){ # starting pop - n columns - survival rates - Fecundity - YY surv. rates - suppression rates - sd of sup. rates
  # Matrix Math w/ stochasticity & YY male surv from stock to spawn
  a <- sapply(1:n,function(x) rpois(1,Fc %*% st[x,1:3]))
  a2 <- t(sapply(1:nrow(st),function(x) matrix(rbinom(9,st[x,1:9],prob = S),1,9,byrow = T)))
  M <- cbind(a,a2[,1],apply(a2[,2:3],1,sum),0,a2[,4],apply(a2[,5:6],1,sum),a2[,7],a2[,8],a2[,9]) # apply(a2[,8:9],1,sum)
  colnames(M) <- c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad")
  # YY influence
  m <- apply(M[,c("YY0","YY1","YYad")],1,sum)
  ni <- apply(M[,c("M1","Mad")],1,sum)
  k <- apply(M[,c("F1","Fad")],1,sum)
  k[which(k > (m + ni))] <- (m + ni)[which(k > (m + ni))] # house keeping - says: if there are more females than males...adjust math so it works...i.e., pick n_males instead of n_females
  tmp <- rhyper(nn = n,m = m,n = ni,k = k)# Number of YY males contributing (hypergeometric = similar to binomial, but withM replacement)
  YYcont <- tmp/apply(M[,c("F1","Fad")],1,sum) # YY male contribution to eggs
  YYcont[!is.finite(YYcont)] <- 0# deal with nans (from 0/0)
  M[,"F0"] <- M[,"F0"] - round(YYcont*M[,"F0"]) # Take YY male progeny from F0
  M[,"M0"] <- M[,"M0"] + round(YYcont*M[,"F0"]) # Add to M0
  male0 <- round(M[,"F0"]/2) # natural Parents - male
  M[,"M0"] <- M[,"M0"] + male0
  M[,"F0"] <- M[,"F0"] - male0
  # Suppress
  tmp <- beta.mom(1 - t(CE(x = t(M[,1:6]),top =rep(sup,2),b = SuppressionRateDecay)),sdd)
  tmp[[1]][tmp[[1]] < 0] <- 0; tmp[[2]][tmp[[1]] < 0] <- 0
  tmp <- matrix(mapply(rbeta,"n" = 1,"shape1" = tmp[[1]],"shape2" = tmp[[2]]),n,6)
  M[,1:6] <- round(tmp * M[,1:6])
  M[,7:9] <- cbind(YYN,M[,7],apply(M[,8:9],1,sum))
  M[M < 0] <- 0
  M
}

out <- array(0,dim = c(n,9,z),dimnames = list(NULL,c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad"),paste0("T",1:z))) #,"Fn","Mn","YYn","N","Sex"
out[,,1] <- cbind(matrix(rep(st/2,n,each = n),n,3),matrix(rep(st/2,n,each = n),n,3),cbind(rep(YYN,n),0,0))
# Add Stoch
out[,c(1,4),1] <- apply(out[,c(1,4),1],2,rpois,n=n)
out[,1:6,1] <- round(f_suppress(out[,1:6,1],sdd,n = n,sup = sup) * out[,1:6,1])

for(k in 1:(z-1)) out[,,k+1] <- f(st = out[,,k],n = n,S = S,Fc = Fc, YYS = SYY,sup = sup,sdd = sdd,YYN = YYN)

# out[,,3] <- f(out[,,2],n,S,Fc,YYS)
# st <- out[,,2]

#############
### plots ###
#############
dim(out)
tots <- apply(out[,1:3,],c(1,3),sum)
erad <- sapply(1:nrow(tots),function(x) which(tots[x,] == 0)[1L]) - 1
any(is.na(erad)) # if true the z not big enough
hist(erad,breaks = 100)
median(erad)
min(erad)
# barplot(tots,beside = T)
asd <- apply(tots,2,median)
plot(1:z,asd,type = "l",xlim = c(1,25),ylim = c(0,2000),xaxs = "i")
invisible(sapply(1:n,function(x) lines(1:z,tots[x,])))
lines(1:z,asd,col = "red",lwd = 2)
abline(h = 0,lty = 2)

########################
### Sensitivity Anal ###
########################
#### One way - Suppression ###
tmp1 <- seq(.05,.95,length.out = 60)
tmp <- rep(tmp1,each = 3)
tmp <- split(tmp,rep(1:length(tmp1),each = 3))
#
outSup <- array(0,dim = c(n,9,z,length(tmp)),dimnames = list(NULL,c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad"),paste0("T",1:z),tmp1)) #,"Fn","Mn","YYn","N","Sex"
dim(outSup);dimnames(outSup)

outSup[,,1,] <- cbind(matrix(rep(st/2,n,each = n),n,3),matrix(rep(st/2,n,each = n),n,3),cbind(rep(YYN,n),0,0))
for(ii in seq_along(tmp)) outSup[,1:6,1,ii] <-round(f_suppress(mat = outSup[,1:6,1,ii],sdd=sdd,n = n,sup = tmp[[ii]]) * outSup[,1:6,1,ii])
for(ii in seq_along(tmp)){
for(k in 1:(z-1)) outSup[,,k+1,ii] <- f(st = outSup[,,k,ii],n = n,S = S,Fc = Fc, YYS = SYY,sup = tmp[[ii]],sdd = sdd,YYN = YYN)
}
outSup[,,1,1]
outSup[,,100,19]
outSup[,,100,]
# save(outSup,file = "code/objects/outSup.Rdata")
load("code/objects/outSup.Rdata")
# outSupput
tots <- apply(outSup[,1:3,,],c(1,3,4),sum)
tots2 <- apply(tots,c(1,3),function(x) which(x == 0)[1L]-1)
aa <- apply(tots2,2,median)
range(aa)
plot(as.numeric(names(aa)),aa,xaxt = "n",xlim = c(0,1),type = "l",las = 2)
axis(1,seq(0,1,.1))

aa
hist(tots2[,])
dim(tots2)

#### One-way Stocking ###
n <- 500
z <- 50
tmp <- round(seq(500,9000,length.out = 60))
#
outStock <- array(0,dim = c(n,9,z,length(tmp)),dimnames = list(NULL,c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad"),paste0("T",1:z),tmp)) #,"Fn","Mn","YYn","N","Sex"
dim(out);dimnames(out)
for(ii in seq_along(tmp)) outStock[,,1,ii] <- cbind(matrix(rep(st/2,n,each = n),n,3),matrix(rep(st/2,n,each = n),n,3),cbind(rep(tmp[ii],n),0,0))
for(ii in seq_along(tmp)) outStock[,1:6,1,ii] <-round(f_suppress(mat = outStock[,1:6,1,ii],sdd=sdd,n = n,sup = sup) * outStock[,1:6,1,ii])
for(ii in 1:length(tmp)){
  for(k in 1:(z-1)) outStock[,,k+1,ii] <- f(st = outStock[,,k,ii],n = n,S = S,Fc = Fc, YYS = SYY,sup = sup,sdd = sdd,YYN = tmp[ii])
}
outStock[,,2,2]
outStock[1:10,,100,66]
# save(outStock,file = "code/objects/outStock.Rdata")
load("code/objects/outStock.Rdata")
# output
tots <- apply(outStock[,1:3,,],c(1,3,4),sum)
tots2 <- apply(tots,c(1,3),function(x) which(x == 0)[1L]-1)
aa <- apply(tots2,2,median,na.rm = T)
range(aa)
plot(as.numeric(names(aa)),aa,xaxt = "n",xlim = c(0,7000),type = "l",las = 2)
axis(1,seq(0,tmp[length(tmp)],1000))

#### One-way Fecundity Age 1###
n <- 10
# tmp <- lapply(1:60,function(x) cbind(0,seq(1.1,5,length.out = 60),seq(3.6,14.2,length.out = 60))[x,])
tmp <- lapply(1:60,function(x) cbind(0,seq(1.1,5,length.out = 60),7.2)[x,])
tmp <- lapply(tmp,round,2)

outFec <- array(0,dim = c(n,9,z,length(tmp)),dimnames = list(NULL,c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad"),paste0("T",1:z),tmp)) #,"Fn","Mn","YYn","N","Sex"
outFec[,,1,] <- cbind(matrix(rep(st/2,n,each = n),n,3),matrix(rep(st/2,n,each = n),n,3),cbind(rep(YYN,n),0,0))
# suppress first matrices
for(ii in seq_along(tmp)) outFec[,1:6,1,ii] <-round(f_suppress(mat = outFec[,1:6,1,ii],sdd=sdd,n = n,sup = sup) * outFec[,1:6,1,ii])
# Run it
for(ii in 1:length(tmp)){
  for(k in 1:(z-1)) outFec[,,k+1,ii] <- f(st = outFec[,,k,ii],n = n,S = S,Fc = tmp[[ii]], YYS = SYY,sup = sup,sdd = sdd,YYN = YYN)
}

tots <- apply(outFec[,1:3,,],c(1,3,4),sum)
tots2 <- apply(tots,c(1,3),function(x) which(x == 0)[1L]-1)
aa <- apply(tots2,2,median,na.rm = T)
range(aa)
plot(1:60,aa,xaxt = "n",xlim = c(0,7000),type = "l",las = 2)

dim(outFec)
outFec[,,10,60]

#### One-way Fecundity Age 2 ###
n <- 10
# tmp <- lapply(1:60,function(x) cbind(0,seq(1.1,5,length.out = 60),seq(3.6,14.2,length.out = 60))[x,])
tmp <- lapply(1:60,function(x) cbind(0,2.5,seq(3.6,14.2,length.out = 60))[x,])
tmp <- lapply(tmp,round,2)

outFec2 <- array(0,dim = c(n,9,z,length(tmp)),dimnames = list(NULL,c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad"),paste0("T",1:z),tmp)) #,"Fn","Mn","YYn","N","Sex"
outFec2[,,1,] <- cbind(matrix(rep(st/2,n,each = n),n,3),matrix(rep(st/2,n,each = n),n,3),cbind(rep(YYN,n),0,0))
# suppress first matrices
for(ii in seq_along(tmp)) outFec2[,1:6,1,ii] <-round(f_suppress(mat = outFec2[,1:6,1,ii],sdd=sdd,n = n,sup = sup) * outFec2[,1:6,1,ii])
# Run it
for(ii in 1:length(tmp)){
  for(k in 1:(z-1)) outFec2[,,k+1,ii] <- f(st = outFec2[,,k,ii],n = n,S = S,Fc = tmp[[ii]], YYS = SYY,sup = sup,sdd = sdd,YYN = YYN)
}

tots <- apply(outFec[,1:3,,],c(1,3,4),sum)
tots2 <- apply(tots,c(1,3),function(x) which(x == 0)[1L]-1)
aa <- apply(tots2,2,median,na.rm = T)
range(aa)
plot(1:60,aa,xaxt = "n",xlim = c(0,7000),type = "l",las = 2)

dim(outFec)
outFec[,,1,60]
#### two-way suppression ###
n <- 100
nn <- 10
tmpSu <- seq(.05,.95,length.out = nn)
tmpSt <- round(tmpSu*1e4)
# tmpSt <- round(seq(500,9000,length.out = nn))
asd <- expand.grid(tmpSt,tmpSu)
# asd <- asplit(asd,2)
nrow(asd)
# Starting point
out2way <- array(0,dim = c(n,9,z,nn**2),dimnames = list(NULL,c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad"),paste0("T",1:z),paste(asd[[1]],asd[[2]],sep = "-"))) #,"Fn","Mn","YYn","N","Sex"
# laydown first matrices with varying yys
for(ii in 1:nn**2) out2way[,,1,ii] <- cbind(matrix(rep(st/2,n,each = n),n,3),matrix(rep(st/2,n,each = n),n,3),cbind(rep(asd[ii,1],n),0,0))
# suppress T_0 with varying suppression rates
for(ii in 1:nn**2) out2way[,1:6,1,ii] <-round(f_suppress(mat = out2way[,1:6,1,ii],sdd=sdd,n = n,sup = asd[ii,2]) * out2way[,1:6,1,ii])
# Make the magic happen
for(ii in 1:nn**2){
  for(k in 1:(z-1)) out2way[,,k+1,ii] <- f(st = out2way[,,k,ii],n = n,S = S,Fc = Fc, YYS = SYY,sup = asd[ii,2],sdd = sdd,YYN = asd[ii,1])
}



out2way[,,5,8]
tots <- apply(out2way[,1:3,,],c(1,3,4),sum)
tots2 <- apply(tots,c(1,3),function(x) which(x == 0)[1L]-1)
aa <- apply(tots2,2,median,na.rm = T)
range(aa)
m <- matrix(names(aa),10,10)
m <- matrix(aa,10,10)
m[is.na(m)] <- "NA"
image(x = tmpSt,y=tmpSu,z = matrix(aa,10,10),las = 2,oldstyle = T,yaxt = "n",bty = "n")
lapply(1:nn, function(x) text(tmpSt,tmpSu[x],m[,x],cex = .8))
axis(2,seq(0,1,.2),seq(0,100,20),las = 2)
axis(1,tmpSt,tmpSt)

asd[1,tmpSu,m]
text(tmpSt,tmpSu,m)

axis(1,tmpSt)
contour(m,nlevels = 20,drawlabels = F)

plot(as.numeric(names(aa)),aa,xaxt = "n",xlim = c(0,7000),type = "l",las = 2)
axis(1,seq(0,tmp[length(tmp)],1000))
