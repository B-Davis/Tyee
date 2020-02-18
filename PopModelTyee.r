source("code/functions.r")
##################
### Parameters ###
##################

n <- 100 # rows per age class - different threads of stochasticity
z <- 100 # number of matrices - years into future
Sex <- .5 # proportion of females - natural population
YYS <- .18 # YY survival from age 0 to spawn
YYN <- 3000
# YYS0 <- 0
#
m <- matrix(c(0,.425,0,1.5,0,.75,2.3,0,.035),ncol = 3) # Male survival matrix
mM <- matrix(c(0,.425,0,0,0,.75,0,0,.035),ncol = 3) # Male survival matrix
mF <- matrix(c(0,.425,0,3,0,.75,4.6,0,.035),ncol = 3) # Female
mYY <- matrix(c(0,.18,0,0,0,.25,0,0,.09),ncol = 3)
st <- c(3000,2000,1000) # Starting abundance - 0,1,ad
# stF <- matrix(rep(c(3000,2000,1000)/2,n,each = n),n,3)
# stM <- matrix(rep(c(3000,2000,1000)/2,n,each = n),n,3)
# # stYY <- cbind(rbinom(n = n,size = YYN,prob = YYS),0,0) # Starting abundance - 0,1,ad
# stYY <- cbind(rep(YYN,n),0,0) # Starting abundance - 0,1,ad
ST <- cbind(stF,stM,stYY,NA,NA,NA,NA,NA)
S <- c(.425,.75,.035,.425,.75,.035,.18,.18,.18)
SYY <- c(.18,.18,.18)
Fc <- c(mF[1,])
###################
### Suppression ###
###################
# sup <- c(.2,.2,.2)
sup <- c(.5,.4,.3) # vector == proprtion of fish suppressed, age 0 ,1, and 2
# sup <- c(.3,.5,.4) # vector == proprtion of fish suppressed, age 0 ,1, and 2
sdd <- .024 # standard deviation of suppression rates
# sapply(1 - sup,function(x) rbeta(n,beta.mom(x,sdd)[1],beta.mom(x,sdd)[2]))
CE <- function(x,top = .03,b) (1-exp(-b*x))*top # capture efficiency as a f() of abundance
SuppressionRateDecay <- .5 # try b = 0.1, .03,.008, .001, and .0005
# curve(CE(x,top = .3,SuppressionRateDecay),0,1000,lwd = 3)




#
# f <- function(st,n,S,Fc){
# a <- rpois(n,Fc %*% st)
# b <- matrix(rbinom(n*3,st,prob = S),n,3,byrow = T)
# cbind(a,b[,1],apply(b[,2:3],1,sum))
# }
# f(st,n,S,Fc)

st <- out[,1:9,1]
YYS <- SYY

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
hist(erad)
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
# One way - by parapmer
# Suppression
tmp1 <- seq(.05,.95,.02)
tmp <- rep(sup1,each = 3)
tmp <- split(sup,rep(1:length(sup1),each = 3))
#
out <- array(0,dim = c(n,9,z,length(tmp)),dimnames = list(NULL,c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad"),paste0("T",1:z),tmp1)) #,"Fn","Mn","YYn","N","Sex"
dim(out);dimnames(out)

out[,,1,] <- cbind(matrix(rep(st/2,n,each = n),n,3),matrix(rep(st/2,n,each = n),n,3),cbind(rep(YYN,n),0,0))
for(ii in seq_along(sup)) out[,1:6,1,ii] <-round(f_suppress(mat = out[,1:6,1,ii],sdd=sdd,n = n,sup = tmp[[ii]]) * out[,1:6,1,ii])
for(ii in seq_along(sup)){
for(k in 1:(z-1)) out[,,k+1,ii] <- f(st = out[,,k,ii],n = n,S = S,Fc = Fc, YYS = SYY,sup = tmp[[ii]],sdd = sdd,YYN = YYN)
}
out[,,1,1]
out[,,100,19]
out[,,100,]

# output
tots <- apply(out[,1:3,,],c(1,3,4),sum)
tots2 <- apply(tots,c(1,3),function(x) which(x == 0)[1L]-1)
aa <- apply(tots2,2,median)
range(aa)
plot(as.numeric(names(aa)),aa,xaxt = "n",xlim = c(0,1),type = "l",las = 2)
axis(1,seq(0,1,.1))

# Stocking
tmp <- seq(500,7000,100)
#
outStock <- array(0,dim = c(n,9,z,length(tmp)),dimnames = list(NULL,c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad"),paste0("T",1:z),tmp)) #,"Fn","Mn","YYn","N","Sex"
dim(out);dimnames(out)
for(ii in seq_along(tmp)) outStock[,,1,ii] <- cbind(matrix(rep(st/2,n,each = n),n,3),matrix(rep(st/2,n,each = n),n,3),cbind(rep(tmp[ii],n),0,0))
for(ii in seq_along(tmp)) outStock[,1:6,1,ii] <-round(f_suppress(mat = outStock[,1:6,1,ii],sdd=sdd,n = n,sup = sup) * outStock[,1:6,1,ii])
for(ii in 1:length(tmp)){
  for(k in 1:(z-1)) outStock[,,k+1,ii] <- f(st = outStock[,,k,ii],n = n,S = S,Fc = Fc, YYS = SYY,sup = sup,sdd = sdd,YYN = tmp[ii])
}

outStock[,,100,19]
# output
tots <- apply(outStock[,1:3,,],c(1,3,4),sum)
tots2 <- apply(tots,c(1,3),function(x) which(x == 0)[1L]-1)
aa <- apply(tots2,2,median)
range(aa)
plot(as.numeric(names(aa)),aa,xaxt = "n",xlim = c(0,1),type = "l",las = 2)
axis(1,seq(0,1,.1))
# two-way suppression -- later
asd <- seq(.1,.9,.1)
asd <- expand.grid(asd,asd,asd)
asd <- asplit(asd,1)
nrow(asd)
# Starting point
out <- array(0,dim = c(n,9,z,nrow(asd)),dimnames = list(NULL,c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad"),paste0("T",1:z))) #,"Fn","Mn","YYn","N","Sex"
out[,,1,] <- cbind(stF,stM,stYY)
for(i in 1:nrow(asd)) out[,1:6,1,i] <-round(f_suppress(mat = out[,1:6,1,i],sdd=sdd,n = n,sup = as.vector(asd[[i]])) * out[,1:6,1,i])
