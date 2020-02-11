source("code/functions.r")
##################
### Parameters ###
##################

n <- 10 # rows per age class - different threads of stochasticity
z <- 50 # number of matrices - years into future
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
stF <- matrix(rep(c(3000,2000,1000)/2,n,each = n),n,3)
stM <- matrix(rep(c(3000,2000,1000)/2,n,each = n),n,3)
stYY <- cbind(rbinom(n = n,size = YYN,prob = YYS),0,0) # Starting abundance - 0,1,ad
ST <- cbind(stF,stM,stYY,NA,NA,NA,NA,NA)
S <- c(.425,.75,.035,.425,.75,.035,1,1,1)
SYY <- c(.18,.18,.18,.18)
Fc <- c(mF[1,])

# Try this again with stoch atart
out <- array(0,dim = c(n,9,z),dimnames = list(NULL,c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad"),paste0("T",1:z))) #,"Fn","Mn","YYn","N","Sex"
out[,,1] <- cbind(stF,stM,stYY)
# Add Stoch
out[,c(1,4),1] <- apply(out[,c(1,4),1],2,rpois,n=n)
out[,1:6,1] <- round(f_suppress(out[,1:6,1],sdd,n = n,sup = sup) * out[,1:6,1])

cbind(rbinom(n = n,size = YYN,prob = YYS),0,0)
#
# f <- function(st,n,S,Fc){
# a <- rpois(n,Fc %*% st)
# b <- matrix(rbinom(n*3,st,prob = S),n,3,byrow = T)
# cbind(a,b[,1],apply(b[,2:3],1,sum))
# }
# f(st,n,S,Fc)

st <- out[,1:9,1]
YYS <- SYY

f <- function(st,n,S,Fc,YYS){
  # Matrix Math w/ stochasticity
  a <- sapply(1:n,function(x) rpois(1,Fc %*% st[x,1:3]))
  a2 <- t(sapply(1:nrow(st),function(x) matrix(rbinom(9,st[x,1:9],prob = S),1,9,byrow = T)))
  M <- cbind(a,a2[,1],apply(a2[,2:3],1,sum),0,a2[,4],apply(a2[,5:6],1,sum),a2[,7],0,0) # apply(a2[,8:9],1,sum)
  colnames(M) <- c("F0","F1","Fad","M0","M1","Mad","YY0","YY1","YYad")
  # YY influence
  m <- apply(M[,c("YY0","YY1","YYad")],1,sum)
  ni <- apply(M[,c("M1","Mad")],1,sum)
  k <- apply(M[,c("F1","Fad")],1,sum)
  k[which(k > (m + ni))] <- (m + ni)[which(k > (m + ni))] # house keeping - says: if there are more females than males...adjust math so it works...i.e., pick n_males instead of n_females
  tmp <- rhyper(nn = n,m = m,n = ni,k = k)# Number of YY males contributing (hypergeometric = similar to binomial, but withM replacement)
  YYcont <- tmp/apply(M[,c("F1","Fad")],1,sum) # YY male contribution to eggs
  YYcont[!is.finite(YYcont)] <- 0# deal with nans (from 0/0)
  # M[,"YY0",i+1] <- stYY[,1] # Add YY males
  M[,"F0"] <- M[,"F0"] - round(YYcont*M[,"F0"]) # Take YY male progeny from F0
  M[,"M0"] <- M[,"M0"] + round(YYcont*M[,"F0"]) # Add to M0
  male0 <- round(M[,"F0"]/2) # natural Parents - male
  M[,"M0"] <- M[,"M0"] + male0
  M[,"F0"] <- M[,"F0"] - male0
  # Suppress
  M[,1:6] <- round(f_suppress(M[,1:6],sdd,n = n,sup = sup) * M[,1:6])
  # Add YY males
  M[,"YY0"] <- YYN
  # YY surv to spawn
  a3 <- t(sapply(1:nrow(M),function(x) matrix(rbinom(4,c(YYN,M[x,7:9]),prob = YYS),1,4,byrow = T)))
  M[,7:9] <- cbind(a3[,1],a3[,2],apply(a3[,3:4],1,sum))
  M[M < 0] <- 0
  M
}

for(j in 1:(z-1)) out[,,j+1] <- f(st = out[,,j],n = n,S = S,Fc = Fc, YYS = SYY)

f(out[,,1],n,S,Fc,YYS)
