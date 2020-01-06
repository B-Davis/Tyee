source("code/functions.r")
##################
### Parameters ###
##################
m <- matrix(c(0,.425,0,1.5,0,.75,2.3,0,.035),ncol = 3)
st <- c(3000,2000,1000) # Starting abundance - 0,1,ad
n <- 50 # rows per age class - different threads of stochasticity
z <- 100 # number of matrices - years into future
###################
### Suppression ###
###################
sup <- c(.3,.5,.4) # vector == proprtion of fish suppressed, age 0 ,1, and 2
sdd <- .024 # standard deviation of suppression rates
# sapply(1 - sup,function(x) rbeta(n,beta.mom(x,sdd)[1],beta.mom(x,sdd)[2]))
CE <- function(x,top,b) (1-exp(-b*x))*top # capture efficiency as a f() of abundance
SuppressionRateDecay <- .0005 # try b = .003,.008, .001, and .0005
curve(CE(x,top = .3,SuppressionRateDecay),0,6000,lwd = 3)
################
### YY Males ###
################
# Sex Ratio
Sex <- .5 # proportion of females - natural population
#
Nstock0 <- 1000
Nstock1 <- 3000 # stocking YY males age 1
Nstock2 <- 2000
# Kennedy et al.
YYS <- .18 # YYmale survival
YYfitness <- .58 # spawning Fitness

# messing
Nstock0*YYS



#############
### Model ###
#############
# # for diagnoses
mat <- m
vec <- matrix(st,1,3)
# vec <- structure(c(12, 5, 16, 12, 2, 11, 6, 6, 20, 7, 2, 1, 3, 1, 1, 
#                    3, 2, 1, 4, 4, 2, 0, 0, 1, 0, 2, 0, 1, 1, 1), .Dim = c(10L, 3L
#                    ), .Dimnames = list(NULL, c("0", "1", "ad")))
dsts <- list("rpois","rpois","rbinom","rbinom","rbinom")
g <- function(n,vec,mat) {
  args <- list(list("n" = n,"lambda" = vec[,2]*mat[1,2]),
               list("n" = n,"lambda" = vec[,3]*mat[1,3]),
               list("n" = n,"size" = vec[,1],"prob" = mat[2,1] ),
               list("n" = n,"size" = vec[,2],"prob" = mat[3,2] ),
               list("n" = n,"size" = vec[,3],"prob" = mat[3,3] ))
  tmp <- mapply(dsts,args,FUN = do.call) # columns = "age" 0 (from 1s), 0 (from ADs), 1s, ADs (2s), ADs
  tmp <- cbind(sapply(list(1:2,4:5),function(x) apply(tmp[,x],1,sum)),tmp[,3])[,c(1,3,2)]
  
  tmp2 <- 1-sapply(1:3,function(x) CE(tmp[,x],top = sup[x],b = SuppressionRateDecay)) # Capture efficiency base on curve
  a <- apply(tmp2,c(1,2),function(x) beta.mom(x,sd = sdd)[1])
  b <- apply(tmp2,c(1,2),function(x) beta.mom(x,sd = sdd)[2])
  tmp2 <- matrix(rbeta(tmp2,a,b),n,3)
  # round(tmp %*% diag(sup)) # removal
  round(tmp*tmp2)
}
out <- array(dim = c(n,3,z),dimnames = list(NULL,c("0","1","ad"),paste0("T",1:z)))
out[,,1] <- g(n,matrix(st,1,3),m)

for(i in 1:(z-1)) {
  out[,,i+1] <- g(n,out[,,i],m)
}
#############
### plots ###
#############
# by age
age <- 0
plot(1:z,out[1,age + 1,],type = "l",xlim = c(1,50),ann = F,xaxs = "i")
sapply(2:n,function(x) lines(1:100,out[x,age + 1,]))
lines(apply(out[,age + 1,],2,mean),col = "red")
# all
tmp <- apply(out,c(1,3),sum)
plot(1:z,tmp[1,],type = "l",xlim = c(0,100),xaxs = "i")
invisible(sapply(2:n,function(x) lines(1:100,tmp[x,])))
lines(apply(tmp,2,mean),col = "red")
#
apply(tmp,2,median)




########################################### Fecundity worksheet
# Scratch
1.5 * 600
# N * SexR * F * Sur
600*.5*200*.015 # if F = 1.5 thrn s0 == .015 when F = 200 ---------- Need to add another parameter to account for 1 year olds not as fecund or mature as adults.
1.5/(600/200*2) # do some algebra and figure out formula to get 1.5
600 = 2*200*1.5*x
x = 2*200*1.5


# parameters for F = 200
.023 # Egg/age0 Survival
0.6521739 # penalty Fec for being 1


# Proof
Sex <- .5
N <- 3000
Fec <- 200
S0 <- .023
Penal <- 0.6521739

N*2.3
N*1.5
N*Sex*Fec*S0
N*Sex*Fec*S0*Penal
(N*Sex*Fec*S0)/N
(N*Sex*Fec*S0*Penal)/N
#
