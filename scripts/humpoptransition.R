## Global human population switched from density enhancement to limitation in the 1950s
## Corey Bradshaw
## Flinders University 
## February 2024

# required R libraries
library(plotrix)
library(boot)
library(tmvnsim)
library(wCorr)
library(truncnorm)
library(orcutt)
library(lmtest)

# functions
# Hildreth-Lu (named for Clifford Hildreth and John Y. Lu, is a method for adjusting
# a linear model in response to the presence of serial correlation in the error term)
hildreth.lu.func <- function(r, model){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- 2:n
  y <- y[t]-r*y[t-1]
  x <- x[t]-r*x[t-1]
  
  return(lm(y~x))
}

hildreth.lu.order.func <- function(r, model, order){
  x <- model.matrix(model)[,-1]
  y <- model.response(model.frame(model))
  n <- length(y)
  t <- (order+1):n
  y <- y[t]-r*y[-c((n-(order-1)):n)]
  x <- x[t]-r*x[-c((n-(order-1)):n)]
  
  return(lm(y~x))
}

# geometric mean
gmMean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

## sample avoiding consecutive values of order o
NonSeqSample <- function(x, size, ord=1, replace) {
  # x = vector to sample
  # size = size of resample
  # ord = order to consider (1 = no consecutive, 2 = no 2nd-order consecutive ...)
  # replace = with (TRUE) or without (FALSE) replacement
  vsmp <- sort(sample(x=x, size=size, replace=F))
  diff.vsmp <- diff(vsmp)
  vsmp <- vsmp[-(which(diff.vsmp <= ord)+1)]
  new.size <- size - length(vsmp)
  
  while(new.size >= 1) {
    vsmp.add <- sort(sample(x=x, size=new.size+1, replace=F))
    vsmp <- sort(c(vsmp, vsmp.add))
    diff.vsmp <- diff(vsmp)
    vsmp <- vsmp[-(which(diff.vsmp <= ord)+1)]
    new.size <- size - length(vsmp)
  }
  return(vsmp)
}

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

## functions
delta.IC <- function(x) x - min(x) ## where x is a vector of an IC
weight.IC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.IC(AIC.vec); wAIC.vec <- weight.IC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

## Historical estimates of world population
## census.gov/data/tables/time-series/demo/international-programs/historical-est-worldpop.html
hpop <- read.csv('Npre1950.csv', sep=",", header=T)
head(hpop)
hpop$pop.md <- apply(hpop[,c(2,3)], MARGIN=1, median, na.rm=T)
head(hpop)
hpop.out <- data.frame(hpop[,1], hpop[,4], hpop[,3], hpop[,2])
colnames(hpop.out) <- c("year", "popMD", "popUP", "popLO")
head(hpop.out)

## UN population
## data.un.org/Data.aspx?d=POP&f=tableCode%3a1
UNpop <- read.csv('UNpop.csv', sep=",", header=T)
head(UNpop)
UNtot <- subset(UNpop, Area=="Total" & Sex=="Both Sexes")
head(UNtot)
UNtot2022 <- subset(UNtot, Year==2022)
head(UNtot2022)
sum(UNtot2022$Value)/10^9

popdat <- read.csv("worldpophist.csv")
head(popdat)
popdat$r <- c(NA, log(popdat$pop[2:length(popdat$pop)] / popdat$pop[1:(length(popdat$pop)-1)]))
head(popdat)
popdat$rpcap <- popdat$r/popdat$pop
head(popdat)

par(mfrow=c(1,3))
plot(popdat$year, popdat$pop, type="l")
plot(popdat$year, popdat$r, type="l")
plot(popdat$year, popdat$rpcap, type="l")
par(mfrow=c(1,1))

## interpolate yearly
yrintp <- seq(-10000, 2021, 1)
popintp <- approx(popdat$year, popdat$pop, xout = yrintp)

popdatintp <- data.frame(yrintp, popintp$y)
colnames(popdatintp) <- c("year", "pop")
popdatintp$r <- c(NA, log(popdatintp$pop[2:length(popdatintp$pop)] / popdatintp$pop[1:(length(popdatintp$pop)-1)]))
popdatintp$rpcap <- popdatintp$r/popdatintp$pop
head(popdatintp)

popdatintp$popdiff <- c(NA,diff(popdatintp$pop))
plot(popdatintp$year[2:(length(popdatintp$year))], (diff(popdatintp$pop)), type="l")

# to avoid pseudo-replication of interpolated r values, take mean r from wide (> annual) intervals
intervals.vec <- diff(popdat$year)
len.gtOne <- length(which(intervals.vec > 1))

mean.r <- rep(NA,len.gtOne)
for (i in 1:len.gtOne) {
  mean.r[i] <- mean(popdatintp$r[1:intervals.vec[i]], na.rm=T)
}

popdat$mean.r <- NA
popdat$mean.r[2:dim(popdat)[1]] <- c(mean.r, popdat$r[(length(mean.r)+2):dim(popdat)[1]])

mean(popdat$mean.r[which(popdat$year >= -8000 & popdat$year <= 0)])

plot(popdat$pop[1:(dim(popdat)[1])-1], popdat$mean.r[2:dim(popdat)[1]], pch=19, xlab="Nt", ylab="r", xlim=c(0,14e9))
popdat.phase1 <- subset(popdat, year < 1951)
points(popdat.phase1$pop[1:(dim(popdat.phase1)[1])-1], popdat.phase1$mean.r[2:dim(popdat.phase1)[1]], col="red", pch=19)
rN.phase1 <- data.frame(popdat.phase1$pop[1:(dim(popdat.phase1)[1])-1], popdat.phase1$mean.r[2:dim(popdat.phase1)[1]])
colnames(rN.phase1) <- c("Nt","r")
st.phase2 <- which(popdat$pop == popdat.phase1$pop[dim(popdat.phase1)[1]]) + 1
rN.phase2 <- data.frame(c(popdat.phase1$pop[dim(popdat.phase1)[1]], popdat$pop[st.phase2:(dim(popdat)[1]-1)]), popdat$mean.r[st.phase2:dim(popdat)[1]]) 
colnames(rN.phase2) <- c("Nt","r")

# Ricker logistic model
plot(rN.phase2$Nt, rN.phase2$r, xlim=c(min(rN.phase1$Nt, na.rm=T),14e9), ylim=c(min(rN.phase1$r), max(rN.phase2$r)) , col="black", xlab="Nt", ylab="r", pch=19)
points(rN.phase1$Nt, rN.phase1$r, pch=19, col="red")
fitRicker.phase1 <- lm(rN.phase1$r ~ rN.phase1$Nt)
ablineclip(fitRicker.phase1, x1=0, x2=3e9, col="red", lty=2, lwd=2)
fitRicker.phase2 <- lm(rN.phase2$r ~ rN.phase2$Nt)
abline(fitRicker.phase2, col="black", lty=2, lwd=2, xpd=F)
abline(h=0,col="black", lty=3)
k2Ricker <- as.numeric(-coef(fitRicker.phase2)[1]/coef(fitRicker.phase2)[2])
abline(v=k2Ricker, lty=3)
round(k2Ricker/10^9, 4)

range(popdatintp$pop)
range(popdatintp$r, na.rm=T)
range(popdatintp$rpcap,na.rm=T)
range(popdat$rpcap[which(popdat$year > 1899)],na.rm=T)


sub3br.13 <- subset(popdatintp, pop <= 3e009 & r < 0.013)
range(sub3br.13$year)

# remove 1950 (post-war anomaly)
sub3br.13p <- sub3br.13[-which(sub3br.13$year > 1941 & sub3br.13$year < 1951), ]
tail(sub3br.13p)
range(sub3br.13p$year)

plot((popdat$pop[which(popdat$year >= 1900)]), popdat$r[which(popdat$year >= 1900)], pch=19, xlab="Nt", ylab="r")
sub3b <- subset(popdat, pop >= 3e09)
plot((sub3b$pop[which(sub3b$year >= 1900)]), sub3b$r[which(sub3b$year >= 1900)], pch=19, xlab="Nt", ylab="r")
plot(sub3b$pop[1:(length(sub3b$pop)-1)], sub3b$r[2:length(sub3b$pop)], pch=19, xlab="Nt", ylab="r")
fitlin <- lm(sub3b$r[2:length(sub3b$pop)] ~ sub3b$pop[1:(length(sub3b$pop)-1)])
abline(fitlin, lty=2, col="red")
k <- as.numeric(-coef(fitlin)[1]/coef(fitlin)[2])
round(k/10^9, 4)


# N children 0-1 1950-2021 (United Nations Population Division)
N01obs <- read.csv("Nchild0-1.csv", header=T)
head(N01obs)

# calculate time to stability
# linear
pastC2 <- sqrt((8.04531145 - 3.12668672)^2 + (0.02185709 - 0.00876466)^2)
pastC2rate <- pastC2 / (2022-1961)

futC2 <- sqrt((11.89 - 8.04531145)^2 + (0.00876466 - 0)^2)
stabilityT <- round(futC2/pastC2rate, 0)
2022 + stabilityT

futC2 <- sqrt((11.55 - 8.04531145)^2 + (0.00876466 - 0)^2)
stabilityT <- round(futC2/pastC2rate, 0)
2022 + stabilityT

futC2 <- sqrt((12.26 - 8.04531145)^2 + (0.00876466 - 0)^2)
stabilityT <- round(futC2/pastC2rate, 0)
2022 + stabilityT

# log
pastC2 <- sqrt((2.08508949 - 1.13997389)^2 + (0.02185709 - 0.00876466)^2)
pastC2rate <- pastC2 / (2022-1961)

futC2 <- sqrt((2.888 - 2.08508949)^2 + (0.00876466 - 0)^2)
stabilityT <- round(futC2/pastC2rate, 0)
2022 + stabilityT

# predict N range at max depensation
intpre1950 <- 0.0008685
intpre1950lo <- -0.0002628
intpre1950up <- 0.002000
slppre1950 <- 0.003508
slppre1950lo <- 0.002791
slppre1950up <- 0.004225

ypre1950mx <- 0.0106954
xpre1950mx <- (ypre1950mx - intpre1950)/slppre1950
xpre1950mx

x.vec <- seq(1, 3, 0.01)
iter <- 10000
y.vec.mat <- matrix(data=NA, nrow=iter, ncol=length(x.vec))
for (i in 1:iter) {
  y.vec.mat[i,] <- runif(1, intpre1950lo, intpre1950up) + runif(1, slppre1950lo, slppre1950up) * x.vec
}
y.vec.lo <- apply(y.vec.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
y.vec.up <- apply(y.vec.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)
y.vec.md <- apply(y.vec.mat, MARGIN=2, median, na.rm=T)

plot(x.vec, y.vec.lo, type="l", lty=2, pch=19, col="red", ylim=c(min(y.vec.lo), max(y.vec.up)))
lines(x.vec, y.vec.up, type="l", lty=2, col="red")
lines(x.vec, y.vec.md, type="l", lty=2, col="black")
dat.out <- data.frame(x.vec, y.vec.md, y.vec.up, y.vec.lo)
colnames(dat.out) <- c("x", "md", "up", "lo")
dat.out

x.vec[which.min(abs(y.vec.md - ypre1950mx))]
x.vec[which.min(abs(y.vec.up - ypre1950mx))]
x.vec[which.min(abs(y.vec.lo - ypre1950mx))]


## population x temperature anomaly (metoffice.gov.uk/hadobs/hadcrut5/data/HadCRUT.5.0.2.0/download.html)
popXta <- read.csv("popXtempanom.csv", header=T)
head(popXta)
popXta.pre1950 <- subset(popXta, year < 1950)
popXta.1950.1961 <- subset(popXta, year >= 1950 & year <= 1961)
popXta.post1961 <- subset(popXta, year > 1961)

dat.use <- popXta.pre1950
#dat.use <- popXta.1950.1961
#dat.use <- popXta.post1961

iter <- 10000
itdiv <- iter/10

R2.vec <- ER.vec <- pHU.vec <- R2HU.vec <- rep(NA, iter)
for (i in 1:iter) {
  ta.samp <- runif(dim(dat.use)[1], dat.use$anomLO, dat.use$anomUP)
  R2.vec[i] <- linreg.ER(dat.use$pop, ta.samp)[2]
  ER.vec[i] <- linreg.ER(dat.use$pop, ta.samp)[1]
  
  popXta.fit <- lm(ta.samp ~ dat.use$pop)
  popXta.fit.orc <- cochrane.orcutt(popXta.fit, convergence = 5, max.iter=1000)
  popXta.fit.hl <- hildreth.lu.order.func(popXta.fit.orc$rho, popXta.fit, order=1)
  popXta.fit.summ <- summary(popXta.fit.hl)
  pHU.vec[i] <- popXta.fit.summ$coefficients[8]
  R2HU.vec[i] <- popXta.fit.summ$adj.r.squared
  
  if (i %% itdiv==0) print(i)

} # end i

R2lo <- quantile(R2.vec, probs=0.025, na.rm=T)
R2up <- quantile(R2.vec, probs=0.975, na.rm=T)
ERlo <- quantile(ER.vec, probs=0.025, na.rm=T)
ERup <- quantile(ER.vec, probs=0.975, na.rm=T)

print(c(R2lo, R2up))
print(c(ERlo, ERup))

pHUlo <- quantile(pHU.vec, probs=0.025, na.rm=T)
pHUup <- quantile(pHU.vec, probs=0.975, na.rm=T)
R2HUlo <- quantile(R2HU.vec, probs=0.025, na.rm=T)
R2HUup <- quantile(R2HU.vec, probs=0.975, na.rm=T)

print(c(R2HUlo, R2HUup))
print(c(pHUlo, pHUup))


## regional Ricker logistic fits
popreg <- read.csv("popregions.csv", header=T)
head(popreg)

## scale
# SUB-SAHARAN AFRICA
SSA <- popreg[,c(1,2)]
SSA$Nsc <- as.numeric(scale(SSA[,2], scale=T, center=F))
SSA$r <- c(log(SSA$Nsc[2:dim(popreg)[1]] / SSA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(SSA)
SSA.use <- SSA[which(SSA$year >= 2010),]
plot(SSA.use$Nsc, SSA.use$r, xlab="N", ylab="r", pch=19)
SSA.rick <- lm(SSA.use$r ~ SSA.use$Nsc)
abline(SSA.rick, lty=2, col="red")
summary(SSA.rick)

# NORTH AFRICA AND WESTERN ASIA
NAWA <- popreg[,c(1,4)]
NAWA$Nsc <- as.numeric(scale(NAWA[,2], scale=T, center=F))
NAWA$r <- c(log(NAWA$Nsc[2:dim(popreg)[1]] / NAWA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(NAWA)
NAWA.use <- NAWA[which(NAWA$year >= 1975),]
plot(NAWA.use$Nsc, NAWA.use$r, xlab="N", ylab="r", pch=19)
NAWA.rick <- lm(NAWA.use$r ~ NAWA.use$Nsc)
abline(NAWA.rick, lty=2, col="red")
summary(NAWA.rick)

# CENTRAL AND SOUTH ASIA
head(popreg)
CSA <- popreg[,c(1,6)]
CSA$Nsc <- as.numeric(scale(CSA[,2], scale=T, center=F))
CSA$r <- c(log(CSA$Nsc[2:dim(popreg)[1]] / CSA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(CSA)
CSA.use <- CSA[which(CSA$year >= 1983),]
plot(CSA.use$Nsc, CSA.use$r, xlab="N", ylab="r", pch=19)
CSA.rick <- lm(CSA.use$r ~ CSA.use$Nsc)
abline(CSA.rick, lty=2, col="red")
summary(CSA.rick)

# EAST AND SOUTH-EASTERN ASIA
head(popreg)
ESEA <- popreg[,c(1,8)]
ESEA$Nsc <- as.numeric(scale(ESEA[,2], scale=T, center=F))
ESEA$r <- c(log(ESEA$Nsc[2:dim(popreg)[1]] / ESEA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(ESEA)
ESEA.use <- ESEA[which(ESEA$year >= 1963),]
plot(ESEA.use$Nsc, ESEA.use$r, xlab="N", ylab="r", pch=19)
ESEA.rick <- lm(ESEA.use$r ~ ESEA.use$Nsc)
abline(ESEA.rick, lty=2, col="red")
summary(ESEA.rick)

# EAST AND SOUTH-EASTERN ASIA EXCLUDING CHINA
china <- read.csv("china1950-2021.csv", header=T)
head(china)
ESEAexclCHN <- ESEA
head(ESEAexclCHN)
ESEAexclCHN$ESEA <- ESEA$ESEA - china$N
ESEAexclCHN$Nsc <- as.numeric(scale(ESEAexclCHN[,2], scale=T, center=F))
ESEAexclCHN$r <- c(log(ESEAexclCHN$Nsc[2:dim(popreg)[1]] / ESEAexclCHN$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(ESEAexclCHN)
plot(ESEAexclCHN$Nsc, ESEAexclCHN$r, xlab="N", ylab="r", pch=19)
ESEAexclCHN.use <- ESEAexclCHN[which(ESEAexclCHN$year >= 1958),]
plot(ESEAexclCHN.use$Nsc, ESEAexclCHN.use$r, xlab="N", ylab="r", pch=19)
ESEAexclCHN.rick <- lm(ESEAexclCHN.use$r ~ ESEAexclCHN.use$Nsc)
abline(ESEAexclCHN.rick, lty=2, col="red")
summary(ESEAexclCHN.rick)

# CHINA
CHN <- china
head(china)
CHN$Nsc <- as.numeric(scale(CHN[,2], scale=T, center=F))
CHN$r <- c(log(CHN$Nsc[2:dim(popreg)[1]] / CHN$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(CHN)
plot(CHN$Nsc, CHN$r, xlab="N", ylab="r", pch=19)
CHN.use <- CHN[which(CHN$year >= 1962),]
plot(CHN.use$Nsc, CHN.use$r, xlab="N", ylab="r", pch=19)
CHN.rick <- lm(CHN.use$r ~ CHN.use$Nsc)
abline(CHN.rick, lty=2, col="red")
summary(CHN.rick)

# LATIN AMERICA AND CARIBBEAN
head(popreg)
LAC <- popreg[,c(1,10)]
LAC$Nsc <- as.numeric(scale(LAC[,2], scale=T, center=F))
LAC$r <- c(log(LAC$Nsc[2:dim(popreg)[1]] / LAC$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(LAC)
LAC.use <- LAC[which(LAC$year >= 1960),]
plot(LAC.use$Nsc, LAC.use$r, xlab="N", ylab="r", pch=19)
LAC.rick <- lm(LAC.use$r ~ LAC.use$Nsc)
abline(LAC.rick, lty=2, col="red")
summary(LAC.rick)

# OCEANIA
head(popreg)
OC <- popreg[,c(1,12)]
OC$Nsc <- as.numeric(scale(OC[,2], scale=T, center=F))
OC$r <- c(log(OC$Nsc[2:dim(popreg)[1]] / OC$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(OC)
OC.use <- OC[which(OC$year >= 1996),]
plot(OC.use$Nsc, OC.use$r, xlab="N", ylab="r", pch=19)
OC.rick <- lm(OC.use$r ~ OC.use$Nsc)
abline(OC.rick, lty=2, col="red")
summary(OC.rick)

# EUROPE AND NORTH AMERICA
head(popreg)
EUNA <- popreg[,c(1,14)]
EUNA$Nsc <- as.numeric(scale(EUNA[,2], scale=T, center=F))
EUNA$r <- c(log(EUNA$Nsc[2:dim(popreg)[1]] / EUNA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(EUNA)
EUNA.use <- EUNA[which(EUNA$year >= 1957),]
plot(EUNA.use$Nsc, EUNA.use$r, xlab="N", ylab="r", pch=19)
EUNA.rick <- lm(EUNA.use$r ~ EUNA.use$Nsc)
abline(EUNA.rick, lty=2, col="red")
summary(EUNA.rick)


ESEAzero <- c(2547, 2664)
EUNAzero <- c(1186, 1230)
LACzero <- c(821.7, 837.7)
CSAzero <- c(2772, 2848)
OCzero <- c(21.90, 23.04)
NAWAzero <- c(906.4, 1150)
SSAzero <- c(3726, 4984)

sum(ESEAzero[1], EUNAzero[1], LACzero[1], CSAzero[1], OCzero[1], NAWAzero[1], SSAzero[1])
sum(ESEAzero[2], EUNAzero[2], LACzero[2], CSAzero[2], OCzero[2], NAWAzero[2], SSAzero[2])
