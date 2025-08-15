## Global human population ended self-facilitation in the 1950s
## Corey Bradshaw
## Flinders University 
## February 2024 / updated August 2025

# required R libraries
library(plotrix)
library(boot)
library(tmvnsim)
library(wCorr)
library(truncnorm)
library(orcutt)
library(lmtest)
library(performance)
library(sjPlot)
library(dismo)
library(gbm)

# source files
source("new_lmer_AIC_tables3.R")
source("r.squared.R")

# functions
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
## https://data.un.org/Data.aspx?d=POP&f=tableCode%3a1
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

## top 10 countries with highest fertility
# Niger
NER <- read.csv("NER.csv", header=T)
head(NER)
NER$Nsc <- as.numeric(scale(NER[,2], scale=T, center=F))
NER$r <- c(log(NER$Nsc[2:dim(popreg)[1]] / NER$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(NER)
plot(NER$Nsc, NER$r, xlab="N", ylab="r", pch=19)
NER.use <- NER[which(NER$year >= 2015),]
plot(NER.use$Nsc, NER.use$r, xlab="N", ylab="r", pch=19)
NER.rick <- lm(NER.use$r ~ NER.use$Nsc)
abline(NER.rick, lty=2, col="red")
summary(NER.rick)

# Democratic Republic of Congo
COD <- read.csv("COD.csv", header=T)
head(COD)
COD$Nsc <- as.numeric(scale(COD[,2], scale=T, center=F))
COD$r <- c(log(COD$Nsc[2:dim(popreg)[1]] / COD$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(COD)
plot(COD$Nsc, COD$r, xlab="N", ylab="r", pch=19)
COD.use <- COD[which(COD$year >= 2013),]
plot(COD.use$Nsc, COD.use$r, xlab="N", ylab="r", pch=19)
COD.rick <- lm(COD.use$r ~ COD.use$Nsc)
abline(COD.rick, lty=2, col="red")
summary(COD.rick)

# Mali
MLI <- read.csv("MLI.csv", header=T)
head(MLI)
MLI$Nsc <- as.numeric(scale(MLI[,2], scale=T, center=F))
MLI$r <- c(log(MLI$Nsc[2:dim(popreg)[1]] / MLI$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(MLI)
plot(MLI$Nsc, MLI$r, xlab="N", ylab="r", pch=19)
MLI.use <- MLI[which(MLI$year >= 2004),]
plot(MLI.use$Nsc, MLI.use$r, xlab="N", ylab="r", pch=19)
MLI.rick <- lm(MLI.use$r ~ MLI.use$Nsc)
abline(MLI.rick, lty=2, col="red")
summary(MLI.rick)

# Chad
TCD <- read.csv("TCD.csv", header=T)
head(TCD)
TCD$Nsc <- as.numeric(scale(TCD[,2], scale=T, center=F))
TCD$r <- c(log(TCD$Nsc[2:dim(popreg)[1]] / TCD$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(TCD)
plot(TCD$Nsc, TCD$r, xlab="N", ylab="r", pch=19)
TCD.use <- TCD[which(TCD$year >= 2003),]
plot(TCD.use$Nsc, TCD.use$r, xlab="N", ylab="r", pch=19)
TCD.rick <- lm(TCD.use$r ~ TCD.use$Nsc)
abline(TCD.rick, lty=2, col="red")
summary(TCD.rick)

# Angola
AGO <- read.csv("AGO.csv", header=T)
head(AGO)
AGO$Nsc <- as.numeric(scale(AGO[,2], scale=T, center=F))
AGO$r <- c(log(AGO$Nsc[2:dim(popreg)[1]] / AGO$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(AGO)
plot(AGO$Nsc, AGO$r, xlab="N", ylab="r", pch=19)
AGO.use <- AGO[which(AGO$year >= 2010),]
plot(AGO.use$Nsc, AGO.use$r, xlab="N", ylab="r", pch=19)
AGO.rick <- lm(AGO.use$r ~ AGO.use$Nsc)
abline(AGO.rick, lty=2, col="red")
summary(AGO.rick)

# Nigeria
NGA <- read.csv("NGA.csv", header=T)
head(NGA)
NGA$Nsc <- as.numeric(scale(NGA[,2], scale=T, center=F))
NGA$r <- c(log(NGA$Nsc[2:dim(popreg)[1]] / NGA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(NGA)
plot(NGA$Nsc, NGA$r, xlab="N", ylab="r", pch=19)
NGA.use <- NGA[which(NGA$year >= 2010),]
plot(NGA.use$Nsc, NGA.use$r, xlab="N", ylab="r", pch=19)
NGA.rick <- lm(NGA.use$r ~ NGA.use$Nsc)
abline(NGA.rick, lty=2, col="red")
summary(NGA.rick)

# Burundi
BDI <- read.csv("BDI.csv", header=T)
head(BDI)
BDI$Nsc <- as.numeric(scale(BDI[,2], scale=T, center=F))
BDI$r <- c(log(BDI$Nsc[2:dim(popreg)[1]] / BDI$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(BDI)
plot(BDI$Nsc, BDI$r, xlab="N", ylab="r", pch=19)
BDI.use <- BDI[which(BDI$year >= 2008),]
plot(BDI.use$Nsc, BDI.use$r, xlab="N", ylab="r", pch=19)
BDI.rick <- lm(BDI.use$r ~ BDI.use$Nsc)
abline(BDI.rick, lty=2, col="red")
summary(BDI.rick)

# Burkina Faso
BFA <- read.csv("BFA.csv", header=T)
head(BFA)
BFA$Nsc <- as.numeric(scale(BFA[,2], scale=T, center=F))
BFA$r <- c(log(BFA$Nsc[2:dim(popreg)[1]] / BFA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(BFA)
plot(BFA$Nsc, BFA$r, xlab="N", ylab="r", pch=19)
BFA.use <- BFA[which(BFA$year >= 2003),]
plot(BFA.use$Nsc, BFA.use$r, xlab="N", ylab="r", pch=19)
BFA.rick <- lm(BFA.use$r ~ BFA.use$Nsc)
abline(BFA.rick, lty=2, col="red")
summary(BFA.rick)

# Gambia
GMB <- read.csv("GMB.csv", header=T)
head(GMB)
GMB$Nsc <- as.numeric(scale(GMB[,2], scale=T, center=F))
GMB$r <- c(log(GMB$Nsc[2:dim(popreg)[1]] / GMB$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(GMB)
plot(GMB$Nsc, GMB$r, xlab="N", ylab="r", pch=19)
GMB.use <- GMB[which(GMB$year >= 2008),]
plot(GMB.use$Nsc, GMB.use$r, xlab="N", ylab="r", pch=19)
GMB.rick <- lm(GMB.use$r ~ GMB.use$Nsc)
abline(GMB.rick, lty=2, col="red")
summary(GMB.rick)

# Uganda
UGA <- read.csv("UGA.csv", header=T)
head(UGA)
UGA$Nsc <- as.numeric(scale(UGA[,2], scale=T, center=F))
UGA$r <- c(log(UGA$Nsc[2:dim(popreg)[1]] / UGA$Nsc[1:(dim(popreg)[1]-1)]), NA)
head(UGA)
plot(UGA$Nsc, UGA$r, xlab="N", ylab="r", pch=19)
UGA.use <- UGA[which(UGA$year >= 2016),]
plot(UGA.use$Nsc, UGA.use$r, xlab="N", ylab="r", pch=19)
UGA.rick <- lm(UGA.use$r ~ UGA.use$Nsc)
abline(UGA.rick, lty=2, col="red")
summary(UGA.rick)



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


########################################################################
## linear models for examining contribution of per-capita consumption ##
########################################################################

## & population size to temperature anomaly
# import data
TAconN <- read.table("consump.csv", sep=",", header=T)
head(TAconN)

# plot
par(mfrow=c(1,3))
plot(TAconN$pcEconsum, TAconN$TaMED, pch=19, xlab="per-capita consumption", ylab="temperature anomaly")
plot(TAconN$pop, TAconN$TaMED, pch=19, xlab="population size", ylab="temperature anomaly")
plot(TAconN$pop, TAconN$pcEconsum, pch=19, xlab="population size", ylab="per-capita consumption")
par(mfrow=c(1,1))

# models
mod1 <- "TaMED~pop+pcEconsum"
mod2 <- "TaMED~pop"
mod3 <- "TaMED~pcEconsum"
mod4 <- "TaMED~1"

## model vector
mod.vec <- c(mod1,mod2,mod3,mod4)
length(mod.vec)
length(unique(mod.vec))

## define n.mod
n.mod <- length(mod.vec)

# model fitting and logLik output loop
Modnum <- length(mod.vec)
LL.vec <- SaveCount <- AICc.vec <- BIC.vec <- k.vec <- terml <- Rm <- Rc <- rep(0,Modnum)
mod.list <- summ.fit <- coeffs <- coeffs.se <- term.labs <- coeffs.st <- list()
mod.num <- seq(1,Modnum,1)

for(i in 1:Modnum) {
  fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=TAconN, na.action=na.omit)
  assign(paste("fit",i,sep=""), fit)
  mod.list[[i]] <- fit
  print(i)
}

sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
row.names(sumtable) <- mod.vec
summary.table <- sumtable[order(sumtable[,7],decreasing=F),1:9]
summary.table

## saturated residual diagnostic
i <- 1
fit.sat <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=TAconN, na.action=na.omit)

check_model(fit.sat, detrend=F)
diagsat <- check_model(fit.sat, detrend=F)
diagsat$VIF$x
diagsat$VIF$y
plot_model(fit.sat, show.values=T, vline.color = "purple")


## resampling within temperature anomaly confidence interval
iter <- 10000
itdiv <- iter/10

mod1 <- "Ta.it~pop+con"
mod2 <- "Ta.it~pop"
mod3 <- "Ta.it~con"
mod4 <- "Ta.it~1"

## model vector
mod.vec <- c(mod1,mod2,mod3,mod4)

topmod.vec <- wAICc.mod1.vec <- wAICc.mod2.vec <- wAICc.mod3.vec <- wAICc.mod4.vec <-
  DE.mod1.vec <- DE.mod2.vec <- DE.mod3.vec <- DE.mod4.vec <- rep(NA,iter)

for (j in 1:iter) {
  Ta.it <- runif(dim(TAconN)[1],min=TAconN$TaLO, max=TAconN$TaUP)
  TAconN.it1 <- data.frame(TAconN$pcEconsum, TAconN$pop, Ta.it)
  colnames(TAconN.it1) <- c("con", "pop", "Ta.it")
  
  # resample iterated dataset with replacement
  TAconN.it <- TAconN.it1[sample(1:30, replace=T), ]
  
  mod.list <- list()
  for(i in 1:length(mod.vec)) {
    fit <- glm(as.formula(mod.vec[i]),family=gaussian(link="identity"), data=TAconN.it, na.action=na.omit)
    assign(paste("fit",i,sep=""), fit)
    mod.list[[i]] <- fit
  }
  sumtable <- aicW(mod.list, finite = TRUE, null.model = NULL, order = F)
  row.names(sumtable) <- mod.vec
  summary.table <- sumtable[order(sumtable[,7],decreasing=F),1:9]
  summary.table
  
  topmod.vec[j] <- which(mod.vec == row.names(summary.table)[1])
                  
  wAICc.mod1.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[1]),5]
  wAICc.mod2.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[2]),5]
  wAICc.mod3.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[3]),5]
  wAICc.mod4.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[4]),5]
  
  DE.mod1.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[1]),9]
  DE.mod2.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[2]),9]
  DE.mod3.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[3]),9]
  DE.mod4.vec[j] <- summary.table[which(row.names(summary.table) == mod.vec[4]),9]
  
  if (j %% itdiv==0) print(j) 
  
} # end j

topmod.header <- mod.vec[as.numeric(attr(table(topmod.vec), "names"))]
topmod.table <- as.data.frame(table(topmod.vec)/iter)
topmod.table[,1] <- topmod.header
colnames(topmod.table) <- c("model", "%top")
topmod.table

mod1.AICmed <- median(wAICc.mod1.vec, na.rm=T)
mod2.AICmed <- median(wAICc.mod2.vec, na.rm=T)
mod3.AICmed <- median(wAICc.mod3.vec, na.rm=T)
mod4.AICmed <- median(wAICc.mod4.vec, na.rm=T)

mod1.AIClo <- quantile(wAICc.mod1.vec, probs=0.025, na.rm=T)
mod2.AIClo <- quantile(wAICc.mod2.vec, probs=0.025, na.rm=T)
mod3.AIClo <- quantile(wAICc.mod3.vec, probs=0.025, na.rm=T)
mod4.AIClo <- quantile(wAICc.mod4.vec, probs=0.025, na.rm=T)

mod1.AICup <- quantile(wAICc.mod1.vec, probs=0.975, na.rm=T)
mod2.AICup <- quantile(wAICc.mod2.vec, probs=0.975, na.rm=T)
mod3.AICup <- quantile(wAICc.mod3.vec, probs=0.975, na.rm=T)
mod4.AICup <- quantile(wAICc.mod4.vec, probs=0.975, na.rm=T)

mod1.AICmed.st <- mod1.AICmed/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod2.AICmed.st <- mod2.AICmed/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod3.AICmed.st <- mod3.AICmed/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod4.AICmed.st <- mod4.AICmed/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
AICmed.st.vec <- c(mod1.AICmed.st,mod2.AICmed.st,mod3.AICmed.st,mod4.AICmed.st)

mod1.AIClo.st <- mod1.AIClo/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod2.AIClo.st <- mod2.AIClo/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod3.AIClo.st <- mod3.AIClo/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod4.AIClo.st <- mod4.AIClo/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
AIClo.st.vec <- c(mod1.AIClo.st,mod2.AIClo.st,mod3.AIClo.st,mod4.AIClo.st)

mod1.AICup.st <- mod1.AICup/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod2.AICup.st <- mod2.AICup/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod3.AICup.st <- mod3.AICup/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
mod4.AICup.st <- mod4.AICup/sum(c(mod1.AICmed, mod2.AICmed, mod3.AICmed, mod4.AICmed))
AICup.st.vec <- c(mod1.AICup.st,mod2.AICup.st,mod3.AICup.st,mod4.AICup.st)

topmod.index <- c(which(topmod.table$model == mod.vec[1]),which(topmod.table$model == mod.vec[2]),
                  which(topmod.table$model == mod.vec[3]))
topmod.table$AICmdst <- round(AICmed.st.vec[topmod.index],4)
topmod.table$AIClost <- round(AIClo.st.vec[topmod.index],4)
topmod.table$AICupst <- round(AICup.st.vec[topmod.index],4)

mod1.DEmed <- median(DE.mod1.vec, na.rm=T)
mod2.DEmed <- median(DE.mod2.vec, na.rm=T)
mod3.DEmed <- median(DE.mod3.vec, na.rm=T)
mod4.DEmed <- median(DE.mod4.vec, na.rm=T)
DEmed.st.vec <- c(mod1.DEmed,mod2.DEmed,mod3.DEmed,mod4.DEmed)

topmod.table$DEmdst <- round(DEmed.st.vec[topmod.index], 1)
topmod.sort <- topmod.table[order(topmod.table[,2],decreasing=T),]
topmod.sort

## boosted regression tree (median TA)
brt.fit <- gbm.step(TAconN, gbm.x = attr(TAconN, "names")[c(2:3)], gbm.y = attr(TAconN, "names")[3],
                    family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.0001,
                    bag.fraction=0.75, tree.complexity = 2)
summary(brt.fit)
D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) / brt.fit$cv.statistics$deviance.mean
gbm.plot(brt.fit)
gbm.plot.fits(brt.fit)

brt.CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
brt.CV.cor.se <- 100 * brt.fit$cv.statistics$correlation.se
print(c(brt.CV.cor, brt.CV.cor.se))


# resampled BRT loop
biter <- 1000
eq.sp.points <- 100

# create storage arrays
val.arr <- pred.arr <- array(data = NA, dim = c(eq.sp.points, 2, biter),
                             dimnames=list(paste("x",1:eq.sp.points,sep=""),
                             attr(TAconN, "names")[c(2:3)], paste("b",1:biter,sep="")))

# create storage vectors
D2.vec <- CV.cor.vec <- CV.cor.se.vec <- N.ri <- E.ri <- rep(NA,biter)

for (b in 1:biter) {
  # resample data among years
  resamp.sub <- sort(sample(x = 1:dim(TAconN)[1], size = dim(TAconN)[1], replace=TRUE))
  dat.resamp <- TAconN[resamp.sub,]
  dat.resamp$TA.resamp <- runif(dim(dat.resamp)[1],min=dat.resamp$TaLO, max=dat.resamp$TaUP)
  
  # boosted regression tree
  brt.fit <- gbm.step(dat.resamp, gbm.x = attr(dat.resamp, "names")[c(2:3)], gbm.y = attr(dat.resamp, "names")[8],
                      family="gaussian", max.trees=100000, tolerance = 0.0001, learning.rate = 0.001, bag.fraction=0.75,
                      tree.complexity = 2, silent=T, tolerance.method = "auto")
  summ.fit <- summary(brt.fit)
  
  length(summ.fit[[1]])
  
  if (length(summ.fit[[1]]) == 2) {
    # variable relative importance
    E.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(2:3)][1])]
    N.ri[b] <- summ.fit$rel.inf[which(summ.fit$var == attr(dat.resamp, "names")[c(2:3)][2])]

    D2 <- 100 * (brt.fit$cv.statistics$deviance.mean - brt.fit$self.statistics$mean.resid) /
                brt.fit$cv.statistics$deviance.mean
    D2.vec[b] <- D2
    CV.cor <- 100 * brt.fit$cv.statistics$correlation.mean
    CV.cor.vec[b] <- CV.cor
    CV.cor.se <- 100 *brt.fit$cv.statistics$correlation.se
    CV.cor.se.vec[b] <- CV.cor.se
    
    RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.points, ncol=2)
    ## output average predictions
    for (p in 1:2) {
      RESP.val[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,1]
      RESP.pred[,p] <- plot.gbm(brt.fit, i.var=p, continuous.resolution=eq.sp.points, return.grid=T)[,2]
    }
    RESP.val.dat <- as.data.frame(RESP.val)
    colnames(RESP.val.dat) <- brt.fit$var.names
    RESP.pred.dat <- as.data.frame(RESP.pred)
    colnames(RESP.pred.dat) <- brt.fit$var.names
    
    val.arr[, , b] <- as.matrix(RESP.val.dat)
    pred.arr[, , b] <- as.matrix(RESP.pred.dat)
    
    print(b)
  }
  
  if (length(summ.fit[[1]]) != 2) {
    b <- b+1
    print(b)
  }
  
} # end b

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- 5
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | 
                                  pred.update[,,z] > (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  }
  print(k)
}

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)

val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

par(mfrow=c(1,2)) 
plot(val.med[,1],pred.med[,1],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) TA (higher→)", xlab="(←higher) E (lower→)")
lines(val.med[,1], pred.lo[,1], type="l", lty=2, col="red")
lines(val.med[,1], pred.up[,1], type="l", lty=2, col="red")

plot(val.med[,2],pred.med[,2],type="l",ylim=c(min(pred.lo[,2]),max(pred.up[,2])), lwd=2, ylab="(←lower) TA (higher→)",  xlab="(←lower) N (higher→)" )
lines(val.med[,2], pred.lo[,2], type="l", lty=2, col="red")
lines(val.med[,2], pred.up[,2], type="l", lty=2, col="red")
par(mfrow=c(1,1)) 

# kappa method for output vectors
D2.update <- D2.vec[1:biter]
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]
E.ri.update <- E.ri[1:biter]
N.ri.update <- N.ri[1:biter]

for (k in 1:kappa.n) {
  D2.mean <- mean(D2.update, na.rm=T); D2.sd <- sd(D2.update, na.rm=T)
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  E.mean <- mean(E.ri.update, na.rm=T); E.sd <- sd(E.ri.update, na.rm=T)
  N.mean <- mean(N.ri.update, na.rm=T); N.sd <- sd(N.ri.update, na.rm=T)

  for (u in 1:biter) {
    D2.update[u] <- ifelse((D2.update[u] < (D2.mean-kappa*D2.sd) | D2.update[u] > (D2.mean+kappa*D2.sd)), NA, D2.update[u])
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] > (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] > (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    
    E.ri.update[u] <- ifelse((E.ri.update[u] < (E.mean-kappa*E.sd) | E.ri.update[u] > (E.mean+kappa*E.sd)), NA, E.ri.update[u])
    N.ri.update[u] <- ifelse((N.ri.update[u] < (N.mean-kappa*N.sd) | N.ri.update[u] > (N.mean+kappa*N.sd)), NA, N.ri.update[u])
  }
  
  print(k)
}

D2.med <- median(D2.update, na.rm=TRUE)
D2.lo <- quantile(D2.update, probs=0.025, na.rm=TRUE)
D2.up <- quantile(D2.update, probs=0.975, na.rm=TRUE)
print(c(D2.lo,D2.med,D2.up))

CV.cor.med <- median(CV.cor.update, na.rm=TRUE)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=TRUE)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=TRUE)
print(c(CV.cor.lo,CV.cor.med,CV.cor.up))

E.ri.lo <- quantile(E.ri.update, probs=0.025, na.rm=TRUE)
E.ri.med <- median(E.ri.update, na.rm=TRUE)
E.ri.up <- quantile(E.ri.update, probs=0.975, na.rm=TRUE)

N.ri.lo <- quantile(N.ri.update, probs=0.025, na.rm=TRUE)
N.ri.med <- median(N.ri.update, na.rm=TRUE)
N.ri.up <- quantile(N.ri.update, probs=0.975, na.rm=TRUE)

ri.lo <- c(E.ri.lo,N.ri.lo)
ri.med <- c(E.ri.med,N.ri.med)
ri.up <- c(E.ri.up,N.ri.up)

ri.out <- as.data.frame(cbind(ri.lo,ri.med,ri.up))
colnames(ri.out) <- c("ri.lo","ri.med","ri.up")
rownames(ri.out) <- attr(TAconN, "names")[c(2:3)]
ri.sort <- ri.out[order(ri.out[,2],decreasing=T),1:3]
ri.sort


############################################################################
# age structure # UN Population Division (World Population Prospects 2022) #
# both sexes population by single age (1950-2021)                          #
############################################################################

astruct <- read.csv("agestructure1950_2021.csv", header=T)
head(astruct)

totpop <- 1000*apply(astruct[,-1], MARGIN=1, sum, na.rm=T)
age.vec <- seq(0.5,100.5,1)
age.mn <- rep(NA,dim(astruct)[1])
for (i in 1:dim(astruct)[1]) {
  age.mn[i] <- sum((1000*astruct[i, 2:102]) * age.vec)/totpop[i]
} # end a
plot(astruct$year, age.mn, type="l", xlab="year", ylab="mean age (years)")
plot(totpop, age.mn, type="l", ylab="mean age (years)", xlab="human population size")

pl15 <- (1000*apply(astruct[,2:16], MARGIN=1, sum, na.rm=T))/totpop

age.mn.out <- data.frame(astruct$year, totpop/10^9, age.mn, pl15)
colnames(age.mn.out) <- c("year", "Ntot", "ageMN", "young")

plot(age.mn.out$Ntot, age.mn.out$young, type="l", xlab="human population size", ylab="proportion < 15 years")
age.mn.out$year[which(age.mn.out$young == max(age.mn.out$young))]


###################################################
## assume different errors on population estimates
###################################################
head(popdat)

popPhase1 <- subset(popdat, year >= 1800 & year < 1950)
popPhase2 <- subset(popdat, year >= 1950 & year < 1962)
popPhase3 <- subset(popdat, year >= 1962)

## set uncertainties
pcUncertPhase1 <- 0.05 # 5% uncertainty
pcUncertPhase2 <- 0.02 # 2% uncertainty
pcUncertPhase3 <- 0.01 # 1% uncertainty

# calculate stochastic time series of r based on phase-specific certainties & estimate slope of
# relationship between r and Nt
iter <- 10000
Phase1rN.slope <- Phase2rN.slope <- Phase3rN.slope <- Phase3K <- rep(NA, iter)
for (i in 1:iter) {
  # Phase 1: facilitation
  Phase1pop.rsmp <- runif(length(popPhase1$pop), min = popPhase1$pop - (popPhase1$pop * pcUncertPhase1),
        max = popPhase1$pop + (popPhase1$pop * pcUncertPhase1))
  Phase1r.rsmp <- c(NA, log(Phase1pop.rsmp[2:length(Phase1pop.rsmp)] / Phase1pop.rsmp[1:(length(Phase1pop.rsmp)-1)]))
  Phase1rN.fit <- lm(Phase1r.rsmp ~ Phase1pop.rsmp)  
  #summary(Phase1rN.fit)
  Phase1rN.slope[i] <- as.numeric(coef(Phase1rN.fit)[2])
  
  # Phase 2: transition
  Phase2pop.rsmp <- runif(length(popPhase2$pop), min = popPhase2$pop - (popPhase2$pop * pcUncertPhase2),
                          max = popPhase2$pop + (popPhase2$pop * pcUncertPhase2))
  Phase2r.rsmp <- c(NA, log(Phase2pop.rsmp[2:length(Phase2pop.rsmp)] / Phase2pop.rsmp[1:(length(Phase2pop.rsmp)-1)]))
  Phase2rN.fit <- lm(Phase2r.rsmp ~ Phase2pop.rsmp)  
  #summary(Phase2rN.fit)
  Phase2rN.slope[i] <- as.numeric(coef(Phase2rN.fit)[2])
  
  # Phase 3: negative
  Phase3pop.rsmp <- runif(length(popPhase3$pop), min = popPhase3$pop - (popPhase3$pop * pcUncertPhase3),
                          max = popPhase3$pop + (popPhase3$pop * pcUncertPhase3))
  Phase3r.rsmp <- c(NA, log(Phase3pop.rsmp[2:length(Phase3pop.rsmp)] / Phase3pop.rsmp[1:(length(Phase3pop.rsmp)-1)]))
  Phase3rN.fit <- lm(Phase3r.rsmp ~ Phase3pop.rsmp)
  #summary(Phase3rN.fit)
  Phase3rN.slope[i] <- as.numeric(coef(Phase3rN.fit)[2])
  Phase3K[i] <- as.numeric(-coef(Phase3rN.fit)[1]/coef(Phase3rN.fit)[2]) / 10^9
  
}

# calculate 95% confidence intervals
Phase1rN.slope.lo <- quantile(Phase1rN.slope, probs=0.025, na.rm=T)
Phase1rN.slope.up <- quantile(Phase1rN.slope, probs=0.975, na.rm=T)
Phase2rN.slope.lo <- quantile(Phase2rN.slope, probs=0.025, na.rm=T)
Phase2rN.slope.up <- quantile(Phase2rN.slope, probs=0.975, na.rm=T)
Phase3rN.slope.lo <- quantile(Phase3rN.slope, probs=0.025, na.rm=T)
Phase3rN.slope.up <- quantile(Phase3rN.slope, probs=0.975, na.rm=T)

print(c(Phase1rN.slope.lo, Phase1rN.slope.up))
print(c(Phase2rN.slope.lo, Phase2rN.slope.up))
print(c(Phase3rN.slope.lo, Phase3rN.slope.up))

Phase3rN.K.lo <- quantile(Phase3K, probs=0.025, na.rm=T)
Phase3rN.K.up <- quantile(Phase3K, probs=0.975, na.rm=T)
print(c(Phase3rN.K.lo, Phase3rN.K.up))
