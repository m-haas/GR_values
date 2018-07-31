#clear all variables
rm(list=ls(all.names=TRUE))

#CODE to estimate GR-values using MC Simulations

# LIBRARIES
#library(nlstools)
#nlsfit() for emr() function
# FUNCTIONS
fmd <- function(mag,mbin){
  mi <- seq(min(round(mag/mbin)*mbin), max(round(mag/mbin)*mbin), mbin)
  nbm <- length(mi)
  cumnbmag <- numeric(nbm)
  nbmag <- numeric(nbm)
  for(i in 1:nbm) cumnbmag[i] <- length(which(mag > mi[i]-mbin/2))
  cumnbmagtmp <- c(cumnbmag,0)
  nbmag <- abs(diff(cumnbmagtmp))
  res <- list(m=mi, cum=cumnbmag, noncum=nbmag)
  return(res)
}
#Maximum Curvature (MAXC) [e.g., Wiemer & Wyss, 2000]
maxc <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  Mc <- FMD$m[which(FMD$noncum == max(FMD$noncum))[1]]
  return(list(Mc=Mc))
}
#Goodness-of-fit test (GFT) [Wiemer & Wyss, 2000]
gft <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  McBound <- maxc(mag,mbin)$Mc
  Mco <- McBound-0.4+(seq(15)-1)/10
  R <- numeric(15)
  for(i in 1:15){
    indmag <- which(mag > Mco[i]-mbin/2)
    b <- log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    a <- log10(length(indmag))+b*Mco[i]
    FMDcum_model <- 10^(a-b*FMD$m)
    indmi <- which(FMD$m >= Mco[i])
    R[i] <- sum(abs(FMD$cum[indmi]-FMDcum_model[indmi]))/sum(FMD$cum[indmi])*100
    #in Wiemer&Wyss [2000]: 100-R
  }
  indGFT <- which(R <= 5)
  #95% confidence
  if(length(indGFT) != 0){
    Mc <- Mco[indGFT[1]]
    best <- "95%"
  } else{
    indGFT <- which(R <= 10)
    #90% confidence
    if(length(indGFT) != 0){
      Mc <- Mco[indGFT[1]]
      best <- "90%"
    } else{
      Mc <- McBound
      best <- "MAXC"
    }
  }
  return(list(Mc=Mc, best=best, Mco=Mco, R=R))
}
#Mc by b-val Stability (MBS) [Cao & Gao, 2002]
#Modification with Shi & Bolt [1982] uncertainty [Woesner & Wiemer, 2005]
mbs <- function(mag,mbin){
  McBound <- maxc(mag,mbin)$Mc
  Mco <- McBound-0.7+(seq(20)-1)/10
  bi <- numeric(20); unc <- numeric(20)
  for(i in 1:20){
    indmag <- which(mag > Mco[i]-mbin/2)
    nbev <- length(indmag)
    bi[i] <- log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    unc[i] <- 2.3*bi[i]^2*sqrt(sum((mag[indmag]-
                                      mean(mag[indmag]))^2)/(nbev*(nbev-1)))
  }
  bave <- numeric(15)
  for(i in 1:15) bave[i] <- mean(bi[i:i+5])
  dbi_old <- abs(diff(bi))
  indMBS_old <- which(dbi_old <= 0.03)
  dbi <- abs(bave[1:15]-bi[1:15])
  indMBS <- which(dbi <= unc[1:15])
  Mc <- Mco[indMBS[1]]
  return(list(Mc=Mc, Mco=Mco, bi=bi, unc=unc, bave=bave))
}
#Entire Magnitude Range method (EMR) [Woesner & Wiemer, 2005]
emr <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  nbm <- length(FMD$m)
  McMAXC <- maxc(mag,mbin)$Mc
  mu <- abs(McMAXC/2); sig <- abs(McMAXC/4)
  if(mu > 1)mu <- abs(McMAXC/10); sig <- abs(McMAXC/20)
  McBound <- McMAXC
  Mco <- McBound-0.3+(seq(9)-1)/10
  params <- numeric(9*4); dim(params) <- c(9,4)
  #a, b, mu, sigma
  prob <- numeric(9)
  savedmodel <- numeric(9*nbm); dim(savedmodel) <- c(9,nbm)
  for(i in 1:9){
    indmag <- which(mag > Mco[i]-mbin/2)
    nbev <- length(indmag)
    b <- log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    a <- log10(length(indmag))+b*Mco[i]
    cumN <- 10^(a-b*FMD$m)
    params[i,1] <- a; params[i,2] <- b
    cumNtmp <- 10^(a-b*(max(FMD$m)+mbin))
    cumNtmp <- c(cumN, cumNtmp)
    N <- abs(diff(cumNtmp))
    data <- data.frame(N=N, m=FMD$m, Nd=FMD$noncum)
    indLow <- which(FMD$m < Mco[i]); indHigh <- which(FMD$m >= Mco[i])
    dataTest <- data.frame(N=data$N[indLow], m=data$m[indLow], Nd=data$Nd[indLow])
    dataTmp <- data.frame(N=data$N[indHigh], m=data$m[indHigh], Nd=data$Nd[indHigh])
    checkNo0 <- which(dataTest$Nd != 0)
    dataTest <- data.frame(N=dataTest$N[checkNo0], m=dataTest$m[checkNo0],
                           Nd=dataTest$Nd[checkNo0])
    #Nmax <- max(dataTmp$Nd)
    Nmax <- max(dataTest$Nd)
    #Nmax <- dataTest$Nd[length(dataTest$Nd)]
    Mmintmp <- min(dataTest$m)
    dataTest$Nd <- dataTest$Nd/Nmax
    dataTest$m <- dataTest$m-Mmintmp
    data4fit <- data.frame(N=dataTest$Nd, m=dataTest$m)
    #non-linear least squares fit
    nlsfit <- nls(N~pnorm(m, mean=mean, sd=sd), data=data4fit,
                  start=list(mean=mu, sd=sig), control=list(maxiter=100, warnOnly = TRUE))
    params[i,3] <- coef(nlsfit)["mean"]; params[i,4] <- coef(nlsfit)["sd"]
    dataTest$N <- pnorm(dataTest$m, mean=coef(nlsfit)["mean"],
                        sd=coef(nlsfit)["sd"])*Nmax
    dataTest$m <- dataTest$m+Mmintmp
    dataTest$Nd <- dataTest$Nd*Nmax
    dataPred <- data.frame(N=c(dataTest$N, dataTmp$N), m=c(dataTest$m, dataTmp$m),
                           Nd=c(dataTest$Nd, dataTmp$Nd))
    dataPred$N <- round(dataPred$N)
    savedmodel[i,c(checkNo0,indHigh)] <- dataPred$N
    #CHECK EMR METHOD#
    #pdf(paste(wd,"plot_NonCumModel_",Mco[i],".pdf", sep=""))
    #plot(dataPred$m, dataPred$Nd, pch=18, xlab="Magnitude",
    #ylab="Cumulative Number", log="y")
    #points(dataPred$m, dataPred$N, pch=1)
    #abline(v=Mco[i], lty="dashed")
    #legend("topright", c("Data","EMR model"), cex=0.8, lty=c(0,0), pch=c(18,1))
    #dev.off()
    #write.table(dataPred, file=paste(wd, "file_NonCumModel_",Mco[i],
    #".txt", sep=""))
    #Logarithm to the basis of 10 of Poisson probability density
    probtmp <- numeric(nbm)
    CheckNo0 <- which(dataPred$N != 0)
    Pmodel <- dataPred$N[CheckNo0]; Pdata <- dataPred$Nd[CheckNo0]
    probtmp[CheckNo0] <- 1/log(10)*(-Pmodel+Pdata*log(Pmodel)-lgamma(Pdata+1))
    prob[i] <- -sum(probtmp)
  }
  indbestfit <- which(prob == min(prob, na.rm=TRUE))
  res <- list(Mc=Mco[indbestfit], a=params[indbestfit,1], b=params[indbestfit,2],
              mu=params[indbestfit,3], sigma=params[indbestfit,4],
              model=savedmodel[indbestfit,], Mco=Mco, prob=prob)
  return(res)
}

library(stats4)

setwd("/home/mhaas/PhD/oq-deserve/seismicity_rev")
data <- read.csv('catalogue.csv')
data<-data[which(data$magnitude>=3.5),]

#truncation of variation at trc*sigma
trc <- 3
#nr of simulations
n <- 1000
#magnitude bin
mbin <- 0.1

#use MC simulation on magnitude uncertainty from conversion and 0.1 assignment for Mw
bvals <- c()
avals <- c()
Mcs <- c()
#create breaks for histogram
Mmax.theoretical <- max(data$magnitude)+trc*max(data$sigmaMagnitude)
mbreaks <- seq(floor(min(10*data$magnitude))/10,ceiling(10*Mmax.theoretical)/10,mbin)

#matrix for cumulative observed counts
h <- hist(data$magnitude,mbreaks,right=FALSE,include.lowest=TRUE)
counts.cumul <- h$counts
for (i in seq(length(h$counts)-1,1)){counts.cumul[i] <- counts.cumul[i]+counts.cumul[i+1]}
counts <- h$mids 
counts <- rbind(counts,counts.cumul)

#actual sample loop
for (i in seq(n)){
#copy catalog
cat <- data[which(data$year>=1990),]
#sample distribution of magnitude
#idxs<-which(cat$sigmaMagnitude==0)
idxs<-which(cat$magnitude>0)
#uniform for rounding (all)
for(i in idxs){cat$magnitude[i]<-runif(n=1,min=cat$magnitude[i]-0.05,max=cat$magnitude[i]+0.05)}
#gaussian for converted
idxs <- which(cat$sigmaMagnitude>0)
random_nrs <- runif(length(cat$magnitude[idxs]),min=pnorm(-1*trc,mean=0,sd=1),max=pnorm(trc,mean=0,sd=1))
cat$magnitude[idxs] <- qnorm(random_nrs,mean=cat$magnitude[idxs],sd=cat$sigmaMagnitude[idxs])
###################################
#Calculation of Mc (Amorese 2007)
###################################
#Magnitude bin
nbsample <- 200 #Bootstrapping
## READ CATALOG ##
#For a catalog with data listed in columns separated by space or tab
#(Longitude, Latitude, Magnitude, etc...)
#Other formats may require a different R function
#cat <- read.csv(cat_file, header=TRUE)
## COMPUTE Mc ##
Mc_bootstrap <- numeric(nbsample)
#select function: maxc(), gft(), mbs(), emr()
#For mbass(), see algorithm Amorese [2007]
for(i in 1:nbsample) Mc_bootstrap[i] <- maxc(sample(cat$magnitude, replace=TRUE),mbin)$Mc
#when using emr(), the loop may break due to failure of nlsfit(),
#in this case use:
#for(i in 1:nbsample) Mc_bootstrap[i] <- as.numeric(try(emr(sample(mag, replace=TRUE),mbin)$Mc))
Mc <- mean(Mc_bootstrap, na.rm=TRUE)

#Get unchanged catalogue
cat <- data[which(data$year>=1990),]
#sample distribution of magnitude
idxs1 <- which(cat$sigmaMagnitude==0)
idxs2 <- which(cat$sigmaMagnitude>0)
#consider rounding error of 0.1 for Mw
cat$sigmaMagnitude[idxs1]=0.1/2

cat$magnitude <- qnorm(random_nrs,mean=cat$magnitude,sd=cat$sigmaMagnitude)

#reduce to >=MC rounded up to .1 digit
#Mc<-ceiling(Mc*10)/10
idxs <- which(cat$magnitude>=Mc)
cat <- cat[idxs,]

#get b-value
m <- length(idxs)
b<-log10(exp(1))/(sum(cat$magnitude)/m-Mc)

#get a-value
#log10(m)=a-bM
#log10(m)+bM=a
period <- max(cat$year)+1-min(cat$year)
a <- log10(m/period)+b*Mc

#get counts
h <- hist(cat$magnitude,mbreaks,right=FALSE,include.lowest=TRUE)
counts.cumul <- h$counts
for (i in seq(length(h$counts)-1,1)){counts.cumul[i] <- counts.cumul[i]+counts.cumul[i+1]}

#append avals, bvals, Mcs, counts
avals <- c(avals,a)
bvals <- c(bvals,b)
Mcs <- c(Mcs,Mc)
counts <- rbind(counts,counts.cumul)
}#end of MC simulation

#get
q025<-c()
q250<-c()
q500<-c()
q750<-c()
q975<-c()
#for each mbin get 2.5,25,50,75,97.5 quantiles
l<-length(counts[,1])
w <- which(counts[1,]>=median(Mcs))
#exclude 1st row=mbins and 2nd row=mean_cat
#and exclude below Mc
for (i in w){
  q025 <- c(q025,quantile(counts[3:l,i]/period,0.025))
  q250 <- c(q250,quantile(counts[3:l,i]/period,0.25))
  q500 <- c(q500,quantile(counts[3:l,i]/period,0.5))
  q750 <- c(q750,quantile(counts[3:l,i]/period,0.75))
  q975 <- c(q975,quantile(counts[3:l,i]/period,0.975))
}
#dataframe
df <- data.frame(Mw=counts[1,w],q500,q025,q250,q750,q975)
write.csv(df,'GR_obs.csv',row.names=F,append=F)

#estimated
ests <- 10**(median(avals)-median(bvals)*counts[1,w])
ests_a025b975 <- 10**(quantile(avals,0.025)-quantile(bvals,0.975)*counts[1,w])
ests_a975b025 <- 10**(quantile(avals,0.975)-quantile(bvals,0.025)*counts[1,w])
df <- data.frame(Mw=counts[1,w],ests,ests_a025b975,ests_a975b025)
write.csv(df,'GR_est.csv',row.names=F,append=F)

#normal test
#
#test <- dnorm(bvals,mean=mean(bvals),sd=sd(bvals))

#fit distribution
library(MASS)

test <- ecdf(sort(bvals))
test_pdf <- test_cdf
for (i in seq(length(test_cdf)-1,1,-1)){test_pdf[i]<-test_pdf[i+1]-test_pdf[i]}

test_x <- c(1.10,1.16,1.22,1.24,1.30)
h<-hist(bvals,c(1.009,1.14,1.20,1.23,1.26,1.332))
d<-density(bvals)
plot(h)
lines(d)
points(h$mids,c(0.8,4,8.1,12.1,1.1))
test()

df <- data.frame(d$x,d$y)
write.csv(df,'b_dist.csv',row.names=F)
df <- data.frame(h$mids,h$counts,h$breaks)
write.csv(df,'b_hist.csv',row.names=F)

# hist(bvals,c(1.009,1.19,1.255,1.332))
# lines(density(bvals))





# #plot GR-law
# #get counts for bins
# #nr of bins
# mbin <- 0.1
# data <- cat
# data<-data[which(data$magnitude>=mean(Mcs)),]
# #data<-data[which(data$year>=1990),] 
# mbreaks <- seq(floor(min(10*data$magnitude))/10,ceiling(10*max(data$magnitude))/10,mbin)
# #histogram
# h <- hist(data$magnitude,mbreaks,right=FALSE,include.lowest=TRUE)
# 
# counts.cumul <- h$counts
# for (i in seq(length(h$counts)-1,1)){counts.cumul[i] <- counts.cumul[i]+counts.cumul[i+1]}
# 
# #get percentiles for whisker
# 
# 
# #get estimates from GR-law
# ests<-10**(mean(avals)-mean(bvals)*h$mids)
# 
# #theoretical b=1
# b1 <- 10**(mean(avals)-h$mids)
# 
# #plot
# plot(h$mids,counts.cumul/period,log="y")
# lines(h$mids,ests,col='red')
# lines(h$mids,b1,col='blue')






# yc <- c(1985,1970,1925,1900,1900,1900,1900)
# mc <- c(4.0,4.5,5.0,5.5,6.0,6.5,7.0)
# 
# count<-c()
# for (i in seq(length(mc))){
#   idxs1 <- which(data$year > yc[i])
#   idxs2 <- which(data$magnitude >= mc[i])
#   count <- c(count,length(intersect(idxs1,idxs2))) 
# }
# 
# inc_count <- count
# for (i in seq(length(inc_count)-1,1)){
#   inc_count[i]=inc_count[i]-inc_count[i+1]
# }
# df <- data.frame(yc,mc=mc+0.25,inc_count)
# #write.csv(df,'weichert_in',row.names=F)
# 
# #result from fortran weichert
# b <- 1.1
# a <- 6.15
# #a <- log10(1307)+b*4.25
# est <- 10**(a-b*(mc+0.25))
# 
# obs <- count/(2015-yc)
# 
# plot(mc,obs,log="y")
# lines(mc,est,col='red')

# LL <- function(x,y,beta0, beta1, mu, sigma) {
#   R = y - x * beta1 - beta0
#   #
#   R = suppressWarnings(dnorm(R, mu, sigma, log = TRUE))
#   #
#   -sum(R)
# }
# 
# fit <- lm(log10(y)~x,df)
# fit <- mle(LL(mc,count), start = list(beta0 = 5, beta1 = 1, mu = 0, sigma=1))

