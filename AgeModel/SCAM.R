#Rcode for SCAM.R
require(PBSmodelling)
source("read.admb.R")
A=read.admb("SCAM")
isd=grep("sd_noF_sbt",A$fit$names)
sbf0 = A$fit$est[isd]
sbf0.ci = 1.96* A$fit$std[isd]

par(mfcol=c(2, 2), las=1)

plot(A$age, A$la, xlab="Age", ylab="Length (cm)")
plot(A$age, A$wa, xlab="Age", ylab="Weight (kg)")
plot(A$age, A$fa, xlab="Age", ylab="Relative fecundity")
matplot(A$age, t(exp(A$log_sel)), type="l", xlab="Age", ylab="Fishery selectivity") 

#Biomass & Fishing mortality
matplot(A$yrs,cbind(A$bt, A$sbt, sbf0,sbf0-sbf0.ci,sbf0+sbf0.ci  ),type="l", xlab="Year", ylab="Biomass (t)", col=c(1, 1, 2, 2, 2), lty=c(1, 2, 3, 3, 3))
legend("topright", c("total biomass", "spawning biomass", "spawning biomass with no F"), lty=c(1, 2, 3), col=c(1, 1, 2), bty="n")
plot(A$yrs, A$sbt/A$bo, type="l", xlab="Year", ylab="Depletion or Mortality", ylim=c(0, max(A$sbt/A$bo)))
lines(A$yr, A$ft, lty=2)
legend("topright", c("spawning depletion", "fishing mortality"), lty=c(1, 2), col=1, bty="n")

#Survey time series
plot(A$iyr,A$it, xlab="Year", ylab="Relative abundance")
q=exp(mean(log(A$it)-log(A$pit)))
lines(A$iyr, A$pit*q)
legend("topright", c("Observed", "Predicted"), lty=c(-1, 1), pch=c(1, -1), bty="n")


#Stock recruitment
st=seq(0, max(A$sbt), length=100)
rt=A$kappa*A$ro*st/(A$bo+(A$kappa-1)*st)
plot(A$sbt[1:length(A$rt)],A$rt,xlim=c(0,max(A$sbt)),ylim=c(0,max(A$rt)), 
	xlab="Spawning biomass (t)", ylab="Age-1 recruits")
lines(st, rt, type="l")

#Residuals
plot(A$yr, log(A$obs_ct)-log(A$ct),xlab= "Year", ylab="Catch residuals", type="h")
plot(A$iyr, A$epsilon, xlab="Year", ylab="Survey residual", type="h")
plot(log(A$it), log(A$pit)+log(q), xlab="Observed log(It)", ylab="Predicted log(It)")
abline(lm(log(A$pit)+log(q)~log(A$it)+0))
abline(0, 1, col=2, lty=2) 

plot(A$yr[-1], A$delta, type="h", xlab="Year", ylab="Recruitment residual (delta)") 

#Catch-age residuals Fishery
plotBubbles(t(A$P),A$yr,A$age, hide0=TRUE, ylab="Age", size=0.075, main="Fishery")
plotBubbles(t(A$P-A$Phat),A$yr,A$age, hide0=TRUE, ylab="Age", size=0.075)

#Catch-age residuals Survey
plotBubbles(t(A$Q),A$yr,A$age, hide0=TRUE, ylab="Age", size=0.075, main="Survey")
plotBubbles(t(A$Q-A$Qhat),A$yr,A$age, hide0=TRUE, ylab="Age", size=0.075)

runsim<-function(nsim=10)
{
	itheta = NULL
	theta = scan("SCAM.ctl", nlines=1, skip=7)
	for(i in 1:nsim)
	{
		arg = paste("./SCAM -nox -sim", 2*i+101)
		system(arg)
		
		P <- read.fit("SCAM")
		print(P$est[1:4])
		itheta <- rbind(itheta, P$est[1:4])
		
	}
	lr = log2(t(t(itheta)/theta))
	boxplot(lr)
}


