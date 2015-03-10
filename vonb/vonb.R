# vonb.R


linf <- 100
k    <- 0.2
to   <- -0.5
cv   <- 0.08


m  <- 1.5*k
a  <- 0:50

lx <- exp(-m*a)
sa <- plogis(a,5, 3.0)
nobs = 100
aged <- rmultinom(1,nobs,lx*sa)
data <- NULL
for(j in a)
{
	n   <- aged[j+1,]
	
	if(n > 0)
	{
		age <- rep(j,length=n)
		
		len <- linf*(1-exp(-k*(age-to)))
		len <- len + rnorm(n,0,cv*len)
		
		data<- rbind(data,cbind(age,round(len,1)))
	}
}
write(nobs,file="SimVonb.dat")
write.table(data,file="SimVonb.dat",append=TRUE,row.names=FALSE,col.names=FALSE)
write(999,file="SimVonb.dat",append=TRUE)