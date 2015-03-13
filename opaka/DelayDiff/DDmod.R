library(ggplot2)
library(dplyr)
library(reshape)

source("read.admb.R")
A <- read.admb("DDmod")


df <- data.frame(year=A$year,
                 bt = A$bt,
                 rt = A$rt,
                 nt = A$nt,
                 epsilon=A$epsilon,
                 nu = A$nu,
                 delta= A$delta,
                 psi  = A$psi,
                 wt = A$wt,
                 what=A$what,
                 ft = A$ft,
                 fdev = A$fdev,
                 cpue = A$cpue,
                 yt = A$yt)

p <- ggplot(df,aes(year,bt)) + geom_line()
p <- p + labs(x="Year",y="Biomass (Mlb)")
print(p + theme_bw())

p <- ggplot(df,aes(year,cpue))+geom_point()
p <- p + geom_line(data=df,aes(year,yt))
p <- p + labs(x="Year",y="CPUE")
print(p + theme_bw())

p <- ggplot(df,aes(year,wt))+geom_point()
p <- p + geom_line(data=df,aes(year,what))
p <- p + labs(x="Year",y="Average Weight (lb)")
print(p + theme_bw())

p <- ggplot(df,aes(year,ft))+geom_line()
p <- p + labs(x="Year",y="Fishing mortality rate")
print(p + theme_bw())


mdf <- melt(df,id.var="year")
ssdf <- mdf %>% subset(variable %in% c("epsilon","nu","delta","psi"))
p <- ggplot(ssdf,aes(year,value)) + geom_point() 
p <- p + facet_wrap(~variable)
print(p + theme_bw())

runSim<- function(n=10)
{
        system("rm SimPars.rep")
        for(iter in 1:n)
        {
                arg = paste("./DDMod -nox -est -sim",iter)
                system(arg);
        }
}

