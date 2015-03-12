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

mdf <- melt(df,id.var="year")
ssdf <- mdf %>% subset(variable %in% c("epsilon","nu","delta","fdev"))
p <- ggplot(ssdf,aes(year,value)) + geom_point() 
p <- p + facet_wrap(~variable)
print(p + theme_bw())
