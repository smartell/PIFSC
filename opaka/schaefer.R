require(ggplot2)

source("read.admb.R")
A <- read.admb("schaefer")

df <- data.frame(year=A$year,bt=A$bt,epsilon=A$epsilon,
                 cpue=A$cpue,yt=A$yt,ut=A$ut)

# plot biomass
p <- ggplot(df,aes(year,bt)) + geom_line()
p <- p + ylim(c(0,max(A$bt)))
print(p)

# plot cpue 
p <- ggplot(df,aes(year,cpue)) + geom_point()
p <- p + geom_line(data=df,aes(year,yt))
print(p)

# plot residuals
p <- ggplot(df,aes(year,epsilon)) + geom_point()
print(p)

# plot exploitation rate
p <- ggplot(df,aes(year,ut)) + geom_line()
print(p)
