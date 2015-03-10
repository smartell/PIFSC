
# Melissa, this is a makefile

all: vonb run

vonb: vonb.tpl
	admb vonb.tpl

run: vonb.cpp
	./vonb

mcmc: 
	./vonb -mcmc -mcsave 1
	./vonb -mceval

clean:
	rm vonb vonb.[bcehlmprs]* sims variance admodel.* fmin.log

