EXEC=DDmod
DATA=DDmod.dat
#RSCRIP FOR AUTMATING GRAPHS
RSCRIPT='source(file.path("DDmod.R"))'

all: $(EXEC) run plots

$(EXEC): $(EXEC).tpl
	admb $(EXEC).tpl

run: $(EXEC).cpp
	./$(EXEC) -ind $(DATA)

mcmc: 
	./$(EXEC) -ind $(DATA) -mcmc -mcsave 100
	./$(EXEC) -ind $(DATA) -mceval

clean:
	rm $(EXEC) $(EXEC).[bcehlmprs]* sims variance admodel.* fmin.log

plots: $(EXEC).rep
	-@echo $(RSCRIPT) | R --vanilla --slave

# target: help - Display callable targets.
help:
	egrep "^# target:" [Mm]akefile