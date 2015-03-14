EXEC=DDmod
DATA=DDmod.dat
AINP=SimDDmod.pin
CPPS=OperatingModel.cpp
ARGS=
#RSCRIP FOR AUTMATING GRAPHS
RSCRIPT='source(file.path("DDmod.R"))'

all: $(EXEC) run plots

comp: $(EXEC)

$(EXEC): $(EXEC).tpl $(CPPS)
	admb $(EXEC).tpl $(CPPS)

run: $(EXEC).cpp
	./$(EXEC) -ind $(DATA) -ainp $(AINP) $(ARGS)

mse: $(EXEC)
	./$(EXEC) -ind $(DATA) -ainp $(AINP) -mse 999

mcmc: 
	./$(EXEC) -ind $(DATA) -ainp $(AINP) -mcmc -mcsave 100
	./$(EXEC) -ind $(DATA) -ainp $(AINP) -mceval

clean:
	rm $(EXEC) $(EXEC).[bcehlmrs]* sims variance admodel.* fmin.log

plots: $(EXEC).rep
	-@echo $(RSCRIPT) | R --vanilla --slave

data: $(msedirs)
	@echo MSE directories created.

simdirs := $(shell echo 'cat(formatC(1:$(NOSIM), digits=3, flag="0"))' | R --vanilla --slave)
msedirs := $(foreach a, $(S), $(foreach b, $(P), $(addprefix mse_$(a)_$(b)/, $(simdirs))))



