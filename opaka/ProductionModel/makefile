EXEC=Schaefer
DATA=Opaka.dat

all: $(EXEC) run

$(EXEC): $(EXEC).tpl
	admb $(EXEC).tpl

run: $(EXEC).cpp
	./$(EXEC) -ind $(DATA)

mcmc: 
	./$(EXEC) -ind $(DATA) -mcmc -mcsave 100
	./$(EXEC) -ind $(DATA) -mceval

clean:
	rm $(EXEC) $(EXEC).[bcehlmprs]* sims variance admodel.* fmin.log

