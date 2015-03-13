#Makefile for building and running the LSMR model on a Linux or Mac OSX
EXEC=SCAM
TPL=$(EXEC).tpl
DAT=SCAM.dat

all: $(EXEC) $(EXEC).rep

$(EXEC): $(TPL)
	admb $(EXEC)

$(EXEC).rep: $(TPL) $(DAT)
	./$(EXEC) -ind $(DAT)

clean:
	rm -f $(EXEC) *.eva *.log *.rpt *.htp *.cor *.par *.r* *.b* *.p*
	rm -f *.rep *.bar *.psv *.std $(EXEC).cpp admodel.* variance *.o *.ecm *.mc2 *.mcm