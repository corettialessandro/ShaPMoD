IDIR=source/headers
CC=gcc
CFLAGS=-I$(IDIR) -lm -Wall

ODIR=source/obj
SDIR=source

_DEPS = Analysis.h aux_func.h Block_MD.h common.h cells.h ConstTensor.h FirstStep.h Forces.h InitConf.h lattice.h output.h parameters.h SHAKE.h ShellRelaxation.h structs.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = Analysis.o aux_func.o Block_MD.o common.o cells.o ConstTensor.o FirstStep.o Forces.o InitConf.o lattice.o main.o output.o SHAKE.o ShellRelaxation.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

EXEC 	= ShaPMoD

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXEC): $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(EXEC) $(ODIR)/*.o *~ core $(INCDIR)/*~ 
