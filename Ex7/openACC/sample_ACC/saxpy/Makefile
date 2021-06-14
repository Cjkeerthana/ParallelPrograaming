PGCC  = mpicc -O3 -ta=tesla -Minfo=all -acc=noautopar

BIN = saxpy

all: saxpy-mpi.c
	$(PGCC) -o $(BIN) $<

clean:
	$(RM) $(BIN) 
