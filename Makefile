CC = mpicc
CFLAGS =

OBJ = main.o dynamics.o mpi_utils.o euler.o io.o

scm: $(OBJ)
	$(CC) $(CFLAGS) -o scm $(OBJ)

main.o: main.c consts.h dynamics.h mpi_utils.h io.h euler.h
dynamics.o: dynamics.c dynamics.h consts.h
mpi_utils.o: mpi_utils.c mpi_utils.h consts.h
euler.o: euler.c euler.h consts.h
io.o: io.c io.h consts.h

clean:
	rm -f *.o scm
