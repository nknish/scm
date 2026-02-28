CC = mpicc
CFLAGS =

OBJ = main.o dynamics.o mpi_utils.o rk2.o io.o

scm: $(OBJ)
	$(CC) $(CFLAGS) -o scm $(OBJ)

main.o: main.c consts.h dynamics.h mpi_utils.h io.h rk2.h
dynamics.o: dynamics.c dynamics.h consts.h
mpi_utils.o: mpi_utils.c mpi_utils.h consts.h
rk2.o: rk2.c rk2.h consts.h
io.o: io.c io.h consts.h

clean:
	rm -f *.o scm
