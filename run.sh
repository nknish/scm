#!/bin/bash
make clean
make
mpirun -n 4 scm
