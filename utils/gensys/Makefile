# Makefile for GENSYS v0.1
#
# (C) Alfredo Metere
# 2012-08-11
#
# Just a simple Makefile
# ------------------------

F77=gfortran
F77FLAGS= -O3
TARGET=gensys

all:
	$(F77) $(F77FLAGS) -o $(TARGET) gensys.f90

clean:
	rm -rf gensys *~ *.out *.csv *.xyz
