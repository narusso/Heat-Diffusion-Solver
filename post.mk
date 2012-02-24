flatten: flatten.o nrutil.o gauss_elim.o 
FTCS_1D: FTCS_1D.o nrutil.o
FTCS_2D: FTCS_2D.o nrutil.o
FTCS_3D: FTCS_3D.o nrutil.o
heat_2D: heat_2D.o nrutil.o gauss_elim.o 
main: main.o heat_3D.o nrutil.o gauss_elim.o utilities.o
timing: main
	for M in FTCS CN BE; do \
	 for D in `seq 2 4 14`; do time ./main -q -X1 -Y1 -Z1 -x$$D -y$$D -z$$D -n4 -s1 -p0 -t.01 -r0 -b0 -m$$M -O$${M}_$${D}_perf.dat; done; done
