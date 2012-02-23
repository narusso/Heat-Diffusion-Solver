flatten: flatten.o nrutil.o gauss_elim.o 
FTCS_1D: FTCS_1D.o nrutil.o
FTCS_2D: FTCS_2D.o nrutil.o
FTCS_3D: FTCS_3D.o nrutil.o
heat_2D: heat_2D.o nrutil.o gauss_elim.o 
main: main.o heat_3D.o nrutil.o gauss_elim.o utilities.o
timing: main
	for M in FTCS BE CN; do \
	 ./main -q -X1 -Y1 -Z1 -x1 -y1 -z1 -n50 -s100 -p0 -t.1 -r0 -b0 -m$$M -o$${M}_soln.dat -O$${M}_perf.dat; done
