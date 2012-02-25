flatten: flatten.o nrutil.o gauss_elim.o
FTCS_1D: FTCS_1D.o nrutil.o
FTCS_2D: FTCS_2D.o nrutil.o
FTCS_3D: FTCS_3D.o nrutil.o
heat_2D: heat_2D.o nrutil.o gauss_elim.o
main: main.o heat_3D.o nrutil.o gauss_elim.o utilities.o
CMD=./main -q -X1 -Y1 -Z1 -n3 -s1 -p0 -t.01 -r0 -b0
timing: timing_ex timing_im
timing_ex: main
	for M in FTCS; do \
	 for D in `seq 5 1 49`; do $(CMD) -x$$D -y$$D -z$$D -m$$M -Otimings/$${M}_$${D}_perf.dat; done; done
timing_im: main
	for M in CN BE; do \
	 for D in `seq 5 1 15`; do $(CMD) -x$$D -y$$D -z$$D -m$$M -Otimings/$${M}_$${D}_perf.dat; done; done
