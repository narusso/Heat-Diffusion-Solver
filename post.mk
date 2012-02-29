main: main.o heat_3D.o nrutil.o gauss_elim.o utilities.o

CMD=./main -q -X1 -Y1 -Z1 -n3 -s1 -p0 -t.01 -r0 -b0

timing: timing_ex timing_im timing_im_new

.PHONY: FTCS CN BE BEj BEgs BEsor timing_ex timing_im timing_im_new timing data

timing_ex: FTCS
FTCS: main
	 @for D in `seq 5 2 59`; do $(CMD) -x$$D -y$$D -z$$D -m$@ -Otimings/$@_$${D}_perf.dat; done

timing_im: CN BE
CN BE: main
	 @for D in `seq 5 2 15`; do $(CMD) -x$$D -y$$D -z$$D -m$@ -Otimings/$@_$${D}_perf.dat; done

timing_im_new: BEj BEgs BEsor
BEj BEgs BEsor: main
	@for D in `seq 50 10 140`; do $(CMD) -x$$D -y$$D -z$$D -m$@ -Otimings/$@_$${D}_perf.dat; done

data: main
	./main -n800 -s100 -t.02 -x40 -y40 -z40 -r2 -b0 -q -odata/BEgs.1 -mBEgs
