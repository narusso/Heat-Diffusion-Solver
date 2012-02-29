main: main.o heat_3D.o nrutil.o gauss_elim.o utilities.o

CMD=./main -q -X1 -Y1 -Z1 -n3 -s1 -p0 -t.01 -r0 -b0

timing: timing_explicit timing_implicit timing_implicit_iterative

.PHONY: FTCS CN BE BEj BEgs BEsor timing_explicit timing_implicit timing_implicit_iterative timing data

timing_explicit: FTCS
timing_implicit: CN BE
timing_implicit_iterative: BEj BEgs BEsor

FTCS: main
	 @for D in `seq 5 2 59`; do $(CMD) -x$$D -y$$D -z$$D -m$@ -Otimings/$@_$${D}_perf.dat; done

CN BE: main
	 @for D in `seq 5 2 15`; do $(CMD) -x$$D -y$$D -z$$D -m$@ -Otimings/$@_$${D}_perf.dat; done

BEj BEgs: main
	@for D in `seq 320 320`; do $(CMD) -x$$D -y$$D -z$$D -m$@ -Otimings/$@_$${D}_perf.dat; done

BEsor: main
	@for w in `seq 1.12 .02 1.2`; do for D in `seq 280 20 320`; do $(CMD) -w$${w} -x$$D -y$$D -z$$D -m$@ -Otimings/$@-$${w}_$${D}_perf.dat; done; done

data: main
	./main -n800 -s100 -t.02 -x40 -y40 -z40 -r2 -b0 -q -odata/BEgs.1 -mBEgs
