CFLAGS=-Wall -MMD -MP

main : main.o mglin.o addint.o relax.o resid.o rstrc.o interp.o nrutil.o slvsml.o utilities.o

clean:
	rm -f *.o main *~ *.d

SOURCEFILES=$(wildcard *.c)
-include $(SOURCEFILES:.c=.d)
