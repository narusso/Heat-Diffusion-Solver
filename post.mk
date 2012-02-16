flatten: flatten.o nrutil.o gauss_elim.o 
FTCS_1D: FTCS_1D.o nrutil.o
FTCS_2D: FTCS_2D.o nrutil.o
FTCS_3D: FTCS_3D.o nrutil.o
heat_2D: heat_2D.o nrutil.o gauss_elim.o 
main: main.o heat_3D.o nrutil.o gauss_elim.o 
