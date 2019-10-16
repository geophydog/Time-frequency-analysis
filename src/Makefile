BIN=~/Tools/bin/time-fre

OBJ = time-fre.o sacio.o

$(BIN) : $(OBJ)
	cc -o $(BIN) $(OBJ) -lm -lfftw3

$(OBJ) : sacio.h

clean : 
	rm -f time-fre $(OBJ)
