OBJ = time-fre.o sacio.o

time-fre : $(OBJ)
	cc -o time-fre $(OBJ) -lm -lfftw3

$(OBJ) : sacio.h

clean : 
	rm -f time-fre $(OBJ)
