COMPILER    = g++ ${OPTFLAG}
LINKER      = g++ ${LDFLAGS}

all: clean BRIDES

BRIDES: main.cpp
	g++ -O3 main.cpp -o brides.exe -fopenmp -lm 
	
clean:
	rm -rf *.o
	rm -rf brides.exe