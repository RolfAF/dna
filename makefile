all: dna
	./dna;
dna: main.cpp
	gcc -fopenmp main.cpp -o dna -lstdc++;

