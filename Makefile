all: main

main: main.cpp
	g++ -O3 -fopenmp main.cpp -o main