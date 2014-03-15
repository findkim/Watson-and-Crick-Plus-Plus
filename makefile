all: main

main: main.o Sequence.o
	g++ main.o Sequence.o -o main

main.o: main.cpp
	g++ -c main.cpp

Sequence.o: Sequence.cpp
	g++ -c Sequence.cpp

clean:
	rm -f *.o main
