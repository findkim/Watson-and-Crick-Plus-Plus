all: main

main: main.o Sequence.o Alignment.o
	g++ main.o Sequence.o Alignment.o -o main

main.o: main.cpp
	g++ -c main.cpp

Sequence.o: Sequence.cpp
	g++ -c Sequence.cpp

Alignment.o: Alignment.cpp
	g++ -c Alignment.cpp

clean:
	rm -f *.o main
