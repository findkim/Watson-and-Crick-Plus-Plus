all: main

main: main.o Sequence.o ExtractSequence.o AlignmentCIGAR.o NeighborJoining.o PairAlignment.o
	g++ main.o Sequence.o ExtractSequence.o AlignmentCIGAR.o NeighborJoining.o PairAlignment.o -o main

main.o: main.cpp
	g++ -c main.cpp

Sequence.o: Sequence.cpp
	g++ -c Sequence.cpp

ExtractSequence.o: ExtractSequence.cpp
	g++ -c ExtractSequence.cpp

AlignmentCIGAR.o: AlignmentCIGAR.cpp
	g++ -c AlignmentCIGAR.cpp

NeighborJoining.o: NeighborJoining.cpp
	g++ -c NeighborJoining.cpp

PairAlignment.o: PairAlignment.cpp
	g++ -c PairAlignment.cpp

clean:
	rm -f *.o main
