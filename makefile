all: main

main: main.o Sequence.o ExtractSequence.o AlignmentCIGAR.o
	g++ main.o Sequence.o ExtractSequence.o AlignmentCIGAR.o -o main

main.o: main.cpp
	g++ -c main.cpp

Sequence.o: Sequence.cpp Sequence.h
	g++ -c Sequence.cpp

ExtractSequence.o: ExtractSequence.cpp ExtractSequence.h
	g++ -c ExtractSequence.cpp

AlignmentCIGAR.o: AlignmentCIGAR.cpp
	g++ -c AlignmentCIGAR.cpp


clean:
	rm -f *.o main
