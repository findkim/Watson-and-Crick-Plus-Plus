all: main

main: main.o Sequence.o ExtractSequence.o AlignmentCIGAR.o CodonFrequency.o MinMax.o
	g++ main.o Sequence.o ExtractSequence.o AlignmentCIGAR.o CodonFrequency.o MinMax.o -o main

main.o: main.cpp
	g++ -c main.cpp

Sequence.o: Sequence.cpp Sequence.h
	g++ -c Sequence.cpp

ExtractSequence.o: ExtractSequence.cpp ExtractSequence.h
	g++ -c ExtractSequence.cpp

AlignmentCIGAR.o: AlignmentCIGAR.cpp
	g++ -c AlignmentCIGAR.cpp

CodonFrequency.o: CodonFrequency.cpp CodonFrequency.h
	g++ -g -c CodonFrequency.cpp #-std=c++0x
	# -std=c++0x used to compile struct like initialization for codon vector

MinMax.o: MinMax.cpp MinMax.h
	g++ -c MinMax.cpp


clean:
	rm -f *.o main
