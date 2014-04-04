all: main

main: main.o Sequence.o ExtractSequence.o CodonFrequency.o
	g++ main.o Sequence.o ExtractSequence.o CodonFrequency.o -o main

main.o: main.cpp
	g++ -c main.cpp

Sequence.o: Sequence.cpp Sequence.h
	g++ -c Sequence.cpp

ExtractSequence.o: ExtractSequence.cpp ExtractSequence.h
	g++ -c ExtractSequence.cpp

CodonFrequency.o: CodonFrequency.cpp CodonFrequency.h
	g++ -c CodonFrequency.cpp -std=c++0x
	# -std=c++0x used to compile struct like initialization for codon vector

MinMax.o: MinMax.cpp MinMax.h
	g++ -c MinMax.cpp

clean:
	rm -f *.o main
