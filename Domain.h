/////////////////////////////////////////////////////////////
//
//
//	Domain.h
//
//	Header for domiain class to be used as member of Sequence.
//	
//	Domain represents an area of Protein that is responsible for creating a
//	3D structure in protein i.e. alpha helix, beta sheet
//
//	Created by
//	Sean Howard, for use with Kim Ngo and Xuanyi Li
//
//////////////////////////////////////////////////////////////////


#ifndef DOMAIN_H
#define DOMAIN_H

#include <string>


class Domain{

	public:
		//Input takes a ID for Sequence, Domain type, start and stop position
		Domain(std::string, std::string, int, int);

		//Accessors for Domain memebers
		int getDomainLength();
		int getDomainStart();
		int getDomainEnd();
		std::string getID();
		std::string getDomainType();

	private:
		//Name of Sequence Domain is on (e.g. Hsap|27)
		std::string sequenceID;
		
		//Type of Domain. One letter identifier, most will
		//be A or B, representing an Alpha Helix or Beta Sheet
		std::string domainType;

		//Start and End Positions of the Domain in the larger Sequence
		int domainStart;
		int domainEnd; 



};	

#endif





































