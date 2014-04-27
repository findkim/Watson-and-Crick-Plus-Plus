//////////////////////////////////////////////////
//
//
//	Domain.cpp
//
//	Implementaion of Domain class, to be used as an object member of 
//	Sequence class.
//
//	Represents a domain on a protein, area responsible for folding of macro
//	structures
//
//	Created by,
//	Sean Howard, for use with Kim Ngo and Xuanyi Li
//
/////////////////////////////////////////////////////


#include "Domain.h"
#include <string>

//Constructor feeds in members
Domain::Domain(std::string ID, std::string dType, int start, int end): sequenceID(ID), domainType(dType), domainStart(start), domainEnd(end)
{

}	

//Basic Accessor Functions
int Domain::getDomainLength()
{
	return domainEnd-domainStart;
}

int Domain::getDomainStart()
{
	return domainStart;
}

int Domain::getDomainEnd()
{
	return domainEnd;
}

std::string Domain::getID()
{
	return sequenceID;
}

std::string Domain::getDomainType()
{
	return domainType;
}






















