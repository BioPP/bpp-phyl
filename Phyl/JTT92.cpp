//
// File: JTT92.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Jan 21 14:09:43 2004
//

#include "JTT92.h"

using namespace std;

/******************************************************************************/

JTT92::JTT92(const ProteicAlphabet * alpha) : ProteinSubstitutionModel(alpha)
{
  #include "JTT92ExchangeabilityCode.cpp"
	#include "JTT92FrequenciesCode.cpp"
	updateMatrices();	
}

/******************************************************************************/

string JTT92::getName() const {
	return "Jones et al. (1992) protein substitution model";
}

/******************************************************************************/

