//
// File: DSO78.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue Oct 05 18:48:19 2004
//

#include "DSO78.h"

using namespace std;

/******************************************************************************/

DSO78::DSO78(const ProteicAlphabet * alpha) : ProteinSubstitutionModel(alpha)
{
  #include "DSO78ExchangeabilityCode.cpp"
	#include "DSO78FrequenciesCode.cpp"
	updateMatrices();
}

/******************************************************************************/

string DSO78::getName() const {
	return "Dayhoff et al. (1978) protein substitution model";
}

/******************************************************************************/

