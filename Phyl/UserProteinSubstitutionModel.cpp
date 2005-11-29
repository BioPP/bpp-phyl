//
// File: UserProteinSubstitutionModel.cpp
// Created by: Julien Dutheil
// Created on: Wed Aug 26 16:27 2005
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "UserProteinSubstitutionModel.h"

// From Utils:
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/StringTokenizer.h>

// From NumCalc:
#include <NumCalc/VectorTools.h>
using namespace VectorFunctions;

// From the STL:
#include <iostream>
#include <fstream>
using namespace std;

/******************************************************************************/

UserProteinSubstitutionModel::UserProteinSubstitutionModel(const ProteicAlphabet * alpha, const string & path) : 
	AbstractSubstitutionModel(alpha),
	ProteinSubstitutionModel(alpha),
	_path(path)
{
	
	readFromFile();
	updateMatrices();	
}

/******************************************************************************/

string UserProteinSubstitutionModel::getName() const
{
	return "User model from file '" + _path + "'";
}

/******************************************************************************/

void UserProteinSubstitutionModel::readFromFile()
{
 	ifstream in(_path.c_str(), ios::in);
	//Read exchangeability matrix:
	for(unsigned int i = 1; i < 20; i++) {
		string line = FileTools::getNextLine(in);
		StringTokenizer st(line);
		for(unsigned int j = 0; j < i; j++) {
			double s = TextTools::toDouble(st.nextToken());
			_exchangeability(i,j) = _exchangeability(j,i) = s;
		}
	}
	//Read frequencies:
	unsigned int fCount = 0;
	while(in && fCount < 20) {
		string line = FileTools::getNextLine(in);
		StringTokenizer st(line);
		while(st.hasMoreToken() && fCount < 20) {
			_freq[fCount] = TextTools::toDouble(st.nextToken());
			fCount++;
		}
	}
	double sf = sum(_freq);
	if(sf - 1 > 0.000001) {
		ApplicationTools::displayMessage("WARNING!!! Frequencies sum to " + TextTools::toString(sf) + ", frequencies have been scaled.");
		sf *= 1./sf;
	}

	//Now build diagonal of the exchangeability matrix:
	for(unsigned int i = 0; i < 20; i++) {
		double sum = 0;
		for(unsigned int j = 0; j < 20; j++) {
			if(j!=i) sum += _exchangeability(i,j);
		}
		_exchangeability(i,i) = -sum;
	}

	//Closing stream:
	in.close();
}

/******************************************************************************/

