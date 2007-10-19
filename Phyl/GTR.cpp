//
// File: GTR.cpp
// Created by: Julien Dutheil
// Created on: Tue Oct 25 10:21 2005
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

#include "GTR.h"

// From SeqLib:
#include <Seq/SequenceContainerTools.h>

// From the STL:
#include <cmath>

// From NumCalc:
#include <NumCalc/MatrixTools.h>

/******************************************************************************/

IncludingInterval GTR::PI_CONSTRAINT(0, 1);

/******************************************************************************/

GTR::GTR(
	const NucleicAlphabet * alpha,
	double a,
	double b,
	double c,
	double d,
	double e,
	double piA,
	double piC,
	double piG,
	double piT):
	NucleotideSubstitutionModel(alpha)
{
	_parameters.addParameter(Parameter("a", a, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("b", b, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("c", c, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("d", d, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("e", e, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("piA", piA, &PI_CONSTRAINT));
	_parameters.addParameter(Parameter("piC", piC, &PI_CONSTRAINT));
	_parameters.addParameter(Parameter("piG", piG, &PI_CONSTRAINT));
	_parameters.addParameter(Parameter("piT", piT, &PI_CONSTRAINT));
	updateMatrices();
}

/******************************************************************************/
	
void GTR::updateMatrices()
{
	_a = _parameters.getParameter("a")->getValue();
	_b = _parameters.getParameter("b")->getValue();
	_c = _parameters.getParameter("c")->getValue();
	_d = _parameters.getParameter("d")->getValue();
	_e = _parameters.getParameter("e")->getValue();
	_piA = _parameters.getParameter("piA")->getValue();
	_piC = _parameters.getParameter("piC")->getValue();
	_piG = _parameters.getParameter("piG")->getValue();
	_piT = _parameters.getParameter("piT")->getValue();
  _p = 2*(_a*_piC*_piT+_b*_piA*_piT+_c*_piG*_piT+_d*_piA*_piC+_e*_piC*_piG+_piA*_piG);

  _freq[0] = _piA;
  _freq[1] = _piC;
  _freq[2] = _piG;
  _freq[3] = _piT;
	
  // Exchangeability matrix:
	_exchangeability(0,0) = (-_b*_piT-_piG-_d*_piC)/(_piA * _p);
	_exchangeability(1,0) = _d/_p;
	_exchangeability(0,1) = _d/_p;
	_exchangeability(2,0) = 1/_p;
	_exchangeability(0,2) = 1/_p;
	_exchangeability(3,0) = _b/_p;
	_exchangeability(0,3) = _b/_p;
	_exchangeability(1,1) = (-_a*_piT-_e*_piG-_d*_piA)/(_piC * _p);
	_exchangeability(1,2) = _e/_p;
	_exchangeability(2,1) = _e/_p;
	_exchangeability(1,3) = _a/_p;
	_exchangeability(3,1) = _a/_p;
	_exchangeability(2,2) = (-_c*_piT-_e*_piC-_piA)/(_piG * _p);
	_exchangeability(2,3) = _c/_p;
	_exchangeability(3,2) = _c/_p;
	_exchangeability(3,3) = (-_c*_piG-_a*_piC-_b*_piA)/(_piT * _p);

  AbstractReversibleSubstitutionModel::updateMatrices();
}

/******************************************************************************/

string GTR::getName() const { return string("General Time-Reversible"); }

/******************************************************************************/

void GTR::setFreqFromData(const SequenceContainer & data)
{
	map<int, double> freqs = SequenceContainerTools::getFrequencies(data);
	double t = 0;
	for(unsigned int i = 0; i < _size; i++) t += freqs[i];
	setParameterValue("piA", freqs[0] / t);
	setParameterValue("piC", freqs[1] / t);
	setParameterValue("piG", freqs[2] / t);
	setParameterValue("piT", freqs[3] / t);
  updateMatrices();
}

/******************************************************************************/

