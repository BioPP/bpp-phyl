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

using namespace bpp;

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
	NucleotideSubstitutionModel(alpha, "GTR.")
{
	Parameter aP("GTR.a", a, &Parameter::R_PLUS_STAR);
	addParameter_(aP);
	Parameter bP("GTR.b", b, &Parameter::R_PLUS_STAR);
	addParameter_(bP);
	Parameter cP("GTR.c", c, &Parameter::R_PLUS_STAR);
	addParameter_(cP);
	Parameter dP("GTR.d", d, &Parameter::R_PLUS_STAR);
	addParameter_(dP);
	Parameter eP("GTR.e", e, &Parameter::R_PLUS_STAR);
	addParameter_(eP);
	_theta = piG + piC;
  _theta1 = piA / (1. - _theta);
  _theta2 = piG / _theta;
	Parameter thetaP("GTR.theta" , _theta , new IncludingInterval(0.001, 0.999), true); //Avoid numerical errors close to the bounds.
	addParameter_(thetaP);
	Parameter theta1P("GTR.theta1", _theta1, new IncludingInterval(0.001, 0.999), true);
	addParameter_(theta1P);
	Parameter theta2P("GTR.theta2", _theta2, new IncludingInterval(0.001, 0.999), true);
	addParameter_(theta2P);
	updateMatrices();
}

/******************************************************************************/
	
void GTR::updateMatrices()
{
	_a = getParameterValue("a");
	_b = getParameterValue("b");
	_c = getParameterValue("c");
	_d = getParameterValue("d");
	_e = getParameterValue("e");
	_theta  = getParameterValue("theta");
	_theta1 = getParameterValue("theta1");
	_theta2 = getParameterValue("theta2");
  _piA = _theta1 * (1. - _theta);
  _piC = (1. - _theta2) * _theta;
  _piG = _theta2 * _theta;
  _piT = (1. - _theta1) * (1. - _theta);
  _p = 2*(_a*_piC*_piT+_b*_piA*_piT+_c*_piG*_piT+_d*_piA*_piC+_e*_piC*_piG+_piA*_piG);

  freq_[0] = _piA;
  freq_[1] = _piC;
  freq_[2] = _piG;
  freq_[3] = _piT;
	
  // Exchangeability matrix:
	exchangeability_(0,0) = (-_b*_piT-_piG-_d*_piC)/(_piA * _p);
	exchangeability_(1,0) = _d/_p;
	exchangeability_(0,1) = _d/_p;
	exchangeability_(2,0) = 1/_p;
	exchangeability_(0,2) = 1/_p;
	exchangeability_(3,0) = _b/_p;
	exchangeability_(0,3) = _b/_p;
	exchangeability_(1,1) = (-_a*_piT-_e*_piG-_d*_piA)/(_piC * _p);
	exchangeability_(1,2) = _e/_p;
	exchangeability_(2,1) = _e/_p;
	exchangeability_(1,3) = _a/_p;
	exchangeability_(3,1) = _a/_p;
	exchangeability_(2,2) = (-_c*_piT-_e*_piC-_piA)/(_piG * _p);
	exchangeability_(2,3) = _c/_p;
	exchangeability_(3,2) = _c/_p;
	exchangeability_(3,3) = (-_c*_piG-_a*_piC-_b*_piA)/(_piT * _p);

  AbstractReversibleSubstitutionModel::updateMatrices();
}

/******************************************************************************/

void GTR::setFreqFromData(const SequenceContainer & data, unsigned int pseudoCount)
{
	map<int, double> freqs = SequenceContainerTools::getFrequencies(data);
	double t = 0;
	for(unsigned int i = 0; i < size_; i++) t += freqs[i] + pseudoCount;
	_piA = (freqs[0] + pseudoCount) / t;
	_piC = (freqs[1] + pseudoCount) / t;
	_piG = (freqs[2] + pseudoCount) / t;
	_piT = (freqs[3] + pseudoCount) / t;
  vector<string> thetas(3);
  thetas[0] = getNamespace() + "theta";
  thetas[1] = getNamespace() + "theta1";
  thetas[2] = getNamespace() + "theta2";
  ParameterList pl = getParameters().subList(thetas);
	pl[0]->setValue(_piC + _piG);
	pl[1]->setValue(_piA / (_piA + _piT));
	pl[2]->setValue(_piG / (_piC + _piG));
  setParametersValues(pl);
}

/******************************************************************************/

