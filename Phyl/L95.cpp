//
// File: L95.cpp
// Created by: Julien Dutheil
// Created on: Tue Nov 4 11:46 2008
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

#include "L95.h"

// From SeqLib:
#include <Seq/SequenceContainerTools.h>

// From the STL:
#include <cmath>

// From NumCalc:
#include <NumCalc/MatrixTools.h>

using namespace bpp;

/******************************************************************************/

L95::L95(
	const NucleicAlphabet * alpha,
	double beta,
	double gamma,
	double delta,
	double theta):
	NucleotideSubstitutionModel(alpha, "L95.")
{
	Parameter betaP("L95.beta" , beta , &Parameter::R_PLUS_STAR);
	addParameter_(betaP);
	Parameter gammaP("L95.gamma", gamma, &Parameter::R_PLUS_STAR);
	addParameter_(gammaP);
	Parameter deltaP("L95.delta", delta, &Parameter::R_PLUS_STAR);
	addParameter_(deltaP);
	Parameter thetaP("L95.theta" , theta , &Parameter::PROP_CONSTRAINT_EX);
	addParameter_(thetaP);
	updateMatrices();
}

/******************************************************************************/
	
void L95::updateMatrices()
{
	_beta  = getParameterValue("beta");
	_gamma = getParameterValue("gamma");
	_delta = getParameterValue("delta");
	_theta = getParameterValue("theta");

  freq_[0] = _piA = (1. - _theta)/2.;
  freq_[1] = _piC = _theta/2.;
  freq_[2] = _piG = _theta/2;
  freq_[3] = _piT = (1. - _theta)/2.;
	
  // Exchangeability matrix:
	exchangeability_(0,0) = -_gamma*_piT-_piG-_beta*_piC;
	exchangeability_(1,0) = _beta;
	exchangeability_(0,1) = _beta;
	exchangeability_(2,0) = 1.;
	exchangeability_(0,2) = 1.;
	exchangeability_(3,0) = _gamma;
	exchangeability_(0,3) = _gamma;
	exchangeability_(1,1) = -_piT-_delta*_piG-_beta*_piA;
	exchangeability_(1,2) = _delta;
	exchangeability_(2,1) = _delta;
	exchangeability_(1,3) = 1.;
	exchangeability_(3,1) = 1.;
	exchangeability_(2,2) = -_beta*_piT-_delta*_piC-_piA;
	exchangeability_(2,3) = _beta;
	exchangeability_(3,2) = _beta;
	exchangeability_(3,3) = -_beta*_piG-_piC-_gamma*_piA;

  AbstractReversibleSubstitutionModel::updateMatrices();
}

/******************************************************************************/

void L95::setFreqFromData(const SequenceContainer & data)
{
	map<int, double> freqs = SequenceContainerTools::getFrequencies(data);
	double t = 0;
	for(unsigned int i = 0; i < size_; i++) t += freqs[i];
	_piC = freqs[1] / t;
	_piG = freqs[2] / t;
	getParameter_("theta").setValue(_piC + _piG);
  updateMatrices();
}

/******************************************************************************/

