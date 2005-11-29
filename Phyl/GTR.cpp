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
	AbstractSubstitutionModel(alpha),
	NucleotideSubstitutionModel(alpha)
{
	piConstraint = new IncludingInterval(0, 1);
	_parameters.addParameter(Parameter("a", a, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("b", b, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("c", c, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("d", d, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("e", e, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("piA", piA, piConstraint));
	_parameters.addParameter(Parameter("piC", piC, piConstraint));
	_parameters.addParameter(Parameter("piG", piG, piConstraint));
	_parameters.addParameter(Parameter("piT", piT, piConstraint));

	// Frequencies:
	_freq[0] = piA;
	_freq[1] = piC;
	_freq[2] = piG;
	_freq[3] = piT;

	double p = piG*(c*piT+e*piC+  piA)
		       + piC*(b*piT+  piG+d*piA)
					 + piA*(a*piT+e*piG+d*piC)
					 + piT*(c*piG+a*piC+b*piA);
		
	// Exchangeability matrix:
	_exchangeability(0,0) = (-b*piT-  piG-d*piC)/(piA * p);
	_exchangeability(1,0) = d/p;
	_exchangeability(0,1) = d/p;
	_exchangeability(2,0) = 1/p;
	_exchangeability(0,2) = 1/p;
	_exchangeability(3,0) = b/p;
	_exchangeability(0,3) = b/p;
	_exchangeability(1,1) = (-a*piT-e*piG-d*piA)/(piC * p);
	_exchangeability(1,2) = e/p;
	_exchangeability(2,1) = e/p;
	_exchangeability(1,3) = a/p;
	_exchangeability(3,1) = a/p;
	_exchangeability(2,2) = (-c*piT-e*piC-  piA)/(piG * p);
	_exchangeability(2,3) = c/p;
	_exchangeability(3,2) = c/p;
	_exchangeability(3,3) = (-c*piG-a*piC-b*piA)/(piT * p);

	updateMatrices();
}

/******************************************************************************/

GTR::~GTR() { delete piConstraint; }
	
/******************************************************************************/

string GTR::getName() const { return string("General Time-Reversible"); }

/******************************************************************************/

void GTR::setFreqFromData(const SequenceContainer & data)
{
	AbstractSubstitutionModel::setFreqFromData(data);
	// In this model, frequencies may be parameters:
	setParameterValue("piA", _freq[0]);
	setParameterValue("piC", _freq[1]);
	setParameterValue("piG", _freq[2]);
	setParameterValue("piT", _freq[3]);
}

/******************************************************************************/

