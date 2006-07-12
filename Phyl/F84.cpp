//
// File: F84.cpp
// Created by: Julien Dutheil
// Created on: Tue May 23 11:13 2005
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

#include "F84.h"

// From SeqLib:
#include <Seq/SequenceContainerTools.h>

// From the STL:
#include <cmath>

// From NumCalc:
#include <NumCalc/MatrixTools.h>

/******************************************************************************/

F84::F84(
	const NucleicAlphabet * alpha,
	double kappa,
	double piA,
	double piC,
	double piG,
	double piT):
	AbstractSubstitutionModel(alpha),
	NucleotideSubstitutionModel(alpha)
{
	piConstraint = new IncludingInterval(0, 1);
	_parameters.addParameter(Parameter("kappa", kappa, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("piA", piA, piConstraint));
	_parameters.addParameter(Parameter("piC", piC, piConstraint));
	_parameters.addParameter(Parameter("piG", piG, piConstraint));
	_parameters.addParameter(Parameter("piT", piT, piConstraint));

	// Frequences:
	_freq[0] = piA;
	_freq[1] = piC;
	_freq[2] = piG;
	_freq[3] = piT;

	updateMatrices();
}

/******************************************************************************/

F84::~F84() { delete piConstraint; }
	
/******************************************************************************/

void F84::updateMatrices()
{
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	
	_generator(0, 0) = -(                  piC + (1 + kappa / piR) * piG +       piT);
	_generator(1, 1) = -(      piA +                   piG + (1 + kappa / piY) * piT); 
	_generator(2, 2) = -((1 + kappa / piR) * piA +       piC             +       piT);
	_generator(3, 3) = -(      piA + (1 + kappa / piY) * piC +       piG            );

	_generator(1, 0) = piA;
	_generator(3, 0) = piA;
	_generator(0, 1) = piC;
	_generator(2, 1) = piC;
	_generator(1, 2) = piG;
	_generator(3, 2) = piG;
	_generator(0, 3) = piT;
	_generator(2, 3) = piT;
	
	_generator(2, 0) = (1 + kappa / piR) * piA;
	_generator(3, 1) = (1 + kappa / piY) * piC;
	_generator(0, 2) = (1 + kappa / piR) * piG;
	_generator(1, 3) = (1 + kappa / piY) * piT;
	
	// Normalization:
	double r = 1. / (1 - piA * piA - piC * piC - piG * piG - piT*piT + 2. * kappa * (piC * piT / piY + piA * piG / piR));
	MatrixTools::scale(_generator, r);
	
	// Exchangeability:
	_exchangeability(0,0) = _generator(0,0) / piA;
	_exchangeability(0,1) = _generator(0,1) / piC; 
	_exchangeability(0,2) = _generator(0,2) / piG; 
	_exchangeability(0,3) = _generator(0,3) / piT;

	_exchangeability(1,0) = _generator(1,0) / piA; 
	_exchangeability(1,1) = _generator(1,1) / piC; 
	_exchangeability(1,2) = _generator(1,2) / piG; 
	_exchangeability(1,3) = _generator(1,3) / piT; 
	
	_exchangeability(2,0) = _generator(2,0) / piA; 
	_exchangeability(2,1) = _generator(2,1) / piC; 
	_exchangeability(2,2) = _generator(2,2) / piG; 
	_exchangeability(2,3) = _generator(2,3) / piT; 
	
	_exchangeability(3,0) = _generator(3,0) / piA;
	_exchangeability(3,1) = _generator(3,1) / piC; 
	_exchangeability(3,2) = _generator(3,2) / piG; 
	_exchangeability(3,3) = _generator(3,3) / piT;

	// Eigen values:
	_eigenValues[0] = 0;
	_eigenValues[1] = -r * (1+kappa);
	_eigenValues[2] = -r * (1+kappa); 
	_eigenValues[3] = -r;
	
	// Eigen vectors:
	_leftEigenVectors(0,0) = piA;
	_leftEigenVectors(0,1) = piC;
	_leftEigenVectors(0,2) = piG;
	_leftEigenVectors(0,3) = piT;

	_leftEigenVectors(1,0) = 0.;
	_leftEigenVectors(1,1) = piT / piY;
	_leftEigenVectors(1,2) = 0.;
	_leftEigenVectors(1,3) = -piT / piY;

	_leftEigenVectors(2,0) = piG / piR;
	_leftEigenVectors(2,1) = 0.;
	_leftEigenVectors(2,2) = -piG / piR;
	_leftEigenVectors(2,3) = 0.;

	_leftEigenVectors(3,0) = piA*piY / piR;
	_leftEigenVectors(3,1) = -piC;
	_leftEigenVectors(3,2) = piG*piY / piR;
	_leftEigenVectors(3,3) = -piT;

	_rightEigenVectors(0,0) = 1.;
	_rightEigenVectors(0,1) = 0.;
	_rightEigenVectors(0,2) = 1.;
	_rightEigenVectors(0,3) = 1.;
	
	_rightEigenVectors(1,0) = 1.;
	_rightEigenVectors(1,1) = 1.;
	_rightEigenVectors(1,2) = 0.;;
	_rightEigenVectors(1,3) = -piR / piY;

	_rightEigenVectors(2,0) = 1.;
	_rightEigenVectors(2,1) = 0.;
	_rightEigenVectors(2,2) = -piA / piG;
	_rightEigenVectors(2,3) = 1.;

	_rightEigenVectors(3,0) = 1.;
	_rightEigenVectors(3,1) = -piC / piT;
	_rightEigenVectors(3,2) = 0.;
	_rightEigenVectors(3,3) = -piR / piY;
}
	
/******************************************************************************/

double F84::Pij_t(int i, int j, double d) const
{
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (1 - piA * piA - piC * piC - piG * piG - piT*piT + 2. * kappa * (piC * piT / piY + piA * piG / piR));
	double l = r * d;
	double k1 = 1;
	double k2 = kappa + 1;
	double exp1 = exp(-k1*l);
	double exp2 = exp(-k2*l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return piA * (1. + (piY/piR) * exp1) + (piG/piR) * exp2; //A
				case 1 : return piC * (1. -             exp1);                    //C
				case 2 : return piG * (1. + (piY/piR) * exp1) - (piG/piR) * exp2; //G
				case 3 : return piT * (1. -             exp1);                    //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return piA * (1. -             exp1);                    //A
				case 1 : return piC * (1. + (piR/piY) * exp1) + (piT/piY) * exp2; //C
				case 2 : return piG * (1. -             exp1);                    //G
				case 3 : return piT * (1. + (piR/piY) * exp1) - (piT/piY) * exp2; //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return piA * (1. + (piY/piR) * exp1) - (piA/piR) * exp2; //A
				case 1 : return piC * (1. -             exp1);                    //C
				case 2 : return piG * (1. + (piY/piR) * exp1) + (piA/piR) * exp2; //G
				case 3 : return piT * (1. -             exp1);                    //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return piA * (1. -             exp1);                    //A
				case 1 : return piC * (1. + (piR/piY) * exp1) - (piC/piY) * exp2; //C
				case 2 : return piG * (1. -             exp1);                    //G
				case 3 : return piT * (1. + (piR/piY) * exp1) + (piC/piY) * exp2; //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double F84::dPij_dt(int i, int j, double d) const
{
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (1 - piA * piA - piC * piC - piG * piG - piT*piT + 2. * kappa * (piC * piT / piY + piA * piG / piR));
	double l = r * d;
	double k1 = 1;
	double k2 = kappa + 1;
	double exp1 = exp(-k1*l);
	double exp2 = exp(-k2*l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r * (piA * -(piY/piR) * exp1 - (piG/piR) * k2 * exp2); //A
				case 1 : return r * (piC *              exp1);                         //C
				case 2 : return r * (piG * -(piY/piR) * exp1 + (piG/piR) * k2 * exp2); //G
				case 3 : return r * (piT *              exp1);                         //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return r * (piA *              exp1);                         //A
				case 1 : return r * (piC * -(piR/piY) * exp1 - (piT/piY) * k2 * exp2); //C
				case 2 : return r * (piG *              exp1);                         //G
				case 3 : return r * (piT * -(piR/piY) * exp1 + (piT/piY) * k2 * exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r * (piA * -(piY/piR) * exp1 + (piA/piR) * k2 * exp2); //A
				case 1 : return r * (piC *              exp1);                         //C
				case 2 : return r * (piG * -(piY/piR) * exp1 - (piA/piR) * k2 * exp2); //G
				case 3 : return r * (piT *              exp1);                         //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r * (piA *              exp1);                         //A
				case 1 : return r * (piC * -(piR/piY) * exp1 + (piC/piY) * k2 * exp2); //C
				case 2 : return r * (piG *              exp1);                         //G
				case 3 : return r * (piT * -(piR/piY) * exp1 - (piC/piY) * k2 * exp2); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double F84::d2Pij_dt2(int i, int j, double d) const
{
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (1 - piA * piA - piC * piC - piG * piG - piT*piT + 2. * kappa * (piC * piT / piY + piA * piG / piR));
	double r_2 = r * r;
	double l = r * d;
	double k1 = 1;
	double k2 = kappa + 1;
	double k2_2 = k2 * k2;
	double exp1 = exp(-k1*l);
	double exp2 = exp(-k2*l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r_2 * (piA * (piY/piR) * exp1 + (piG/piR) * k2_2 * exp2); //A
				case 1 : return r_2 * (piC *           - exp1);                           //C
				case 2 : return r_2 * (piG * (piY/piR) * exp1 - (piG/piR) * k2_2 * exp2); //G
				case 3 : return r_2 * (piT *           - exp1);                           //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return r_2 * (piA *           - exp1);                           //A
				case 1 : return r_2 * (piC * (piR/piY) * exp1 + (piT/piY) * k2_2 * exp2); //C
				case 2 : return r_2 * (piG *           - exp1);                           //G
				case 3 : return r_2 * (piT * (piR/piY) * exp1 - (piT/piY) * k2_2 * exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r_2 * (piA * (piY/piR) * exp1 - (piA/piR) * k2_2 * exp2); //A
				case 1 : return r_2 * (piC *           - exp1);                           //C
				case 2 : return r_2 * (piG * (piY/piR) * exp1 + (piA/piR) * k2_2 * exp2); //G
				case 3 : return r_2 * (piT *           - exp1);                           //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r_2 * (piA *           - exp1);                           //A
				case 1 : return r_2 * (piC * (piR/piY) * exp1 - (piC/piY) * k2_2 * exp2); //C
				case 2 : return r_2 * (piG *           - exp1);                           //G
				case 3 : return r_2 * (piT * (piR/piY) * exp1 + (piC/piY) * k2_2 * exp2); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

RowMatrix<double> F84::getPij_t(double d) const
{
	RowMatrix<double> p(_size, _size);
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (1 - piA * piA - piC * piC - piG * piG - piT*piT + 2. * kappa * (piC * piT / piY + piA * piG / piR));
	double l = r * d;
	double k1 = 1;
	double k2 = kappa + 1;
	double exp1 = exp(-k1*l);
	double exp2 = exp(-k2*l);

	//A
	p(0, 0) = piA * (1. + (piY/piR) * exp1) + (piG/piR) * exp2; //A
	p(0, 1) = piC * (1. -             exp1);                    //C
	p(0, 2) = piG * (1. + (piY/piR) * exp1) - (piG/piR) * exp2; //G
	p(0, 3) = piT * (1. -             exp1);                    //T, U

	//C
	p(1, 0) = piA * (1. -             exp1);                    //A
	p(1, 1) = piC * (1. + (piR/piY) * exp1) + (piT/piY) * exp2; //C
	p(1, 2) = piG * (1. -             exp1);                    //G
	p(1, 3) = piT * (1. + (piR/piY) * exp1) - (piT/piY) * exp2; //T, U

	//G
	p(2, 0) = piA * (1. + (piY/piR) * exp1) - (piA/piR) * exp2; //A
	p(2, 1) = piC * (1. -             exp1);                    //C
	p(2, 2) = piG * (1. + (piY/piR) * exp1) + (piA/piR) * exp2; //G
	p(2, 3) = piT * (1. -             exp1);                    //T, U

	//T, U
	p(3, 0) = piA * (1. -             exp1);                    //A
	p(3, 1) = piC * (1. + (piR/piY) * exp1) - (piC/piY) * exp2; //C
	p(3, 2) = piG * (1. -             exp1);                    //G
	p(3, 3) = piT * (1. + (piR/piY) * exp1) + (piC/piY) * exp2; //T, U

	return p;
}

RowMatrix<double> F84::getdPij_dt(double d) const
{
	RowMatrix<double> p(_size, _size);
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (1 - piA * piA - piC * piC - piG * piG - piT*piT + 2. * kappa * (piC * piT / piY + piA * piG / piR));
	double l = r * d;
	double k1 = 1;
	double k2 = kappa + 1;
	double exp1 = exp(-k1*l);
	double exp2 = exp(-k2*l);

	//A
	p(0, 0) = r * (piA * -(piY/piR) * exp1 - (piG/piR) * k2 * exp2); //A
	p(0, 1) = r * (piC *              exp1);                         //C
	p(0, 2) = r * (piG * -(piY/piR) * exp1 + (piG/piR) * k2 * exp2); //G
	p(0, 3) = r * (piT *              exp1);                         //T, U

	//C
	p(1, 0) = r * (piA *             exp1);                         //A
	p(1, 1) = r * (piC * -(piR/piY) * exp1 - (piT/piY) * k2 * exp2); //C
	p(1, 2) = r * (piG *             exp1);                         //G
	p(1, 3) = r * (piT * -(piR/piY) * exp1 + (piT/piY) * k2 * exp2); //T, U

	//G
	p(2, 0) = r * (piA * -(piY/piR) * exp1 + (piA/piR) * k2 * exp2); //A
	p(2, 1) = r * (piC *              exp1);                         //C
	p(2, 2) = r * (piG * -(piY/piR) * exp1 - (piA/piR) * k2 * exp2); //G
	p(2, 3) = r * (piT *              exp1);                         //T, U

	//T, U
	p(3, 0) = r * (piA *              exp1);                         //A
	p(3, 1) = r * (piC * -(piR/piY) * exp1 + (piC/piY) * k2 * exp2); //C
	p(3, 2) = r * (piG *              exp1);                         //G
	p(3, 3) = r * (piT * -(piR/piY) * exp1 - (piC/piY) * k2 * exp2); //T, U

	return p;
}

RowMatrix<double> F84::getd2Pij_dt2(double d) const
{
	RowMatrix<double> p(_size, _size);
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (1 - piA * piA - piC * piC - piG * piG - piT*piT + 2. * kappa * (piC * piT / piY + piA * piG / piR));
	double r_2 = r * r;
	double l = r * d;
	double k1 = 1;
	double k2 = kappa + 1;
	double k2_2 = k2 * k2;
	double exp1 = exp(-k1*l);
	double exp2 = exp(-k2*l);

	//A
	p(0, 0) = r_2 * (piA * (piY/piR) * exp1 + (piG/piR) * k2_2 * exp2); //A
	p(0, 1) = r_2 * (piC *           - exp1);                           //C
	p(0, 2) = r_2 * (piG * (piY/piR) * exp1 - (piG/piR) * k2_2 * exp2); //G
	p(0, 3) = r_2 * (piT *           - exp1);                           //T, U

	//C
	p(1, 0) = r_2 * (piA *           - exp1);                           //A
	p(1, 1) = r_2 * (piC * (piR/piY) * exp1 + (piT/piY) * k2_2 * exp2); //C
	p(1, 2) = r_2 * (piG *           - exp1);                           //G
	p(1, 3) = r_2 * (piT * (piR/piY) * exp1 - (piT/piY) * k2_2 * exp2); //T, U

	//G
	p(2, 0) = r_2 * (piA * (piY/piR) * exp1 - (piA/piR) * k2_2 * exp2); //A
	p(2, 1) = r_2 * (piC *           - exp1);                           //C
	p(2, 2) = r_2 * (piG * (piY/piR) * exp1 + (piA/piR) * k2_2 * exp2); //G
	p(2, 3) = r_2 * (piT *           - exp1);                           //T, U

	//T, U
	p(3, 0) = r_2 * (piA *           - exp1);                           //A
	p(3, 1) = r_2 * (piC * (piR/piY) * exp1 - (piC/piY) * k2_2 * exp2); //C
	p(3, 2) = r_2 * (piG *           - exp1);                           //G
	p(3, 3) = r_2 * (piT * (piR/piY) * exp1 + (piC/piY) * k2_2 * exp2); //T, U

	return p;
}

/******************************************************************************/

string F84::getName() const { return string("Felsenstein (1984)"); }

/******************************************************************************/

void F84::setFreqFromData(const SequenceContainer & data)
{
	AbstractSubstitutionModel::setFreqFromData(data);
	// In this model, frequencies may be parameters:
	setParameterValue("piA", _freq[0]);
	setParameterValue("piC", _freq[1]);
	setParameterValue("piG", _freq[2]);
	setParameterValue("piT", _freq[3]);
}

/******************************************************************************/

