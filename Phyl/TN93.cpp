//
// File: TN93.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thu Jan 22 10:26:51 2004
//

#include "TN93.h"

// From SeqLib:
#include <Seq/SequenceContainerTools.h>

// From the STL:
#include <cmath>

// From the MTL:
#include <mtl/mtl.h>

/******************************************************************************/

TN93::TN93(
	const NucleicAlphabet * alpha,
	double kappa1,
	double kappa2,
	double piA,
	double piC,
	double piG,
	double piT):
	NucleotideSubstitutionModel(alpha)
{
	piConstraint = new IncludingInterval(0, 1);
	_parameters.addParameter(Parameter("kappa1", kappa1, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("kappa2", kappa2, &Parameter::R_PLUS));
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

TN93::~TN93() { delete piConstraint; }
	
/******************************************************************************/

void TN93::updateMatrices()
{
	double kappa1 = _parameters.getParameter("kappa1") -> getValue();
	double kappa2 = _parameters.getParameter("kappa2") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	
	_generator(0, 0) = -(                    piC + kappa1*piG +        piT);
	_generator(1, 1) = -(       piA +                     piG + kappa2*piT); 
	_generator(2, 2) = -(kappa1*piA +        piC              +        piT);
	_generator(3, 3) = -(       piA + kappa2*piC +        piG             );

	_generator(0, 1) = piA;
	_generator(0, 3) = piA;
	_generator(1, 0) = piC;
	_generator(1, 2) = piC;
	_generator(2, 1) = piG;
	_generator(2, 3) = piG;
	_generator(3, 0) = piT;
	_generator(3, 2) = piT;
	
	_generator(0, 2) = kappa1 * piA;
	_generator(1, 3) = kappa2 * piC;
	_generator(2, 0) = kappa1 * piG;
	_generator(3, 1) = kappa2 * piT;
	
	// Normalization:
	double r = 1. / (2. * (piA * piC + piC * piG + piA * piT + piG * piT + kappa2 * piC * piT + kappa1 * piA * piG));
	scale(_generator, r);
	
	// Eigen values:
	_eigenValues[0] = 0;
	_eigenValues[1] = -r;
	_eigenValues[2] = -r * (kappa2 * piY + piR);
	_eigenValues[3] = -r * (kappa1 * piR + piY); 
	
	// Eigen vectors:
	//TODO!!!

}
	
/******************************************************************************/

double TN93::Pij_t(int i, int j, double d) const
{
	double kappa1 = _parameters.getParameter("kappa1") -> getValue();
	double kappa2 = _parameters.getParameter("kappa2") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (2. * (piA * piC + piC * piG + piA * piT + piG * piT + kappa2 * piC * piT + kappa1 * piA * piG));
	double l = r * d;
	double k1 = kappa2 * piY + piR;
	double k2 = kappa1 * piR + piY;
	double exp1 = exp(-l);
	double exp22 = exp(-k2 * l);
	double exp21 = exp(-k1 * l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return piA * (1. + (piY/piR) * exp1) + (piG/piR) * exp22; //A
				case 1 : return piC * (1. -             exp1);                     //C
				case 2 : return piG * (1. + (piY/piR) * exp1) - (piG/piR) * exp22; //G
				case 3 : return piT * (1. -             exp1);                     //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return piA * (1. -             exp1);                     //A
				case 1 : return piC * (1. + (piR/piY) * exp1) + (piT/piY) * exp21; //C
				case 2 : return piG * (1. -             exp1);                     //G
				case 3 : return piT * (1. + (piR/piY) * exp1) - (piT/piY) * exp21; //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return piA * (1. + (piY/piR) * exp1) - (piA/piR) * exp22; //A
				case 1 : return piC * (1. -             exp1);                     //C
				case 2 : return piG * (1. + (piY/piR) * exp1) + (piA/piR) * exp22; //G
				case 3 : return piT * (1. -             exp1);                     //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return piA * (1. -             exp1);                     //A
				case 1 : return piC * (1. + (piR/piY) * exp1) - (piC/piY) * exp21; //C
				case 2 : return piG * (1. -             exp1);                     //G
				case 3 : return piT * (1. + (piR/piY) * exp1) + (piC/piY) * exp21; //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double TN93::dPij_dt(int i, int j, double d) const
{
	double kappa1 = _parameters.getParameter("kappa1") -> getValue();
	double kappa2 = _parameters.getParameter("kappa2") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (2. * (piA * piC + piC * piG + piA * piT + piG * piT + kappa2 * piC * piT + kappa1 * piA * piG));
	double l = r * d;
	double k1 = kappa2 * piY + piR;
	double k2 = kappa1 * piR + piY;
	double exp1 = exp(-l);
	double exp22 = exp(-k2 * l);
	double exp21 = exp(-k1 * l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r * (piA * (piY/piR) * exp1 + (piG/piR) * k2 * exp22); //A
				case 1 : return r * (piC *             exp1);                          //C
				case 2 : return r * (piG * (piY/piR) * exp1 - (piG/piR) * k2 * exp22); //G
				case 3 : return r * (piT *             exp1);                          //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return r * (piA *             exp1);                          //A
				case 1 : return r * (piC * (piR/piY) * exp1 + (piT/piY) * k1 * exp21); //C
				case 2 : return r * (piG *             exp1);                          //G
				case 3 : return r * (piT * (piR/piY) * exp1 - (piT/piY) * k1 * exp21); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r * (piA * (piY/piR) * exp1 - (piA/piR) * k2 * exp22); //A
				case 1 : return r * (piC *             exp1);                          //C
				case 2 : return r * (piG * (piY/piR) * exp1 + (piA/piR) * k2 * exp22); //G
				case 3 : return r * (piT *             exp1);                          //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r * (piA *             exp1);                          //A
				case 1 : return r * (piC * (piR/piY) * exp1 - (piC/piY) * k1 * exp21); //C
				case 2 : return r * (piG *             exp1);                          //G
				case 3 : return r * (piT * (piR/piY) * exp1 + (piC/piY) * k1 * exp21); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double TN93::d2Pij_dt2(int i, int j, double d) const
{
	double kappa1 = _parameters.getParameter("kappa1") -> getValue();
	double kappa2 = _parameters.getParameter("kappa2") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (2. * (piA * piC + piC * piG + piA * piT + piG * piT + kappa2 * piC * piT + kappa1 * piA * piG));
	double r_2 = r * r;
	double l = r * d;
	double k1 = kappa2 * piY + piR;
	double k1_2 = k1 * k1;
	double k2 = kappa1 * piR + piY;
	double k2_2 = k2 * k2;
	double exp1 = exp(-l);
	double exp22 = exp(-k2 * l);
	double exp21 = exp(-k1 * l);

	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r_2 * (piA * (piY/piR) * exp1 + (piG/piR) * k2_2 * exp22); //A
				case 1 : return r_2 * (piC *             exp1);                            //C
				case 2 : return r_2 * (piG * (piY/piR) * exp1 - (piG/piR) * k2_2 * exp22); //G
				case 3 : return r_2 * (piT *             exp1);                            //T, U
			}
		}
		//C
		case 1 : {
			switch(j) {
				case 0 : return r_2 * (piA *             exp1);                            //A
				case 1 : return r_2 * (piC * (piR/piY) * exp1 + (piT/piY) * k1_2 * exp21); //C
				case 2 : return r_2 * (piG *             exp1);                            //G
				case 3 : return r_2 * (piT * (piR/piY) * exp1 - (piT/piY) * k1_2 * exp21); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r_2 * (piA * (piY/piR) * exp1 - (piA/piR) * k2_2 * exp22); //A
				case 1 : return r_2 * (piC *             exp1);                            //C
				case 2 : return r_2 * (piG * (piY/piR) * exp1 + (piA/piR) * k2_2 * exp22); //G
				case 3 : return r_2 * (piT *             exp1);                            //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r_2 * (piA *             exp1);                            //A
				case 1 : return r_2 * (piC * (piR/piY) * exp1 - (piC/piY) * k1_2 * exp21); //C
				case 2 : return r_2 * (piG *             exp1);                            //G
				case 3 : return r_2 * (piT * (piR/piY) * exp1 + (piC/piY) * k1_2 * exp21); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

Matrix TN93::getPij_t(double d) const {
	Matrix p(_size, _size);
	double kappa1 = _parameters.getParameter("kappa1") -> getValue();
	double kappa2 = _parameters.getParameter("kappa2") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (2. * (piA * piC + piC * piG + piA * piT + piG * piT + kappa2 * piC * piT + kappa1 * piA * piG));
	double l = r * d;
	double k1 = kappa2 * piY + piR;
	double k2 = kappa1 * piR + piY;
	double exp1 = exp(-l);
	double exp22 = exp(-k2 * l);
	double exp21 = exp(-k1 * l);

	//A
	p(0, 0) = piA * (1. + (piY/piR) * exp1) + (piG/piR) * exp22; //A
	p(0, 1) = piC * (1. -             exp1);                     //C
	p(0, 2) = piG * (1. + (piY/piR) * exp1) - (piG/piR) * exp22; //G
	p(0, 3) = piT * (1. -             exp1);                     //T, U

	//C
	p(1, 0) = piA * (1. -             exp1);                     //A
	p(1, 1) = piC * (1. + (piR/piY) * exp1) + (piT/piY) * exp21; //C
	p(1, 2) = piG * (1. -             exp1);                     //G
	p(1, 3) = piT * (1. + (piR/piY) * exp1) - (piT/piY) * exp21; //T, U

	//G
	p(2, 0) = piA * (1. + (piY/piR) * exp1) - (piA/piR) * exp22; //A
	p(2, 1) = piC * (1. -             exp1);                     //C
	p(2, 2) = piG * (1. + (piY/piR) * exp1) + (piA/piR) * exp22; //G
	p(2, 3) = piT * (1. -             exp1);                     //T, U

	//T, U
	p(3, 0) = piA * (1. -             exp1);                     //A
	p(3, 1) = piC * (1. + (piR/piY) * exp1) - (piC/piY) * exp21; //C
	p(3, 2) = piG * (1. -             exp1);                     //G
	p(3, 3) = piT * (1. + (piR/piY) * exp1) + (piC/piY) * exp21; //T, U

	return p;
}

Matrix TN93::getdPij_dt(double d) const {
	Matrix p(_size, _size);
	double kappa1 = _parameters.getParameter("kappa1") -> getValue();
	double kappa2 = _parameters.getParameter("kappa2") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (2. * (piA * piC + piC * piG + piA * piT + piG * piT + kappa2 * piC * piT + kappa1 * piA * piG));
	double l = r * d;
	double k1 = kappa2 * piY + piR;
	double k2 = kappa1 * piR + piY;
	double exp1 = exp(-l);
	double exp22 = exp(-k2 * l);
	double exp21 = exp(-k1 * l);

	//A
	p(0, 0) = r * (piA * (piY/piR) * exp1 + (piG/piR) * k2 * exp22); //A
	p(0, 1) = r * (piC *             exp1);                          //C
	p(0, 2) = r * (piG * (piY/piR) * exp1 - (piG/piR) * k2 * exp22); //G
	p(0, 3) = r * (piT *             exp1);                          //T, U

	//C
	p(1, 0) = r * (piA *             exp1);                          //A
	p(1, 1) = r * (piC * (piR/piY) * exp1 + (piT/piY) * k1 * exp21); //C
	p(1, 2) = r * (piG *             exp1);                          //G
	p(1, 3) = r * (piT * (piR/piY) * exp1 - (piT/piY) * k1 * exp21); //T, U

	//G
	p(2, 0) = r * (piA * (piY/piR) * exp1 - (piA/piR) * k2 * exp22); //A
	p(2, 1) = r * (piC *             exp1);                          //C
	p(2, 2) = r * (piG * (piY/piR) * exp1 + (piA/piR) * k2 * exp22); //G
	p(2, 3) = r * (piT *             exp1);                          //T, U

	//T, U
	p(3, 0) = r * (piA *             exp1);                          //A
	p(3, 1) = r * (piC * (piR/piY) * exp1 - (piC/piY) * k1 * exp21); //C
	p(3, 2) = r * (piG *             exp1);                          //G
	p(3, 3) = r * (piT * (piR/piY) * exp1 + (piC/piY) * k1 * exp21); //T, U

	return p;
}

Matrix TN93::getd2Pij_dt2(double d) const {
	Matrix p(_size, _size);
	double kappa1 = _parameters.getParameter("kappa1") -> getValue();
	double kappa2 = _parameters.getParameter("kappa2") -> getValue();
	double piA = _parameters.getParameter("piA") -> getValue();
	double piC = _parameters.getParameter("piC") -> getValue();
	double piG = _parameters.getParameter("piG") -> getValue();
	double piT = _parameters.getParameter("piT") -> getValue();
	double piR = piA + piG;
	double piY = piT + piC;
	double r = 1. / (2. * (piA * piC + piC * piG + piA * piT + piG * piT + kappa2 * piC * piT + kappa1 * piA * piG));
	double r_2 = r * r;
	double l = r * d;
	double k1 = kappa2 * piY + piR;
	double k1_2 = k1 * k1;
	double k2 = kappa1 * piR + piY;
	double k2_2 = k2 * k2;
	double exp1 = exp(-l);
	double exp22 = exp(-k2 * l);
	double exp21 = exp(-k1 * l);

	//A
	p(0, 0) = r_2 * (piA * (piY/piR) * exp1 + (piG/piR) * k2_2 * exp22); //A
	p(0, 1) = r_2 * (piC *             exp1);                            //C
	p(0, 2) = r_2 * (piG * (piY/piR) * exp1 - (piG/piR) * k2_2 * exp22); //G
	p(0, 3) = r_2 * (piT *             exp1);                            //T, U

	//C
	p(1, 0) = r_2 * (piA *             exp1);                            //A
	p(1, 1) = r_2 * (piC * (piR/piY) * exp1 + (piT/piY) * k1_2 * exp21); //C
	p(1, 2) = r_2 * (piG *             exp1);                            //G
	p(1, 3) = r_2 * (piT * (piR/piY) * exp1 - (piT/piY) * k1_2 * exp21); //T, U

	//G
	p(2, 0) = r_2 * (piA * (piY/piR) * exp1 - (piA/piR) * k2_2 * exp22); //A
	p(2, 1) = r_2 * (piC *             exp1);                            //C
	p(2, 2) = r_2 * (piG * (piY/piR) * exp1 + (piA/piR) * k2_2 * exp22); //G
	p(2, 3) = r_2 * (piT *             exp1);                            //T, U

	//T, U
	p(3, 0) = r_2 * (piA *             exp1);                            //A
	p(3, 1) = r_2 * (piC * (piR/piY) * exp1 - (piC/piY) * k1_2 * exp21); //C
	p(3, 2) = r_2 * (piG *             exp1);                            //G
	p(3, 3) = r_2 * (piT * (piR/piY) * exp1 + (piC/piY) * k1_2 * exp21); //T, U

	return p;
}

/******************************************************************************/

string TN93::getName() const { return string("Tamura and Nei (1993)"); }

/******************************************************************************/

void TN93::setFreqFromData(const SequenceContainer & data) {
	AbstractSubstitutionModel::setFreqFromData(data);
	setParameterValue("piA", _freq[0]);
	setParameterValue("piC", _freq[1]);
	setParameterValue("piG", _freq[2]);
	setParameterValue("piT", _freq[3]);
}

/******************************************************************************/
