//
// File: K80.h
// Created by:  <@bogdanof>
// Created on: Tue May 27 15:24:30 2003
//

#include "K80.h"

// From the STL:
#include <cmath>

// From the MTL:
#include <mtl/mtl.h>

/******************************************************************************/

K80::K80(const NucleicAlphabet * alpha, double kappa): NucleotideSubstitutionModel(alpha)
{
	_parameters.addParameter(Parameter("kappa", kappa, &Parameter::R_PLUS));
	updateMatrices();
}

K80::~K80() {}
	
/******************************************************************************/

void K80::updateMatrices() {
	double kappa = _parameters.getParameter("kappa") -> getValue();
	
	// Generator:
	_generator(0, 0) = -1;
	_generator(1, 1) = -1;
	_generator(2, 2) = -1;
	_generator(3, 3) = -1;

	_generator(0, 1) = 1.;
	_generator(0, 3) = 1.;
	_generator(1, 0) = 1.;
	_generator(1, 2) = 1.;
	_generator(2, 1) = 1.;
	_generator(2, 3) = 1.;
	_generator(3, 0) = 1.;
	_generator(3, 2) = 1.;
	
	_generator(0, 2) = kappa;
	_generator(1, 3) = kappa;
	_generator(2, 0) = kappa;
	_generator(3, 1) = kappa;
	
	// Normalization:
	double r = 1. / 2. + kappa;
	scale(_generator, r);

	// Frequences:
	_freq[0] = _freq[1] = _freq[2] = _freq[3] = 1. / 4.;

	// Eigen values:
	_eigenValues[0] = 0;
	_eigenValues[1] = -4. * r;
	_eigenValues[2] = - 2. * (1. + kappa) * r;
	
	// Eigen vectors:
	// todo!
}
	
/******************************************************************************/

double K80::Pij_t(int i, int j, double d) const {
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double k = (kappa + 1.) / 2.;
	double r = 4. / (kappa + 2.);
	double l = r * d;
	double exp1 = exp(-l);
	double exp2 = exp(-k * l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return 0.25 * (1. + exp1) + 0.5 * exp2; //A
				case 1 : return 0.25 * (1. - exp1);              //C
				case 2 : return 0.25 * (1. + exp1) - 0.5 * exp2; //G
				case 3 : return 0.25 * (1. - exp1);              //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return 0.25 * (1. - exp1);              //A
				case 1 : return 0.25 * (1. + exp1) + 0.5 * exp2; //C
				case 2 : return 0.25 * (1. - exp1);              //G
				case 3 : return 0.25 * (1. + exp1) - 0.5 * exp2; //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return 0.25 * (1. + exp1) - 0.5 * exp2; //A
				case 1 : return 0.25 * (1. - exp1);              //C
				case 2 : return 0.25 * (1. + exp1) + 0.5 * exp2; //G
				case 3 : return 0.25 * (1. - exp1);              //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return 0.25 * (1. - exp1);              //A
				case 1 : return 0.25 * (1. + exp1) - 0.5 * exp2; //C
				case 2 : return 0.25 * (1. - exp1);              //G
				case 3 : return 0.25 * (1. + exp1) + 0.5 * exp2; //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double K80::dPij_dt(int i, int j, double d) const {
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double k = (kappa + 1.) / 2.;
	double r = 4. / (kappa + 2.);
	double l = r * d;
	double exp1 = exp(-l);
	double exp2 = exp(-k * l);

	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r * (- exp1 - 0.5 * k * exp2); //A
				case 1 : return r * (  exp1);                  //C
				case 2 : return r * (- exp1 + 0.5 * k * exp2); //G
				case 3 : return r * (  exp1);                  //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return r * (  exp1);                  //A
				case 1 : return r * (- exp1 - 0.5 * k * exp2); //C
				case 2 : return r * (  exp1);                  //G
				case 3 : return r * (- exp1 + 0.5 * k * exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r * (- exp1 + 0.5 * k * exp2); //A
				case 1 : return r * (  exp1);                  //C
				case 2 : return r * (- exp1 - 0.5 * k * exp2); //G
				case 3 : return r * (  exp1);                  //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r * (  exp1);                  //A
				case 1 : return r * (- exp1 + 0.5 * k * exp2); //C
				case 2 : return r * (  exp1);                  //G
				case 3 : return r * (- exp1 - 0.5 * k * exp2); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double K80::d2Pij_dt2(int i, int j, double d) const {
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double k = (kappa + 1.) / 2.;
	double k_2 = k * k;
	double r = 4. / (kappa + 2.);
	double r_2 = r * r;
	double l = r * d;
	double exp1 = exp(-l);
	double exp2 = exp(-k * l);

	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r_2 * (  exp1 + 0.5 * k_2 * exp2); //A
				case 1 : return r_2 * (- exp1);                    //C
				case 2 : return r_2 * (  exp1 - 0.5 * k_2 * exp2); //G
				case 3 : return r_2 * (- exp1);                    //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return r_2 * (- exp1);                    //A
				case 1 : return r_2 * (  exp1 + 0.5 * k_2 * exp2); //C
				case 2 : return r_2 * (- exp1);                    //G
				case 3 : return r_2 * (  exp1 - 0.5 * k_2 * exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r_2 * (  exp1 - 0.5 * k_2 * exp2); //A
				case 1 : return r_2 * (- exp1);                    //C
				case 2 : return r_2 * (  exp1 + 0.5 * k_2 * exp2); //G
				case 3 : return r_2 * (- exp1);                    //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r_2 * (- exp1);                    //A
				case 1 : return r_2 * (  exp1 - 0.5 * k_2 * exp2); //C
				case 2 : return r_2 * (- exp1);                    //G
				case 3 : return r_2 * (  exp1 + 0.5 * k_2 * exp2); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

Matrix K80::getPij_t(double d) const {
	Matrix p(_size, _size);
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double k = (kappa + 1.) / 2.;
	double r = 4. / (kappa + 2.);
	double l = r * d;
	double exp1 = exp(-l);
	double exp2 = exp(-k * l);

	//A
	p(0, 0) = 0.25 * (1. + exp1) + 0.5 * exp2; //A
	p(0, 1) = 0.25 * (1. - exp1);              //C
	p(0, 2) = 0.25 * (1. + exp1) - 0.5 * exp2; //G
	p(0, 3) = 0.25 * (1. - exp1);              //T, U

	//C
	p(1, 0) = 0.25 * (1. - exp1);              //A
	p(1, 1) = 0.25 * (1. + exp1) + 0.5 * exp2; //C
	p(1, 2) = 0.25 * (1. - exp1);              //G
	p(1, 3) = 0.25 * (1. + exp1) - 0.5 * exp2; //T, U

	//G
	p(2, 0) = 0.25 * (1. + exp1) - 0.5 * exp2; //A
	p(2, 1) = 0.25 * (1. - exp1);              //C
	p(2, 2) = 0.25 * (1. + exp1) + 0.5 * exp2; //G
	p(2, 3) = 0.25 * (1. - exp1);              //T, U

	//T, U
	p(3, 0) = 0.25 * (1. - exp1);              //A
	p(3, 1) = 0.25 * (1. + exp1) - 0.5 * exp2; //C
	p(3, 2) = 0.25 * (1. - exp1);              //G
	p(3, 3) = 0.25 * (1. + exp1) + 0.5 * exp2; //T, U

	return p;
}

Matrix K80::getdPij_dt(double d) const {
	Matrix p(_size, _size);
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double k = (kappa + 1.) / 2.;
	double r = 4. / (kappa + 2.);
	double l = r * d;
	double exp1 = exp(-l);
	double exp2 = exp(-k * l);

	p(0, 0) = r * (- exp1 - 0.5 * k * exp2); //A
	p(0, 1) = r * (  exp1);                  //C
	p(0, 2) = r * (- exp1 + 0.5 * k * exp2); //G
	p(0, 3) = r * (  exp1);                  //T, U

	//C
	p(1, 0) = r * (  exp1);                  //A
	p(1, 1) = r * (- exp1 - 0.5 * k * exp2); //C
	p(1, 2) = r * (  exp1);                  //G
	p(1, 3) = r * (- exp1 + 0.5 * k * exp2); //T, U

	//G
	p(2, 0) = r * (- exp1 + 0.5 * k * exp2); //A
	p(2, 1) = r * (  exp1);                  //C
	p(2, 2) = r * (- exp1 - 0.5 * k * exp2); //G
	p(2, 3) = r * (  exp1);                  //T, U

	//T, U
	p(3, 0) = r * (  exp1);                  //A
	p(3, 1) = r * (- exp1 + 0.5 * k * exp2); //C
	p(3, 2) = r * (  exp1);                  //G
	p(3, 3) = r * (- exp1 - 0.5 * k * exp2); //T, U

	return p;
}

Matrix K80::getd2Pij_dt2(double d) const {
	Matrix p(_size, _size);
	double kappa = _parameters.getParameter("kappa") -> getValue();
	double k = (kappa + 1.) / 2.;
	double k_2 = k * k;
	double r = 4. / (kappa + 2.);
	double r_2 = r * r;
	double l = r * d;
	double exp1 = exp(-l);
	double exp2 = exp(-k * l);

	p(0, 0) = r_2 * (  exp1 + 0.5 * k_2 * exp2); //A
	p(0, 1) = r_2 * (- exp1);                    //C
	p(0, 2) = r_2 * (  exp1 - 0.5 * k_2 * exp2); //G
	p(0, 3) = r_2 * (- exp1);                    //T, U

	//C
	p(1, 0) = r_2 * (- exp1);                    //A
	p(1, 1) = r_2 * (  exp1 + 0.5 * k_2 * exp2); //C
	p(1, 2) = r_2 * (- exp1);                    //G
	p(1, 3) = r_2 * (  exp1 - 0.5 * k_2 * exp2); //T, U

	//G
	p(2, 0) = r_2 * (  exp1 - 0.5 * k_2 * exp2); //A
	p(2, 1) = r_2 * (- exp1);                    //C
	p(2, 2) = r_2 * (  exp1 + 0.5 * k_2 * exp2); //G
	p(2, 3) = r_2 * (- exp1);                    //T, U

	//T, U
	p(3, 0) = r_2 * (- exp1);                    //A
	p(3, 1) = r_2 * (  exp1 - 0.5 * k_2 * exp2); //C
	p(3, 2) = r_2 * (- exp1);                    //G
	p(3, 3) = r_2 * (  exp1 + 0.5 * k_2 * exp2); //T, U

	return p;
}

/******************************************************************************/

string K80::getName() const { return string("Kimura (1980)"); }

/******************************************************************************/
