//
// File: K80.cpp
// Created by: Julien Dutheil
// Created on: Tue May 27 15:24:30 2003
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

#include "K80.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <NumCalc/MatrixTools.h>

using namespace std;

/******************************************************************************/

K80::K80(const NucleicAlphabet * alpha, double kappa) :
  //AbstractSubstitutionModel(alpha),
  NucleotideSubstitutionModel(alpha)
{
	_parameters.addParameter(Parameter("kappa", kappa, &Parameter::R_PLUS_STAR));
  _p.resize(_size, _size);
	updateMatrices();
}

/******************************************************************************/

void K80::updateMatrices()
{
	_kappa = _parameters.getParameter("kappa")->getValue();
	_k = (_kappa + 1.) / 2.;
	_r = 4. / (_kappa + 2.);
	
  // Frequences:
	_freq[0] = _freq[1] = _freq[2] = _freq[3] = 1. / 4.;

	// Generator:
	_generator(0, 0) = -2. - _kappa;
	_generator(1, 1) = -2. - _kappa;
	_generator(2, 2) = -2. - _kappa;
	_generator(3, 3) = -2. - _kappa;

	_generator(0, 1) = 1.;
	_generator(0, 3) = 1.;
	_generator(1, 0) = 1.;
	_generator(1, 2) = 1.;
	_generator(2, 1) = 1.;
	_generator(2, 3) = 1.;
	_generator(3, 0) = 1.;
	_generator(3, 2) = 1.;
	
	_generator(0, 2) = _kappa;
	_generator(1, 3) = _kappa;
	_generator(2, 0) = _kappa;
	_generator(3, 1) = _kappa;

	// Normalization:
	MatrixTools::scale(_generator, _r/4);

	// Exchangeability:
	_exchangeability = _generator;
	MatrixTools::scale(_exchangeability, 4.);

	// Eigen values:
	_eigenValues[0] = 0;
	_eigenValues[1] = -_r * (1. + _kappa)/2;
	_eigenValues[2] = -_r * (1. + _kappa)/2;
	_eigenValues[3] = -_r;
	
	// Eigen vectors:
	_leftEigenVectors(0,0) = 1. / 4.;
	_leftEigenVectors(0,1) = 1. / 4.;
	_leftEigenVectors(0,2) = 1. / 4.;
	_leftEigenVectors(0,3) = 1. / 4.;
	_leftEigenVectors(1,0) = 0.;
	_leftEigenVectors(1,1) = 1. / 2.;
	_leftEigenVectors(1,2) = 0.;
	_leftEigenVectors(1,3) = -1. / 2.;
	_leftEigenVectors(2,0) = 1. / 2.;
	_leftEigenVectors(2,1) = 0.;
	_leftEigenVectors(2,2) = -1. / 2.;
	_leftEigenVectors(2,3) = 0.;
	_leftEigenVectors(3,0) = 1. / 4.;
	_leftEigenVectors(3,1) = -1. / 4.;
	_leftEigenVectors(3,2) = 1. / 4.;
	_leftEigenVectors(3,3) = -1. / 4.;

	_rightEigenVectors(0,0) = 1.;
	_rightEigenVectors(0,1) = 0.;
	_rightEigenVectors(0,2) = 1.;
	_rightEigenVectors(0,3) = 1.;
	_rightEigenVectors(1,0) = 1.;
	_rightEigenVectors(1,1) = 1.;
	_rightEigenVectors(1,2) = 0.;
	_rightEigenVectors(1,3) = -1.;
	_rightEigenVectors(2,0) = 1.;
	_rightEigenVectors(2,1) = 0.;
	_rightEigenVectors(2,2) = -1.;
	_rightEigenVectors(2,3) = 1.;
	_rightEigenVectors(3,0) = 1.;
	_rightEigenVectors(3,1) = -1.;
	_rightEigenVectors(3,2) = 0;
	_rightEigenVectors(3,3) = -1.;
}
	
/******************************************************************************/

double K80::Pij_t(int i, int j, double d) const
{
	_l = _r * d;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return 0.25 * (1. + _exp1) + 0.5 * _exp2; //A
				case 1 : return 0.25 * (1. - _exp1);               //C
				case 2 : return 0.25 * (1. + _exp1) - 0.5 * _exp2; //G
				case 3 : return 0.25 * (1. - _exp1);               //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return 0.25 * (1. - _exp1);               //A
				case 1 : return 0.25 * (1. + _exp1) + 0.5 * _exp2; //C
				case 2 : return 0.25 * (1. - _exp1);               //G
				case 3 : return 0.25 * (1. + _exp1) - 0.5 * _exp2; //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return 0.25 * (1. + _exp1) - 0.5 * _exp2; //A
				case 1 : return 0.25 * (1. - _exp1);               //C
				case 2 : return 0.25 * (1. + _exp1) + 0.5 * _exp2; //G
				case 3 : return 0.25 * (1. - _exp1);               //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return 0.25 * (1. - _exp1);               //A
				case 1 : return 0.25 * (1. + _exp1) - 0.5 * _exp2; //C
				case 2 : return 0.25 * (1. - _exp1);               //G
				case 3 : return 0.25 * (1. + _exp1) + 0.5 * _exp2; //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double K80::dPij_dt(int i, int j, double d) const
{
	_l = _r * d;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);

	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return _r/4. * (- _exp1 - 2. * _k * _exp2); //A
				case 1 : return _r/4. * (  _exp1);                   //C
				case 2 : return _r/4. * (- _exp1 + 2. * _k * _exp2); //G
				case 3 : return _r/4. * (  _exp1);                   //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return _r/4. * (  _exp1);                   //A
				case 1 : return _r/4. * (- _exp1 - 2. * _k * _exp2); //C
				case 2 : return _r/4. * (  _exp1);                   //G
				case 3 : return _r/4. * (- _exp1 + 2. * _k * _exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return _r/4. * (- _exp1 + 2. * _k * _exp2); //A
				case 1 : return _r/4. * (  _exp1);                   //C
				case 2 : return _r/4. * (- _exp1 - 2. * _k * _exp2); //G
				case 3 : return _r/4. * (  _exp1);                   //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return _r/4. * (  _exp1);                   //A
				case 1 : return _r/4. * (- _exp1 + 2. * _k * _exp2); //C
				case 2 : return _r/4. * (  _exp1);                   //G
				case 3 : return _r/4. * (- _exp1 - 2. * _k * _exp2); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double K80::d2Pij_dt2(int i, int j, double d) const
{
	double k_2 = _k * _k;
	double r_2 = _r * _r;
	_l = _r * d;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);

	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r_2/4. * (  _exp1 + 2. * k_2 * _exp2); //A
				case 1 : return r_2/4. * (- _exp1);                    //C
				case 2 : return r_2/4. * (  _exp1 - 2. * k_2 * _exp2); //G
				case 3 : return r_2/4. * (- _exp1);                    //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return r_2/4. * (- _exp1);                    //A
				case 1 : return r_2/4. * (  _exp1 + 2. * k_2 * _exp2); //C
				case 2 : return r_2/4. * (- _exp1);                    //G
				case 3 : return r_2/4. * (  _exp1 - 2. * k_2 * _exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r_2/4. * (  _exp1 - 2. * k_2 * _exp2); //A
				case 1 : return r_2/4. * (- _exp1);                    //C
				case 2 : return r_2/4. * (  _exp1 + 2. * k_2 * _exp2); //G
				case 3 : return r_2/4. * (- _exp1);                    //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r_2/4. * (- _exp1);                    //A
				case 1 : return r_2/4. * (  _exp1 - 2. * k_2 * _exp2); //C
				case 2 : return r_2/4. * (- _exp1);                    //G
				case 3 : return r_2/4. * (  _exp1 + 2. * k_2 * _exp2); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

RowMatrix<double> K80::getPij_t(double d) const
{
	_l = _r * d;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);

	//A
	_p(0, 0) = 0.25 * (1. + _exp1) + 0.5 * _exp2; //A
	_p(0, 1) = 0.25 * (1. - _exp1);               //C
	_p(0, 2) = 0.25 * (1. + _exp1) - 0.5 * _exp2; //G
	_p(0, 3) = 0.25 * (1. - _exp1);               //T, U

	//C
	_p(1, 0) = 0.25 * (1. - _exp1);               //A
	_p(1, 1) = 0.25 * (1. + _exp1) + 0.5 * _exp2; //C
	_p(1, 2) = 0.25 * (1. - _exp1);               //G
	_p(1, 3) = 0.25 * (1. + _exp1) - 0.5 * _exp2; //T, U

	//G
	_p(2, 0) = 0.25 * (1. + _exp1) - 0.5 * _exp2; //A
	_p(2, 1) = 0.25 * (1. - _exp1);               //C
	_p(2, 2) = 0.25 * (1. + _exp1) + 0.5 * _exp2; //G
	_p(2, 3) = 0.25 * (1. - _exp1);               //T, U

	//T, U
	_p(3, 0) = 0.25 * (1. - _exp1);               //A
	_p(3, 1) = 0.25 * (1. + _exp1) - 0.5 * _exp2; //C
	_p(3, 2) = 0.25 * (1. - _exp1);               //G
	_p(3, 3) = 0.25 * (1. + _exp1) + 0.5 * _exp2; //T, U

	return _p;
}

RowMatrix<double> K80::getdPij_dt(double d) const
{
	_l = _r * d;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);

	_p(0, 0) = _r/4. * (- _exp1 - 2. * _k * _exp2); //A
	_p(0, 1) = _r/4. * (  _exp1);                   //C
	_p(0, 2) = _r/4. * (- _exp1 + 2. * _k * _exp2); //G
	_p(0, 3) = _r/4. * (  _exp1);                   //T, U

	//C
	_p(1, 0) = _r/4. * (  _exp1);                   //A
	_p(1, 1) = _r/4. * (- _exp1 - 2. * _k * _exp2); //C
	_p(1, 2) = _r/4. * (  _exp1);                   //G
	_p(1, 3) = _r/4. * (- _exp1 + 2. * _k * _exp2); //T, U

	//G
	_p(2, 0) = _r/4. * (- _exp1 + 2. * _k * _exp2); //A
	_p(2, 1) = _r/4. * (  _exp1);                   //C
	_p(2, 2) = _r/4. * (- _exp1 - 2. * _k * _exp2); //G
	_p(2, 3) = _r/4. * (  _exp1);                   //T, U

	//T, U
	_p(3, 0) = _r/4. * (  _exp1);                   //A
	_p(3, 1) = _r/4. * (- _exp1 + 2. * _k * _exp2); //C
	_p(3, 2) = _r/4. * (  _exp1);                   //G
	_p(3, 3) = _r/4. * (- _exp1 - 2. * _k * _exp2); //T, U

	return _p;
}

RowMatrix<double> K80::getd2Pij_dt2(double d) const
{
	double k_2 = _k * _k;
	double r_2 = _r * _r;
	_l = _r * d;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);

	_p(0, 0) = r_2/4. * (  _exp1 + 2. * k_2 * _exp2); //A
	_p(0, 1) = r_2/4. * (- _exp1);                    //C
	_p(0, 2) = r_2/4. * (  _exp1 - 2. * k_2 * _exp2); //G
	_p(0, 3) = r_2/4. * (- _exp1);                    //T, U

	//C
	_p(1, 0) = r_2/4. * (- _exp1);                    //A
	_p(1, 1) = r_2/4. * (  _exp1 + 2. * k_2 * _exp2); //C
	_p(1, 2) = r_2/4. * (- _exp1);                    //G
	_p(1, 3) = r_2/4. * (  _exp1 - 2. * k_2 * _exp2); //T, U

	//G
	_p(2, 0) = r_2/4. * (  _exp1 - 2. * k_2 * _exp2); //A
	_p(2, 1) = r_2/4. * (- _exp1);                    //C
	_p(2, 2) = r_2/4. * (  _exp1 + 2. * k_2 * _exp2); //G
	_p(2, 3) = r_2/4. * (- _exp1);                    //T, U

	//T, U
	_p(3, 0) = r_2/4. * (- _exp1);                    //A
	_p(3, 1) = r_2/4. * (  _exp1 - 2. * k_2 * _exp2); //C
	_p(3, 2) = r_2/4. * (- _exp1);                    //G
	_p(3, 3) = r_2/4. * (  _exp1 + 2. * k_2 * _exp2); //T, U

	return _p;
}

/******************************************************************************/

