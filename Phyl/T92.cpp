//
// File: T92.cpp
// Created by:  Julien Dutheil
// Created on: Mon May 26 14:41:24 2003
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

#include "T92.h"

// From SeqLib:
#include <Seq/SequenceContainerTools.h>

// From NumCalc:
#include <NumCalc/MatrixTools.h>

// From the STL:
#include <cmath>


/******************************************************************************/

T92::T92(const NucleicAlphabet * alpha, double kappa, double theta):
	//AbstractSubstitutionModel(alpha),
	NucleotideSubstitutionModel(alpha)
{
	thetaConstraint = new IncludingInterval(0, 1);
	_parameters.addParameter(Parameter("kappa", kappa, &Parameter::R_PLUS));
	_parameters.addParameter(Parameter("theta", theta, thetaConstraint));
  _p.resize(_size, _size);
	updateMatrices();
}

/******************************************************************************/

T92::~T92()
{
	delete thetaConstraint;
}

/******************************************************************************/

void T92::updateMatrices()
{
	_kappa = _parameters.getParameter("kappa")->getValue();
	_theta = _parameters.getParameter("theta")->getValue();
  cout << "theta=" << _theta << endl;
  _piA = (1 - _theta) / 2;
	_piC = _theta / 2;
	_piG = _theta / 2;
	_piT = (1 - _theta) / 2;
	_k = (_kappa + 1.) / 2.;
	_r = 2./ (1. + 2. * _theta * _kappa - 2. * _theta * _theta * _kappa);

  _freq[0] = _piA;
	_freq[1] = _piC;
	_freq[2] = _piG;
	_freq[3] = _piT;
	
	_generator(0, 0) = -(1. +        _theta * _kappa)/2;
	_generator(1, 1) = -(1. + (1. - _theta) * _kappa)/2;
	_generator(2, 2) = -(1. + (1. - _theta) * _kappa)/2;
	_generator(3, 3) = -(1. +        _theta * _kappa)/2;

	_generator(1, 0) = (1. - _theta)/2;
	_generator(3, 0) = (1. - _theta)/2;
	_generator(0, 1) = _theta/2;
	_generator(2, 1) = _theta/2;
	_generator(1, 2) = _theta/2;
	_generator(3, 2) = _theta/2;
	_generator(0, 3) = (1. - _theta)/2;
	_generator(2, 3) = (1. - _theta)/2;
	
	_generator(2, 0) = _kappa * (1. - _theta)/2;
	_generator(3, 1) = _kappa * _theta/2;
	_generator(0, 2) = _kappa * _theta/2;
	_generator(1, 3) = _kappa * (1. - _theta)/2;
	
	// Normalization:
	MatrixTools::scale(_generator, _r);

	// Exchangeability:
	_exchangeability(0,0) = _generator(0,0) * 2./(1. - _theta);
	_exchangeability(0,1) = _generator(0,1) * 2./_theta; 
	_exchangeability(0,2) = _generator(0,2) * 2./_theta; 
	_exchangeability(0,3) = _generator(0,3) * 2./(1. - _theta);

	_exchangeability(1,0) = _generator(1,0) * 2./(1. - _theta); 
	_exchangeability(1,1) = _generator(1,1) * 2/_theta; 
	_exchangeability(1,2) = _generator(1,2) * 2/_theta; 
	_exchangeability(1,3) = _generator(1,3) * 2/(1. - _theta); 
	
	_exchangeability(2,0) = _generator(2,0) * 2./(1. - _theta); 
	_exchangeability(2,1) = _generator(2,1) * 2/_theta; 
	_exchangeability(2,2) = _generator(2,2) * 2/_theta; 
	_exchangeability(2,3) = _generator(2,3) * 2/(1. - _theta); 
	
	_exchangeability(3,0) = _generator(3,0) * 2./(1. - _theta);
	_exchangeability(3,1) = _generator(3,1) * 2./_theta; 
	_exchangeability(3,2) = _generator(3,2) * 2./_theta; 
	_exchangeability(3,3) = _generator(3,3) * 2./(1. - _theta);

	// Eigen values:
	_eigenValues[0] = 0;
	_eigenValues[1] = _eigenValues[2] = -_r * (1. + _kappa)/2; 
	_eigenValues[3] = -_r;
	
	// Eigen vectors:
	_leftEigenVectors(0,0) = - (_theta - 1.)/2.;
	_leftEigenVectors(0,1) = _theta/2.;
	_leftEigenVectors(0,2) = _theta/2.;
	_leftEigenVectors(0,3) = - (_theta - 1.)/2.;
	
	_leftEigenVectors(1,0) = 0.;
	_leftEigenVectors(1,1) = - (_theta - 1.);
	_leftEigenVectors(1,2) = 0.;
	_leftEigenVectors(1,3) = _theta - 1.;
	
	_leftEigenVectors(2,0) = _theta;
	_leftEigenVectors(2,1) = 0.;
	_leftEigenVectors(2,2) = -_theta;
	_leftEigenVectors(2,3) = 0.;
	
	_leftEigenVectors(3,0) = - (_theta - 1.)/2.;
	_leftEigenVectors(3,1) = - _theta/2.;
	_leftEigenVectors(3,2) = _theta/2.;
	_leftEigenVectors(3,3) = (_theta - 1.)/2.;


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
	_rightEigenVectors(2,2) = (_theta-1.)/_theta;
	_rightEigenVectors(2,3) = 1.;
	
	_rightEigenVectors(3,0) = 1.;
	_rightEigenVectors(3,1) = _theta/(_theta - 1.);
	_rightEigenVectors(3,2) = 0;
	_rightEigenVectors(3,3) = -1.;
}
	
/******************************************************************************/

double T92::Pij_t(int i, int j, double d) const
{
	_l = _r * d;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return _piA * (1. + _exp1) + _theta * _exp2; //A
				case 1 : return _piC * (1. - _exp1);                  //C
				case 2 : return _piG * (1. + _exp1) - _theta * _exp2; //G
				case 3 : return _piT * (1. - _exp1);                  //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return _piA * (1. - _exp1);                         //A
				case 1 : return _piC * (1. + _exp1) + (1. - _theta) * _exp2; //C
				case 2 : return _piG * (1. - _exp1);                         //G
				case 3 : return _piT * (1. + _exp1) - (1. - _theta) * _exp2; //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return _piA * (1. + _exp1) - (1. - _theta) * _exp2; //A
				case 1 : return _piC * (1. - _exp1);                         //C
				case 2 : return _piG * (1. + _exp1) + (1. - _theta) * _exp2; //G
				case 3 : return _piT * (1. - _exp1);                         //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return _piA * (1. - _exp1);                  //A
				case 1 : return _piC * (1. + _exp1) - _theta * _exp2; //C
				case 2 : return _piG * (1. - _exp1);                  //G
				case 3 : return _piT * (1. + _exp1) + _theta * _exp2; //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double T92::dPij_dt(int i, int j, double d) const
{
	_l = _r * d;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);

	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return _r * (_piA * - _exp1 + _theta * -_k * _exp2); //A
				case 1 : return _r * (_piC *   _exp1);                        //C
				case 2 : return _r * (_piG * - _exp1 - _theta * -_k * _exp2); //G
				case 3 : return _r * (_piT *   _exp1);                        //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return _r * (_piA *   _exp1);                               //A
				case 1 : return _r * (_piC * - _exp1 + (1. - _theta) * -_k * _exp2); //C
				case 2 : return _r * (_piG *   _exp1);                               //G
				case 3 : return _r * (_piT * - _exp1 - (1. - _theta) * -_k * _exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return _r * (_piA * - _exp1 - (1. - _theta) * -_k * _exp2); //A
				case 1 : return _r * (_piC *   _exp1);                               //C
				case 2 : return _r * (_piG * - _exp1 + (1. - _theta) * -_k * _exp2); //G
				case 3 : return _r * (_piT *   _exp1);                               //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return _r * (_piA *   _exp1);                        //A
				case 1 : return _r * (_piC * - _exp1 - _theta * -_k * _exp2); //C
				case 2 : return _r * (_piG *   _exp1);                        //G
				case 3 : return _r * (_piT * - _exp1 + _theta * -_k * _exp2); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double T92::d2Pij_dt2(int i, int j, double d) const
{
	double _k2 = _k * _k;
	_l = _r * d;
	double r2 = _r * _r;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);

	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r2 * (_piA *   _exp1 + _theta * _k2 * _exp2); //A
				case 1 : return r2 * (_piC * - _exp1);                        //C
				case 2 : return r2 * (_piG *   _exp1 - _theta * _k2 * _exp2); //G
				case 3 : return r2 * (_piT * - _exp1);                        //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return r2 * (_piA * - _exp1);                               //A
				case 1 : return r2 * (_piC *   _exp1 + (1. - _theta) * _k2 * _exp2); //C
				case 2 : return r2 * (_piG * - _exp1);                               //G
				case 3 : return r2 * (_piT *   _exp1 - (1. - _theta) * _k2 * _exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r2 * (_piA *   _exp1 - (1. - _theta) * _k2 * _exp2); //A
				case 1 : return r2 * (_piC * - _exp1);                               //C
				case 2 : return r2 * (_piG *   _exp1 + (1. - _theta) * _k2 * _exp2); //G
				case 3 : return r2 * (_piT * - _exp1);                               //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r2 * (_piA * - _exp1);                        //A
				case 1 : return r2 * (_piC *   _exp1 - _theta * _k2 * _exp2); //C
				case 2 : return r2 * (_piG * - _exp1);                        //G
				case 3 : return r2 * (_piT *   _exp1 + _theta * _k2 * _exp2); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

RowMatrix<double> T92::getPij_t(double d) const
{
	_l = _r * d;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);

	//A
	_p(0, 0) = _piA * (1. + _exp1) + _theta * _exp2; //A
	_p(0, 1) = _piC * (1. - _exp1);                  //C
	_p(0, 2) = _piG * (1. + _exp1) - _theta * _exp2; //G
	_p(0, 3) = _piT * (1. - _exp1);                  //T, U

	//C
	_p(1, 0) = _piA * (1. - _exp1);                         //A
	_p(1, 1) = _piC * (1. + _exp1) + (1. - _theta) * _exp2; //C
	_p(1, 2) = _piG * (1. - _exp1);                         //G
	_p(1, 3) = _piT * (1. + _exp1) - (1. - _theta) * _exp2; //T, U

	//G
	_p(2, 0) = _piA * (1. + _exp1) - (1. - _theta) * _exp2; //A
	_p(2, 1) = _piC * (1. - _exp1);                         //C
	_p(2, 2) = _piG * (1. + _exp1) + (1. - _theta) * _exp2; //G
	_p(2, 3) = _piT * (1. - _exp1);                         //T, U

	//T, U
	_p(3, 0) = _piA * (1. - _exp1);                  //A
	_p(3, 1) = _piC * (1. + _exp1) - _theta * _exp2; //C
	_p(3, 2) = _piG * (1. - _exp1);                  //G
	_p(3, 3) = _piT * (1. + _exp1) + _theta * _exp2; //T, U

	return _p;
}

RowMatrix<double> T92::getdPij_dt(double d) const
{
	_l = _r * d;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);

	//A
	_p(0, 0) = _r * (_piA * - _exp1 + _theta * -_k * _exp2); //A
	_p(0, 1) = _r * (_piC *   _exp1);                        //C
	_p(0, 2) = _r * (_piG * - _exp1 - _theta * -_k * _exp2); //G
	_p(0, 3) = _r * (_piT *   _exp1);                        //T, U

	//C
	_p(1, 0) = _r * (_piA *   _exp1);                               //A
	_p(1, 1) = _r * (_piC * - _exp1 + (1. - _theta) * -_k * _exp2); //C
	_p(1, 2) = _r * (_piG *   _exp1);                               //G
	_p(1, 3) = _r * (_piT * - _exp1 - (1. - _theta) * -_k * _exp2); //T, U

	//G
	_p(2, 0) = _r * (_piA * - _exp1 - (1. - _theta) * -_k * _exp2); //A
	_p(2, 1) = _r * (_piC *   _exp1);                               //C
	_p(2, 2) = _r * (_piG * - _exp1 + (1. - _theta) * -_k * _exp2); //G
	_p(2, 3) = _r * (_piT *   _exp1);                               //T, U

	//T, U
	_p(3, 0) = _r * (_piA *   _exp1);                        //A
	_p(3, 1) = _r * (_piC * - _exp1 - _theta * -_k * _exp2); //C
	_p(3, 2) = _r * (_piG *   _exp1);                        //G
	_p(3, 3) = _r * (_piT * - _exp1 + _theta * -_k * _exp2); //T, U

	return _p;
}

RowMatrix<double> T92::getd2Pij_dt2(double d) const
{
	double k2 = _k * _k;
	_l = _r * d;
	double r2 = _r * _r;
	_exp1 = exp(-_l);
	_exp2 = exp(-_k * _l);

	//A
	_p(0, 0) = r2 * (_piA *   _exp1 + _theta * k2 * _exp2); //A
	_p(0, 1) = r2 * (_piC * - _exp1);                       //C
	_p(0, 2) = r2 * (_piG *   _exp1 - _theta * k2 * _exp2); //G
	_p(0, 3) = r2 * (_piT * - _exp1);                       //T, U

	//C
	_p(1, 0) = r2 * (_piA * - _exp1);                              //A
	_p(1, 1) = r2 * (_piC *   _exp1 + (1. - _theta) * k2 * _exp2); //C
	_p(1, 2) = r2 * (_piG * - _exp1);                              //G
	_p(1, 3) = r2 * (_piT *   _exp1 - (1. - _theta) * k2 * _exp2); //T, U

	//G
	_p(2, 0) = r2 * (_piA *   _exp1 - (1. - _theta) * k2 * _exp2); //A
	_p(2, 1) = r2 * (_piC * - _exp1);                              //C
	_p(2, 2) = r2 * (_piG *   _exp1 + (1. - _theta) * k2 * _exp2); //G
	_p(2, 3) = r2 * (_piT * - _exp1);                              //T, U

	//T, U
	_p(3, 0) = r2 * (_piA * - _exp1);                       //A
	_p(3, 1) = r2 * (_piC *   _exp1 - _theta * k2 * _exp2); //C
	_p(3, 2) = r2 * (_piG * - _exp1);                       //G
	_p(3, 3) = r2 * (_piT *   _exp1 + _theta * k2 * _exp2); //T, U

	return _p;
}

/******************************************************************************/

string T92::getName() const { return string("Tamura (1992)"); }

/******************************************************************************/

void T92::setFreqFromData(const SequenceContainer & data)
{
	map<int, double> freqs = SequenceContainerTools::getFrequencies(data);
	double f = (freqs[1] + freqs[2]) / (freqs[0] + freqs[1] + freqs[2] + freqs[3]);
	setParameterValue("theta", f);
	updateMatrices();
}

/******************************************************************************/

