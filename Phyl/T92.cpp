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

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/******************************************************************************/

T92::T92(const NucleicAlphabet * alpha, double kappa, double theta):
	NucleotideSubstitutionModel(alpha, "T92.")
{
	Parameter kappaP("T92.kappa", kappa, &Parameter::R_PLUS_STAR);
	addParameter_(kappaP);
  Parameter thetaP("T92.theta", theta, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(thetaP);
  _p.resize(size_, size_);
	updateMatrices();
}

/******************************************************************************/

void T92::updateMatrices()
{
	_kappa = getParameterValue("kappa");
	_theta = getParameterValue("theta");
  _piA = (1 - _theta) / 2;
	_piC = _theta / 2;
	_piG = _theta / 2;
	_piT = (1 - _theta) / 2;
	_k = (_kappa + 1.) / 2.;
	_r = 2./ (1. + 2. * _theta * _kappa - 2. * _theta * _theta * _kappa);

  freq_[0] = _piA;
	freq_[1] = _piC;
	freq_[2] = _piG;
	freq_[3] = _piT;
	
	generator_(0, 0) = -(1. +        _theta * _kappa)/2;
	generator_(1, 1) = -(1. + (1. - _theta) * _kappa)/2;
	generator_(2, 2) = -(1. + (1. - _theta) * _kappa)/2;
	generator_(3, 3) = -(1. +        _theta * _kappa)/2;

	generator_(1, 0) = (1. - _theta)/2;
	generator_(3, 0) = (1. - _theta)/2;
	generator_(0, 1) = _theta/2;
	generator_(2, 1) = _theta/2;
	generator_(1, 2) = _theta/2;
	generator_(3, 2) = _theta/2;
	generator_(0, 3) = (1. - _theta)/2;
	generator_(2, 3) = (1. - _theta)/2;
	
	generator_(2, 0) = _kappa * (1. - _theta)/2;
	generator_(3, 1) = _kappa * _theta/2;
	generator_(0, 2) = _kappa * _theta/2;
	generator_(1, 3) = _kappa * (1. - _theta)/2;
	
	// Normalization:
	MatrixTools::scale(generator_, _r);

	// Exchangeability:
	exchangeability_(0,0) = generator_(0,0) * 2./(1. - _theta);
	exchangeability_(0,1) = generator_(0,1) * 2./_theta; 
	exchangeability_(0,2) = generator_(0,2) * 2./_theta; 
	exchangeability_(0,3) = generator_(0,3) * 2./(1. - _theta);

	exchangeability_(1,0) = generator_(1,0) * 2./(1. - _theta); 
	exchangeability_(1,1) = generator_(1,1) * 2/_theta; 
	exchangeability_(1,2) = generator_(1,2) * 2/_theta; 
	exchangeability_(1,3) = generator_(1,3) * 2/(1. - _theta); 
	
	exchangeability_(2,0) = generator_(2,0) * 2./(1. - _theta); 
	exchangeability_(2,1) = generator_(2,1) * 2/_theta; 
	exchangeability_(2,2) = generator_(2,2) * 2/_theta; 
	exchangeability_(2,3) = generator_(2,3) * 2/(1. - _theta); 
	
	exchangeability_(3,0) = generator_(3,0) * 2./(1. - _theta);
	exchangeability_(3,1) = generator_(3,1) * 2./_theta; 
	exchangeability_(3,2) = generator_(3,2) * 2./_theta; 
	exchangeability_(3,3) = generator_(3,3) * 2./(1. - _theta);

	// Eigen values:
	eigenValues_[0] = 0;
	eigenValues_[1] = eigenValues_[2] = -_r * (1. + _kappa)/2; 
	eigenValues_[3] = -_r;
	
	// Eigen vectors:
	leftEigenVectors_(0,0) = - (_theta - 1.)/2.;
	leftEigenVectors_(0,1) = _theta/2.;
	leftEigenVectors_(0,2) = _theta/2.;
	leftEigenVectors_(0,3) = - (_theta - 1.)/2.;
	
	leftEigenVectors_(1,0) = 0.;
	leftEigenVectors_(1,1) = - (_theta - 1.);
	leftEigenVectors_(1,2) = 0.;
	leftEigenVectors_(1,3) = _theta - 1.;
	
	leftEigenVectors_(2,0) = _theta;
	leftEigenVectors_(2,1) = 0.;
	leftEigenVectors_(2,2) = -_theta;
	leftEigenVectors_(2,3) = 0.;
	
	leftEigenVectors_(3,0) = - (_theta - 1.)/2.;
	leftEigenVectors_(3,1) = - _theta/2.;
	leftEigenVectors_(3,2) = _theta/2.;
	leftEigenVectors_(3,3) = (_theta - 1.)/2.;


	rightEigenVectors_(0,0) = 1.;
	rightEigenVectors_(0,1) = 0.;
	rightEigenVectors_(0,2) = 1.;
	rightEigenVectors_(0,3) = 1.;
	
	rightEigenVectors_(1,0) = 1.;
	rightEigenVectors_(1,1) = 1.;
	rightEigenVectors_(1,2) = 0.;
	rightEigenVectors_(1,3) = -1.;

	rightEigenVectors_(2,0) = 1.;
	rightEigenVectors_(2,1) = 0.;
	rightEigenVectors_(2,2) = (_theta-1.)/_theta;
	rightEigenVectors_(2,3) = 1.;
	
	rightEigenVectors_(3,0) = 1.;
	rightEigenVectors_(3,1) = _theta/(_theta - 1.);
	rightEigenVectors_(3,2) = 0;
	rightEigenVectors_(3,3) = -1.;
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

const Matrix<double> & T92::getPij_t(double d) const
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

const Matrix<double> & T92::getdPij_dt(double d) const
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

const Matrix<double> & T92::getd2Pij_dt2(double d) const
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

void T92::setFreqFromData(const SequenceContainer & data, unsigned int pseudoCount)
{
	map<int, double> freqs = SequenceContainerTools::getFrequencies(data);
	double f = (freqs[1] + freqs[2] + 2 * pseudoCount) / (freqs[0] + freqs[1] + freqs[2] + freqs[3] + 4 * pseudoCount);
	setParameterValue("theta", f);
	updateMatrices();
}

/******************************************************************************/

