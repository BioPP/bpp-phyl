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

using namespace bpp;

/******************************************************************************/

F84::F84(
	const NucleicAlphabet * alpha,
	double kappa,
	double piA,
	double piC,
	double piG,
	double piT):
	NucleotideSubstitutionModel(alpha, "F84.")
{
	Parameter kappaP("F84.kappa", kappa, &Parameter::R_PLUS);
	addParameter_(kappaP);
  _theta = piG + piC;
  _theta1 = piA / (1. - _theta);
  _theta2 = piG / _theta;
	Parameter thetaP("F84.theta" , _theta , &Parameter::PROP_CONSTRAINT_EX);
	addParameter_(thetaP);
	Parameter theta1P("F84.theta1", _theta1, &Parameter::PROP_CONSTRAINT_EX);
	addParameter_(theta1P);
	Parameter theta2P("F84.theta2", _theta2, &Parameter::PROP_CONSTRAINT_EX);
	addParameter_(theta2P);
  _p.resize(size_, size_);
	updateMatrices();
}

/******************************************************************************/

void F84::updateMatrices()
{
	_kappa = getParameterValue("kappa");
  _theta  = getParameterValue("theta");
	_theta1 = getParameterValue("theta1");
	_theta2 = getParameterValue("theta2");
  _piA = _theta1 * (1. - _theta);
  _piC = (1. - _theta2) * _theta;
  _piG = _theta2 * _theta;
  _piT = (1. - _theta1) * (1. - _theta);
	_piR = _piA + _piG;
	_piY = _piT + _piC;
	_r = 1. / (1 - _piA * _piA - _piC * _piC - _piG * _piG - _piT*_piT + 2. * _kappa * (_piC * _piT / _piY + _piA * _piG / _piR));
	_k1 = 1;
	_k2 = _kappa + 1;
	
  // Frequences:
	freq_[0] = _piA;
	freq_[1] = _piC;
	freq_[2] = _piG;
	freq_[3] = _piT;

  //Generator:
	generator_(0, 0) = -(                  _piC + (1 + _kappa / _piR) * _piG +       _piT);
	generator_(1, 1) = -(      _piA +                   _piG + (1 + _kappa / _piY) * _piT); 
	generator_(2, 2) = -((1 + _kappa / _piR) * _piA +       _piC             +       _piT);
	generator_(3, 3) = -(      _piA + (1 + _kappa / _piY) * _piC +       _piG            );

	generator_(1, 0) = _piA;
	generator_(3, 0) = _piA;
	generator_(0, 1) = _piC;
	generator_(2, 1) = _piC;
	generator_(1, 2) = _piG;
	generator_(3, 2) = _piG;
	generator_(0, 3) = _piT;
	generator_(2, 3) = _piT;
	
	generator_(2, 0) = (1 + _kappa / _piR) * _piA;
	generator_(3, 1) = (1 + _kappa / _piY) * _piC;
	generator_(0, 2) = (1 + _kappa / _piR) * _piG;
	generator_(1, 3) = (1 + _kappa / _piY) * _piT;
	
	// Normalization:
	MatrixTools::scale(generator_, _r);
	
	// Exchangeability:
	exchangeability_(0,0) = generator_(0,0) / _piA;
	exchangeability_(0,1) = generator_(0,1) / _piC; 
	exchangeability_(0,2) = generator_(0,2) / _piG; 
	exchangeability_(0,3) = generator_(0,3) / _piT;

	exchangeability_(1,0) = generator_(1,0) / _piA; 
	exchangeability_(1,1) = generator_(1,1) / _piC; 
	exchangeability_(1,2) = generator_(1,2) / _piG; 
	exchangeability_(1,3) = generator_(1,3) / _piT; 
	
	exchangeability_(2,0) = generator_(2,0) / _piA; 
	exchangeability_(2,1) = generator_(2,1) / _piC; 
	exchangeability_(2,2) = generator_(2,2) / _piG; 
	exchangeability_(2,3) = generator_(2,3) / _piT; 
	
	exchangeability_(3,0) = generator_(3,0) / _piA;
	exchangeability_(3,1) = generator_(3,1) / _piC; 
	exchangeability_(3,2) = generator_(3,2) / _piG; 
	exchangeability_(3,3) = generator_(3,3) / _piT;

	// Eigen values:
	eigenValues_[0] = 0;
	eigenValues_[1] = -_r * (1+_kappa);
	eigenValues_[2] = -_r * (1+_kappa); 
	eigenValues_[3] = -_r;
	
	// Eigen vectors:
	leftEigenVectors_(0,0) = _piA;
	leftEigenVectors_(0,1) = _piC;
	leftEigenVectors_(0,2) = _piG;
	leftEigenVectors_(0,3) = _piT;

	leftEigenVectors_(1,0) = 0.;
	leftEigenVectors_(1,1) = _piT / _piY;
	leftEigenVectors_(1,2) = 0.;
	leftEigenVectors_(1,3) = -_piT / _piY;

	leftEigenVectors_(2,0) = _piG / _piR;
	leftEigenVectors_(2,1) = 0.;
	leftEigenVectors_(2,2) = -_piG / _piR;
	leftEigenVectors_(2,3) = 0.;

	leftEigenVectors_(3,0) = _piA*_piY / _piR;
	leftEigenVectors_(3,1) = -_piC;
	leftEigenVectors_(3,2) = _piG*_piY / _piR;
	leftEigenVectors_(3,3) = -_piT;

	rightEigenVectors_(0,0) = 1.;
	rightEigenVectors_(0,1) = 0.;
	rightEigenVectors_(0,2) = 1.;
	rightEigenVectors_(0,3) = 1.;
	
	rightEigenVectors_(1,0) = 1.;
	rightEigenVectors_(1,1) = 1.;
	rightEigenVectors_(1,2) = 0.;;
	rightEigenVectors_(1,3) = -_piR / _piY;

	rightEigenVectors_(2,0) = 1.;
	rightEigenVectors_(2,1) = 0.;
	rightEigenVectors_(2,2) = -_piA / _piG;
	rightEigenVectors_(2,3) = 1.;

	rightEigenVectors_(3,0) = 1.;
	rightEigenVectors_(3,1) = -_piC / _piT;
	rightEigenVectors_(3,2) = 0.;
	rightEigenVectors_(3,3) = -_piR / _piY;
}
	
/******************************************************************************/

double F84::Pij_t(int i, int j, double d) const
{
	_l = _r * d;
	_exp1 = exp(-_k1*_l);
	_exp2 = exp(-_k2*_l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return _piA * (1. + (_piY/_piR) * _exp1) + (_piG/_piR) * _exp2; //A
				case 1 : return _piC * (1. -               _exp1);                       //C
				case 2 : return _piG * (1. + (_piY/_piR) * _exp1) - (_piG/_piR) * _exp2; //G
				case 3 : return _piT * (1. -               _exp1);                       //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return _piA * (1. -               _exp1);                       //A
				case 1 : return _piC * (1. + (_piR/_piY) * _exp1) + (_piT/_piY) * _exp2; //C
				case 2 : return _piG * (1. -               _exp1);                       //G
				case 3 : return _piT * (1. + (_piR/_piY) * _exp1) - (_piT/_piY) * _exp2; //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return _piA * (1. + (_piY/_piR) * _exp1) - (_piA/_piR) * _exp2; //A
				case 1 : return _piC * (1. -               _exp1);                       //C
				case 2 : return _piG * (1. + (_piY/_piR) * _exp1) + (_piA/_piR) * _exp2; //G
				case 3 : return _piT * (1. -               _exp1);                       //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return _piA * (1. -               _exp1);                       //A
				case 1 : return _piC * (1. + (_piR/_piY) * _exp1) - (_piC/_piY) * _exp2; //C
				case 2 : return _piG * (1. -               _exp1);                       //G
				case 3 : return _piT * (1. + (_piR/_piY) * _exp1) + (_piC/_piY) * _exp2; //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double F84::dPij_dt(int i, int j, double d) const
{
	_l = _r * d;
	_exp1 = exp(-_k1*_l);
	_exp2 = exp(-_k2*_l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return _r * (_piA * -(_piY/_piR) * _exp1 - (_piG/_piR) * _k2 * _exp2); //A
				case 1 : return _r * (_piC *                _exp1);                             //C
				case 2 : return _r * (_piG * -(_piY/_piR) * _exp1 + (_piG/_piR) * _k2 * _exp2); //G
				case 3 : return _r * (_piT *                _exp1);                             //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return _r * (_piA *                _exp1);                             //A
				case 1 : return _r * (_piC * -(_piR/_piY) * _exp1 - (_piT/_piY) * _k2 * _exp2); //C
				case 2 : return _r * (_piG *                _exp1);                             //G
				case 3 : return _r * (_piT * -(_piR/_piY) * _exp1 + (_piT/_piY) * _k2 * _exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return _r * (_piA * -(_piY/_piR) * _exp1 + (_piA/_piR) * _k2 * _exp2); //A
				case 1 : return _r * (_piC *                _exp1);                             //C
				case 2 : return _r * (_piG * -(_piY/_piR) * _exp1 - (_piA/_piR) * _k2 * _exp2); //G
				case 3 : return _r * (_piT *                _exp1);                             //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return _r * (_piA *                _exp1);                             //A
				case 1 : return _r * (_piC * -(_piR/_piY) * _exp1 + (_piC/_piY) * _k2 * _exp2); //C
				case 2 : return _r * (_piG *                _exp1);                             //G
				case 3 : return _r * (_piT * -(_piR/_piY) * _exp1 - (_piC/_piY) * _k2 * _exp2); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double F84::d2Pij_dt2(int i, int j, double d) const
{
	double r_2 = _r * _r;
	_l = _r * d;
	double k2_2 = _k2 * _k2;
	_exp1 = exp(-_k1*_l);
	_exp2 = exp(-_k2*_l);
	
	switch(i) {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r_2 * (_piA * (_piY/_piR) * _exp1 + (_piG/_piR) * k2_2 * _exp2); //A
				case 1 : return r_2 * (_piC *             - _exp1);                              //C
				case 2 : return r_2 * (_piG * (_piY/_piR) * _exp1 - (_piG/_piR) * k2_2 * _exp2); //G
				case 3 : return r_2 * (_piT *             - _exp1);                              //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return r_2 * (_piA *             - _exp1);                              //A
				case 1 : return r_2 * (_piC * (_piR/_piY) * _exp1 + (_piT/_piY) * k2_2 * _exp2); //C
				case 2 : return r_2 * (_piG *             - _exp1);                              //G
				case 3 : return r_2 * (_piT * (_piR/_piY) * _exp1 - (_piT/_piY) * k2_2 * _exp2); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r_2 * (_piA * (_piY/_piR) * _exp1 - (_piA/_piR) * k2_2 * _exp2); //A
				case 1 : return r_2 * (_piC *             - _exp1);                              //C
				case 2 : return r_2 * (_piG * (_piY/_piR) * _exp1 + (_piA/_piR) * k2_2 * _exp2); //G
				case 3 : return r_2 * (_piT *             - _exp1);                              //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r_2 * (_piA *             - _exp1);                              //A
				case 1 : return r_2 * (_piC * (_piR/_piY) * _exp1 - (_piC/_piY) * k2_2 * _exp2); //C
				case 2 : return r_2 * (_piG *             - _exp1);                              //G
				case 3 : return r_2 * (_piT * (_piR/_piY) * _exp1 + (_piC/_piY) * k2_2 * _exp2); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

const Matrix<double> & F84::getPij_t(double d) const
{
	_l = _r * d;
	_exp1 = exp(-_k1*_l);
	_exp2 = exp(-_k2*_l);

	//A
	_p(0, 0) = _piA * (1. + (_piY/_piR) * _exp1) + (_piG/_piR) * _exp2; //A
	_p(0, 1) = _piC * (1. -               _exp1);                       //C
	_p(0, 2) = _piG * (1. + (_piY/_piR) * _exp1) - (_piG/_piR) * _exp2; //G
	_p(0, 3) = _piT * (1. -               _exp1);                       //T, U

	//C
	_p(1, 0) = _piA * (1. -               _exp1);                       //A
	_p(1, 1) = _piC * (1. + (_piR/_piY) * _exp1) + (_piT/_piY) * _exp2; //C
	_p(1, 2) = _piG * (1. -               _exp1);                       //G
	_p(1, 3) = _piT * (1. + (_piR/_piY) * _exp1) - (_piT/_piY) * _exp2; //T, U

	//G
	_p(2, 0) = _piA * (1. + (_piY/_piR) * _exp1) - (_piA/_piR) * _exp2; //A
	_p(2, 1) = _piC * (1. -               _exp1);                       //C
	_p(2, 2) = _piG * (1. + (_piY/_piR) * _exp1) + (_piA/_piR) * _exp2; //G
	_p(2, 3) = _piT * (1. -               _exp1);                       //T, U

	//T, U
	_p(3, 0) = _piA * (1. -               _exp1);                       //A
	_p(3, 1) = _piC * (1. + (_piR/_piY) * _exp1) - (_piC/_piY) * _exp2; //C
	_p(3, 2) = _piG * (1. -               _exp1);                       //G
	_p(3, 3) = _piT * (1. + (_piR/_piY) * _exp1) + (_piC/_piY) * _exp2; //T, U

	return _p;
}

const Matrix<double> & F84::getdPij_dt(double d) const
{
	_l = _r * d;
	_exp1 = exp(-_k1*_l);
	_exp2 = exp(-_k2*_l);

	//A
	_p(0, 0) = _r * (_piA * -(_piY/_piR) * _exp1 - (_piG/_piR) * _k2 * _exp2); //A
	_p(0, 1) = _r * (_piC *                _exp1);                             //C
	_p(0, 2) = _r * (_piG * -(_piY/_piR) * _exp1 + (_piG/_piR) * _k2 * _exp2); //G
	_p(0, 3) = _r * (_piT *                _exp1);                             //T, U

	//C
	_p(1, 0) = _r * (_piA *                _exp1);                             //A
	_p(1, 1) = _r * (_piC * -(_piR/_piY) * _exp1 - (_piT/_piY) * _k2 * _exp2); //C
	_p(1, 2) = _r * (_piG *                _exp1);                             //G
	_p(1, 3) = _r * (_piT * -(_piR/_piY) * _exp1 + (_piT/_piY) * _k2 * _exp2); //T, U

	//G
	_p(2, 0) = _r * (_piA * -(_piY/_piR) * _exp1 + (_piA/_piR) * _k2 * _exp2); //A
	_p(2, 1) = _r * (_piC *                _exp1);                             //C
	_p(2, 2) = _r * (_piG * -(_piY/_piR) * _exp1 - (_piA/_piR) * _k2 * _exp2); //G
	_p(2, 3) = _r * (_piT *                _exp1);                             //T, U

	//T, U
	_p(3, 0) = _r * (_piA *                _exp1);                             //A
	_p(3, 1) = _r * (_piC * -(_piR/_piY) * _exp1 + (_piC/_piY) * _k2 * _exp2); //C
	_p(3, 2) = _r * (_piG *                _exp1);                             //G
	_p(3, 3) = _r * (_piT * -(_piR/_piY) * _exp1 - (_piC/_piY) * _k2 * _exp2); //T, U

	return _p;
}

const Matrix<double> & F84::getd2Pij_dt2(double d) const
{
	double r_2 = _r * _r;
	_l = _r * d;
	double k2_2 = _k2 * _k2;
	_exp1 = exp(-_k1*_l);
	_exp2 = exp(-_k2*_l);

	//A
	_p(0, 0) = r_2 * (_piA * (_piY/_piR) * _exp1 + (_piG/_piR) * k2_2 * _exp2); //A
	_p(0, 1) = r_2 * (_piC *             - _exp1);                              //C
	_p(0, 2) = r_2 * (_piG * (_piY/_piR) * _exp1 - (_piG/_piR) * k2_2 * _exp2); //G
	_p(0, 3) = r_2 * (_piT *             - _exp1);                              //T, U

	//C
	_p(1, 0) = r_2 * (_piA *             - _exp1);                              //A
	_p(1, 1) = r_2 * (_piC * (_piR/_piY) * _exp1 + (_piT/_piY) * k2_2 * _exp2); //C
	_p(1, 2) = r_2 * (_piG *             - _exp1);                              //G
	_p(1, 3) = r_2 * (_piT * (_piR/_piY) * _exp1 - (_piT/_piY) * k2_2 * _exp2); //T, U

	//G
	_p(2, 0) = r_2 * (_piA * (_piY/_piR) * _exp1 - (_piA/_piR) * k2_2 * _exp2); //A
	_p(2, 1) = r_2 * (_piC *             - _exp1);                              //C
	_p(2, 2) = r_2 * (_piG * (_piY/_piR) * _exp1 + (_piA/_piR) * k2_2 * _exp2); //G
	_p(2, 3) = r_2 * (_piT *             - _exp1);                              //T, U
 
	//T, U
	_p(3, 0) = r_2 * (_piA *             - _exp1);                              //A
	_p(3, 1) = r_2 * (_piC * (_piR/_piY) * _exp1 - (_piC/_piY) * k2_2 * _exp2); //C
	_p(3, 2) = r_2 * (_piG *             - _exp1);                              //G
	_p(3, 3) = r_2 * (_piT * (_piR/_piY) * _exp1 + (_piC/_piY) * k2_2 * _exp2); //T, U

	return _p;
}

/******************************************************************************/

void F84::setFreqFromData(const SequenceContainer & data, unsigned int pseudoCount)
{
	map<int, double> freqs = SequenceContainerTools::getFrequencies(data);
	double t = 0;
	for(unsigned int i = 0; i < size_; i++) t += freqs[i] + pseudoCount;
	_piA = (freqs[0] + pseudoCount) / t;
	_piC = (freqs[1] + pseudoCount) / t;
	_piG = (freqs[2] + pseudoCount) / t;
	_piT = (freqs[3] + pseudoCount) / t;
  vector<string> thetas(3);
  thetas[0] = getNamespace()+ "theta";
  thetas[1] = getNamespace() + "theta1";
  thetas[2] = getNamespace() + "theta2";
  ParameterList pl = getParameters().subList(thetas);
	pl[0]->setValue(_piC + _piG);
	pl[1]->setValue(_piA / (_piA + _piT));
	pl[2]->setValue(_piG / (_piC + _piG));
  setParametersValues(pl);
}

/******************************************************************************/

