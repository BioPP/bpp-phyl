//
// File: HKY85.cpp
// Created by: Julien Dutheil
// Created on: Thu Jan 22 16:17:39 2004
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

#include "HKY85.h"

// From SeqLib:
#include <Seq/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <cmath>

// From NumCalc:
#include <NumCalc/MatrixTools.h>

using namespace std;

/******************************************************************************/

HKY85::HKY85(
	const NucleicAlphabet * alpha,
	double kappa,
	double piA,
	double piC,
	double piG,
	double piT):
	NucleotideSubstitutionModel(alpha)
{
	_parameters.addParameter(Parameter("kappa", kappa, &Parameter::R_PLUS));
  _theta = piG + piC;
  _theta1 = piA / (1. - _theta);
  _theta2 = piG / _theta;
	_parameters.addParameter(Parameter("theta" , _theta , &Parameter::PROP_CONSTRAINT_EX));
	_parameters.addParameter(Parameter("theta1", _theta1, &Parameter::PROP_CONSTRAINT_EX));
	_parameters.addParameter(Parameter("theta2", _theta2, &Parameter::PROP_CONSTRAINT_EX));
  _p.resize(_size, _size);
	updateMatrices();
}

/******************************************************************************/

void HKY85::updateMatrices()
{
	_kappa  = _parameters.getParameter("kappa")->getValue();
	_theta  = _parameters.getParameter("theta")->getValue();
	_theta1 = _parameters.getParameter("theta1")->getValue();
	_theta2 = _parameters.getParameter("theta2")->getValue();
  _piA = _theta1 * (1. - _theta);
  _piC = (1. - _theta2) * _theta;
  _piG = _theta2 * _theta;
  _piT = (1. - _theta1) * (1. - _theta);
	_piR   = _piA + _piG;
	_piY   = _piT + _piC;
	_k1    = _kappa * _piY + _piR;
	_k2    = _kappa * _piR + _piY;

  _freq[0] = _piA;
  _freq[1] = _piC;
  _freq[2] = _piG;
  _freq[3] = _piT;
	
	_generator(0, 0) = -(                     _piC + _kappa*_piG +        _piT);
	_generator(1, 1) = -(       _piA +                      _piG + _kappa*_piT); 
	_generator(2, 2) = -(_kappa*_piA +        _piC               +       _piT);
	_generator(3, 3) = -(       _piA + _kappa*_piC +        _piG             );

	_generator(1, 0) = _piA;
	_generator(3, 0) = _piA;
	_generator(0, 1) = _piC;
	_generator(2, 1) = _piC;
	_generator(1, 2) = _piG;
	_generator(3, 2) = _piG;
	_generator(0, 3) = _piT;
	_generator(2, 3) = _piT;
	
	_generator(2, 0) = _kappa * _piA;
	_generator(3, 1) = _kappa * _piC;
	_generator(0, 2) = _kappa * _piG;
	_generator(1, 3) = _kappa * _piT;
	
	// Normalization:
	_r = 1. / (2. * (_piA * _piC + _piC * _piG + _piA * _piT + _piG * _piT + _kappa * (_piC * _piT + _piA * _piG)));
	MatrixTools::scale(_generator, _r);
	
	// Exchangeability:
	_exchangeability(0,0) = _generator(0,0) / _piA;
	_exchangeability(0,1) = _generator(0,1) / _piC; 
	_exchangeability(0,2) = _generator(0,2) / _piG; 
	_exchangeability(0,3) = _generator(0,3) / _piT;

	_exchangeability(1,0) = _generator(1,0) / _piA; 
	_exchangeability(1,1) = _generator(1,1) / _piC; 
	_exchangeability(1,2) = _generator(1,2) / _piG; 
	_exchangeability(1,3) = _generator(1,3) / _piT; 
	
	_exchangeability(2,0) = _generator(2,0) / _piA; 
	_exchangeability(2,1) = _generator(2,1) / _piC; 
	_exchangeability(2,2) = _generator(2,2) / _piG; 
	_exchangeability(2,3) = _generator(2,3) / _piT; 
	
	_exchangeability(3,0) = _generator(3,0) / _piA;
	_exchangeability(3,1) = _generator(3,1) / _piC; 
	_exchangeability(3,2) = _generator(3,2) / _piG; 
	_exchangeability(3,3) = _generator(3,3) / _piT;

	// Eigen values:
	_eigenValues[0] = 0;
	_eigenValues[1] = -_r * (_kappa * _piY + _piR);
	_eigenValues[2] = -_r * (_kappa * _piR + _piY); 
	_eigenValues[3] = -_r;
	
	// Eigen vectors:
	_leftEigenVectors(0,0) = _piA;
	_leftEigenVectors(0,1) = _piC;
	_leftEigenVectors(0,2) = _piG;
	_leftEigenVectors(0,3) = _piT;

	_leftEigenVectors(1,0) = 0.;
	_leftEigenVectors(1,1) = _piT / _piY;
	_leftEigenVectors(1,2) = 0.;
	_leftEigenVectors(1,3) = -_piT / _piY;

	_leftEigenVectors(2,0) = _piG / _piR;
	_leftEigenVectors(2,1) = 0.;
	_leftEigenVectors(2,2) = -_piG / _piR;
	_leftEigenVectors(2,3) = 0.;

	_leftEigenVectors(3,0) = _piA*_piY / _piR;
	_leftEigenVectors(3,1) = -_piC;
	_leftEigenVectors(3,2) = _piG*_piY / _piR;
	_leftEigenVectors(3,3) = -_piT;

	_rightEigenVectors(0,0) = 1.;
	_rightEigenVectors(0,1) = 0.;
	_rightEigenVectors(0,2) = 1.;
	_rightEigenVectors(0,3) = 1.;
	
	_rightEigenVectors(1,0) = 1.;
	_rightEigenVectors(1,1) = 1.;
	_rightEigenVectors(1,2) = 0.;;
	_rightEigenVectors(1,3) = -_piR / _piY;

	_rightEigenVectors(2,0) = 1.;
	_rightEigenVectors(2,1) = 0.;
	_rightEigenVectors(2,2) = -_piA / _piG;
	_rightEigenVectors(2,3) = 1.;

	_rightEigenVectors(3,0) = 1.;
	_rightEigenVectors(3,1) = -_piC / _piT;
	_rightEigenVectors(3,2) = 0.;
	_rightEigenVectors(3,3) = -_piR / _piY;
}
	
/******************************************************************************/

double HKY85::Pij_t(int i, int j, double d) const
{
	_l     = _r * d;
	_exp1  = exp(-_l);
	_exp22 = exp(-_k2 * _l);
	_exp21 = exp(-_k1 * _l);
	
	switch(i)
  {
		//A
		case 0 : {
			switch(j) {
				case 0 : return _piA * (1. + (_piY/_piR) * _exp1) + (_piG/_piR) * _exp22; //A
				case 1 : return _piC * (1. -               _exp1);                        //C
				case 2 : return _piG * (1. + (_piY/_piR) * _exp1) - (_piG/_piR) * _exp22; //G
				case 3 : return _piT * (1. -               _exp1);                        //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return _piA * (1. -               _exp1);                        //A
				case 1 : return _piC * (1. + (_piR/_piY) * _exp1) + (_piT/_piY) * _exp21; //C
				case 2 : return _piG * (1. -               _exp1);                        //G
				case 3 : return _piT * (1. + (_piR/_piY) * _exp1) - (_piT/_piY) * _exp21; //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return _piA * (1. + (_piY/_piR) * _exp1) - (_piA/_piR) * _exp22; //A
				case 1 : return _piC * (1. -               _exp1);                        //C
				case 2 : return _piG * (1. + (_piY/_piR) * _exp1) + (_piA/_piR) * _exp22; //G
				case 3 : return _piT * (1. -               _exp1);                        //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return _piA * (1. -               _exp1);                        //A
				case 1 : return _piC * (1. + (_piR/_piY) * _exp1) - (_piC/_piY) * _exp21; //C
				case 2 : return _piG * (1. -               _exp1);                        //G
				case 3 : return _piT * (1. + (_piR/_piY) * _exp1) + (_piC/_piY) * _exp21; //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double HKY85::dPij_dt(int i, int j, double d) const
{
	_l     = _r * d;
	_exp1  = exp(-_l);
	_exp22 = exp(-_k2 * _l);
	_exp21 = exp(-_k1 * _l);
	
	switch(i)
  {
		//A
		case 0 : {
			switch(j) {
				case 0 : return _r * (_piA * -(_piY/_piR) * _exp1 - (_piG/_piR) * _k2 * _exp22); //A
				case 1 : return _r * (_piC *                _exp1);                              //C
				case 2 : return _r * (_piG * -(_piY/_piR) * _exp1 + (_piG/_piR) * _k2 * _exp22); //G
				case 3 : return _r * (_piT *                _exp1);                              //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return _r * (_piA *                _exp1);                              //A
				case 1 : return _r * (_piC * -(_piR/_piY) * _exp1 - (_piT/_piY) * _k1 * _exp21); //C
				case 2 : return _r * (_piG *                _exp1);                              //G
				case 3 : return _r * (_piT * -(_piR/_piY) * _exp1 + (_piT/_piY) * _k1 * _exp21); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return _r * (_piA * -(_piY/_piR) * _exp1 + (_piA/_piR) * _k2 * _exp22); //A
				case 1 : return _r * (_piC *                _exp1);                              //C
				case 2 : return _r * (_piG * -(_piY/_piR) * _exp1 - (_piA/_piR) * _k2 * _exp22); //G
				case 3 : return _r * (_piT *                _exp1);                              //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return _r * (_piA *                _exp1);                              //A
				case 1 : return _r * (_piC * -(_piR/_piY) * _exp1 + (_piC/_piY) * _k1 * _exp21); //C
				case 2 : return _r * (_piG *                _exp1);                              //G
				case 3 : return _r * (_piT * -(_piR/_piY) * _exp1 - (_piC/_piY) * _k1 * _exp21); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

double HKY85::d2Pij_dt2(int i, int j, double d) const
{
	double r_2 = _r * _r;
	_l = _r * d;
	double k1_2 = _k1 * _k1;
	double k2_2 = _k2 * _k2;
	_exp1 = exp(-_l);
	_exp22 = exp(-_k2 * _l);
	_exp21 = exp(-_k1 * _l);
	
	switch(i)
  {
		//A
		case 0 : {
			switch(j) {
				case 0 : return r_2 * (_piA * (_piY/_piR) * _exp1 + (_piG/_piR) * k2_2 * _exp22); //A
				case 1 : return r_2 * (_piC *             - _exp1);                               //C
				case 2 : return r_2 * (_piG * (_piY/_piR) * _exp1 - (_piG/_piR) * k2_2 * _exp22); //G
				case 3 : return r_2 * (_piT *             - _exp1);                               //T, U
			}
		} 
		//C
		case 1 : {
			switch(j) {
				case 0 : return r_2 * (_piA *             - _exp1);                               //A
				case 1 : return r_2 * (_piC * (_piR/_piY) * _exp1 + (_piT/_piY) * k1_2 * _exp21); //C
				case 2 : return r_2 * (_piG *             - _exp1);                               //G
				case 3 : return r_2 * (_piT * (_piR/_piY) * _exp1 - (_piT/_piY) * k1_2 * _exp21); //T, U
			}
		}
		//G
		case 2 : {
			switch(j) {
				case 0 : return r_2 * (_piA * (_piY/_piR) * _exp1 - (_piA/_piR) * k2_2 * _exp22); //A
				case 1 : return r_2 * (_piC *             - _exp1);                               //C
				case 2 : return r_2 * (_piG * (_piY/_piR) * _exp1 + (_piA/_piR) * k2_2 * _exp22); //G
				case 3 : return r_2 * (_piT *             - _exp1);                               //T, U
			}
		}
		//T, U
		case 3 : {
			switch(j) {
				case 0 : return r_2 * (_piA *             - _exp1);                              //A
				case 1 : return r_2 * (_piC * (_piR/_piY) * _exp1 - (_piC/_piY) * k1_2 * _exp21); //C
				case 2 : return r_2 * (_piG *             - _exp1);                              //G
				case 3 : return r_2 * (_piT * (_piR/_piY) * _exp1 + (_piC/_piY) * k1_2 * _exp21); //T, U
			}
		}
	}
	return 0;
}

/******************************************************************************/

RowMatrix<double> HKY85::getPij_t(double d) const
{
  _l = _r * d;
	_exp1 = exp(-_l);
	_exp22 = exp(-_k2 * _l);
	_exp21 = exp(-_k1 * _l);

	//A
	_p(0, 0) = _piA * (1. + (_piY/_piR) * _exp1) + (_piG/_piR) * _exp22; //A
	_p(0, 1) = _piC * (1. -               _exp1);                        //C
	_p(0, 2) = _piG * (1. + (_piY/_piR) * _exp1) - (_piG/_piR) * _exp22; //G
	_p(0, 3) = _piT * (1. -               _exp1);                        //T, U

	//C
	_p(1, 0) = _piA * (1. -               _exp1);                        //A
	_p(1, 1) = _piC * (1. + (_piR/_piY) * _exp1) + (_piT/_piY) * _exp21; //C
	_p(1, 2) = _piG * (1. -               _exp1);                        //G
	_p(1, 3) = _piT * (1. + (_piR/_piY) * _exp1) - (_piT/_piY) * _exp21; //T, U

	//G
	_p(2, 0) = _piA * (1. + (_piY/_piR) * _exp1) - (_piA/_piR) * _exp22; //A
	_p(2, 1) = _piC * (1. -               _exp1);                        //C
	_p(2, 2) = _piG * (1. + (_piY/_piR) * _exp1) + (_piA/_piR) * _exp22; //G
	_p(2, 3) = _piT * (1. -               _exp1);                        //T, U

	//T, U
	_p(3, 0) = _piA * (1. -               _exp1);                        //A
	_p(3, 1) = _piC * (1. + (_piR/_piY) * _exp1) - (_piC/_piY) * _exp21; //C
	_p(3, 2) = _piG * (1. -               _exp1);                        //G
	_p(3, 3) = _piT * (1. + (_piR/_piY) * _exp1) + (_piC/_piY) * _exp21; //T, U

	return _p;
}

RowMatrix<double> HKY85::getdPij_dt(double d) const
{
	_l = _r * d;
	_exp1 = exp(-_l);
	_exp22 = exp(-_k2 * _l);
	_exp21 = exp(-_k1 * _l);

	//A
	_p(0, 0) = _r * (_piA * -(_piY/_piR) * _exp1 - (_piG/_piR) * _k2 * _exp22); //A
	_p(0, 1) = _r * (_piC *                _exp1);                              //C
	_p(0, 2) = _r * (_piG * -(_piY/_piR) * _exp1 + (_piG/_piR) * _k2 * _exp22); //G
	_p(0, 3) = _r * (_piT *                _exp1);                              //T, U

	//C
	_p(1, 0) = _r * (_piA *                _exp1);                              //A
	_p(1, 1) = _r * (_piC * -(_piR/_piY) * _exp1 - (_piT/_piY) * _k1 * _exp21); //C
	_p(1, 2) = _r * (_piG *                _exp1);                              //G
	_p(1, 3) = _r * (_piT * -(_piR/_piY) * _exp1 + (_piT/_piY) * _k1 * _exp21); //T, U

	//G
	_p(2, 0) = _r * (_piA * -(_piY/_piR) * _exp1 + (_piA/_piR) * _k2 * _exp22); //A
	_p(2, 1) = _r * (_piC *                _exp1);                              //C
	_p(2, 2) = _r * (_piG * -(_piY/_piR) * _exp1 - (_piA/_piR) * _k2 * _exp22); //G
	_p(2, 3) = _r * (_piT *                _exp1);                              //T, U

	//T, U
	_p(3, 0) = _r * (_piA *                _exp1);                              //A
	_p(3, 1) = _r * (_piC * -(_piR/_piY) * _exp1 + (_piC/_piY) * _k1 * _exp21); //C
	_p(3, 2) = _r * (_piG *                _exp1);                              //G
	_p(3, 3) = _r * (_piT * -(_piR/_piY) * _exp1 - (_piC/_piY) * _k1 * _exp21); //T, U

	return _p;
}

RowMatrix<double> HKY85::getd2Pij_dt2(double d) const
{
	double r_2 = _r * _r;
	_l = _r * d;
	double k1_2 = _k1 * _k1;
	double k2_2 = _k2 * _k2;
	_exp1 = exp(-_l);
	_exp22 = exp(-_k2 * _l);
	_exp21 = exp(-_k1 * _l);

	//A
	_p(0, 0) = r_2 * (_piA * (_piY/_piR) * _exp1 + (_piG/_piR) * k2_2 * _exp22); //A
	_p(0, 1) = r_2 * (_piC *             - _exp1);                               //C
	_p(0, 2) = r_2 * (_piG * (_piY/_piR) * _exp1 - (_piG/_piR) * k2_2 * _exp22); //G
	_p(0, 3) = r_2 * (_piT *             - _exp1);                               //T, U

	//C
	_p(1, 0) = r_2 * (_piA *             - _exp1);                               //A
	_p(1, 1) = r_2 * (_piC * (_piR/_piY) * _exp1 + (_piT/_piY) * k1_2 * _exp21); //C
	_p(1, 2) = r_2 * (_piG *             - _exp1);                               //G
	_p(1, 3) = r_2 * (_piT * (_piR/_piY) * _exp1 - (_piT/_piY) * k1_2 * _exp21); //T, U

	//G
	_p(2, 0) = r_2 * (_piA * (_piY/_piR) * _exp1 - (_piA/_piR) * k2_2 * _exp22); //A
	_p(2, 1) = r_2 * (_piC *             - _exp1);                               //C
	_p(2, 2) = r_2 * (_piG * (_piY/_piR) * _exp1 + (_piA/_piR) * k2_2 * _exp22); //G
	_p(2, 3) = r_2 * (_piT *             - _exp1);                               //T, U

	//T, U
	_p(3, 0) = r_2 * (_piA *             - _exp1);                               //A
	_p(3, 1) = r_2 * (_piC * (_piR/_piY) * _exp1 - (_piC/_piY) * k1_2 * _exp21); //C
	_p(3, 2) = r_2 * (_piG *             - _exp1);                               //G
	_p(3, 3) = r_2 * (_piT * (_piR/_piY) * _exp1 + (_piC/_piY) * k1_2 * _exp21); //T, U

	return _p;
}

/******************************************************************************/

void HKY85::setFreqFromData(const SequenceContainer & data)
{
	map<int, double> freqs = SequenceContainerTools::getFrequencies(data);
	double t = 0;
	for(unsigned int i = 0; i < _size; i++) t += freqs[i];
	_piA = freqs[0] / t;
	_piC = freqs[1] / t;
	_piG = freqs[2] / t;
	_piT = freqs[3] / t;
	_parameters.getParameter("theta")->setValue(_piC + _piG);
	_parameters.getParameter("theta1")->setValue(_piA / (_piA + _piT));
	_parameters.getParameter("theta2")->setValue(_piG / (_piC + _piG));
  updateMatrices();
}

/******************************************************************************/

