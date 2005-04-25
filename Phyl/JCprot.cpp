//
// File: JCprot.cpp
// Created by:  <Julien.Dutheil@univ-montp2.fr>
// Created on: Tue May 27 16:04:36 2003
//

#include "JCprot.h"

#include <cmath>

/******************************************************************************/

JCprot::JCprot(const ProteicAlphabet * alpha): ProteinSubstitutionModel(alpha)
{
	_parameters = ParameterList(); //no parameters for this model.	
	
	// Frequencies:
	for(unsigned int i = 0; i < 20; i++) _freq[i] = 1. / 20.;	
	updateMatrices();
}

JCprot::~JCprot() {}

/******************************************************************************/
	
void JCprot::updateMatrices()
{
	// Generator:
	for(unsigned int i = 0; i < 20; i++) {
		for(unsigned int j = 0; j < 20; j++) {
			_generator(i, j) = (i == j) ? -1. : 1./19.;
		}
	}
	
	// Eigen values:
	_eigenValues[0] = 0;
	for(unsigned int i = 1; i < 20; i++) _eigenValues[i] = -20. / 19.;
	
	// Eigen vectors:
	// todo!
}
	
/******************************************************************************/

double JCprot::Pij_t(int i, int j, double d) const {
	if(i == j) return 1./20. + 19./20. * exp(- 20./19. * d);
	else       return 1./20. -  1./20. * exp(- 20./19. * d);
}

/******************************************************************************/

double JCprot::dPij_dt(int i, int j, double d) const {
	if(i == j) return -        exp(- 20./19. * d);
	else       return 1./19. * exp(- 20./19. * d);
}

/******************************************************************************/

double JCprot::d2Pij_dt2(int i, int j, double d) const {
	if(i == j) return   20./19.  * exp(- 20./19. * d);
	else       return - 20./361. * exp(- 20./19. * d);
}

/******************************************************************************/

Mat JCprot::getPij_t(double d) const {
	Mat p(_size, _size);
	for(unsigned int i = 0; i < _size; i++) {
		for(unsigned int j = 0; j < _size; j++) {
			p(i,j) = (i==j) ? 1./20. + 19./20. * exp(- 20./19. * d) : 1./20. - 1./20. * exp(- 20./19. * d);
		}
	}
	return p;
}

Mat JCprot::getdPij_dt(double d) const {
	Mat p(_size, _size);
	for(unsigned int i = 0; i < _size; i++) {
		for(unsigned int j = 0; j < _size; j++) {
			p(i,j) = (i==j) ? - exp(- 20./19. * d) : 1./19. * exp(- 20./19. * d);
		}
	}
	return p;
}

Mat JCprot::getd2Pij_dt2(double d) const {
	Mat p(_size, _size);
	for(unsigned int i = 0; i < _size; i++) {
		for(unsigned int j = 0; j < _size; j++) {
			p(i,j) = (i==j) ? 20./19. * exp(- 20./19. * d) : - 20./361. * exp(- 420./19. * d);
		}
	}
	return p;
}

/******************************************************************************/

string JCprot::getName() const {
	return string("Jukes and Cantor (1969) for proteins");
}

/******************************************************************************/
