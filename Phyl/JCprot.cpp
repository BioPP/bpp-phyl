//
// File: JCprot.cpp
// Created by: Julien Dutheil
// Created on: Tue May 27 16:04:36 2003
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

#include "JCprot.h"

using namespace bpp;

#include <cmath>

using namespace std;

/******************************************************************************/

JCprot::JCprot(const ProteicAlphabet * alpha) :
  ProteinSubstitutionModel(alpha, "JC69.")
{
	_p.resize(size_, size_);
	updateMatrices();
}

/******************************************************************************/
	
void JCprot::updateMatrices()
{
	// Frequencies:
	for(unsigned int i = 0; i < 20; i++) freq_[i] = 1. / 20.;

	// Generator:
	for(unsigned int i = 0; i < 20; i++)
  {
		for(unsigned int j = 0; j < 20; j++)
    {
			generator_(i, j) = (i == j) ? -1. : 1./19.;
			exchangeability_(i, j) = generator_(i, j) * 20.;
		}
	}
	
	// Eigen values:
	eigenValues_[0] = 0;
	for(unsigned int i = 1; i < 20; i++) eigenValues_[i] = -20. / 19.;
	
	// Eigen vectors:
	for(unsigned int i = 0; i < 20; i++) leftEigenVectors_(0,i) = 1./20.;
	for(unsigned int i = 1; i < 20; i++) 
		for(unsigned int j = 0; j < 20; j++)
			leftEigenVectors_(i,j) = -1./20.;
	for(unsigned int i = 0; i < 19; i++) leftEigenVectors_(19-i,i) = 19./20.;

	for(unsigned int i = 0; i < 20; i++) rightEigenVectors_(i,0) = 1.;
	for(unsigned int i = 1; i < 20; i++) rightEigenVectors_(19,i) = -1.;
	for(unsigned int i = 0; i < 19; i++) 
		for(unsigned int j = 1; j < 20; j++)
			rightEigenVectors_(i,j) = 0.;
	for(unsigned int i = 1; i < 20; i++) rightEigenVectors_(19-i,i) = 1.;

}
	
/******************************************************************************/

double JCprot::Pij_t(int i, int j, double d) const
{
	if(i == j) return 1./20. + 19./20. * exp(- 20./19. * d);
	else       return 1./20. -  1./20. * exp(- 20./19. * d);
}

/******************************************************************************/

double JCprot::dPij_dt(int i, int j, double d) const
{
	if(i == j) return -        exp(- 20./19. * d);
	else       return 1./19. * exp(- 20./19. * d);
}

/******************************************************************************/

double JCprot::d2Pij_dt2(int i, int j, double d) const
{
	if(i == j) return   20./19.  * exp(- 20./19. * d);
	else       return - 20./361. * exp(- 20./19. * d);
}

/******************************************************************************/

const Matrix<double> & JCprot::getPij_t(double d) const
{
  _exp = exp(- 20./19. * d);
	for(unsigned int i = 0; i < size_; i++)
  {
		for(unsigned int j = 0; j < size_; j++)
    {
			_p(i,j) = (i==j) ? 1./20. + 19./20. * _exp : 1./20. - 1./20. * _exp;
		}
	}
	return _p;
}

const Matrix<double> & JCprot::getdPij_dt(double d) const
{
  _exp = exp(- 20./19. * d);
	for(unsigned int i = 0; i < size_; i++)
  {
		for(unsigned int j = 0; j < size_; j++)
    {
			_p(i,j) = (i==j) ? - _exp : 1./19. * _exp;
		}
	}
	return _p;
}

const Matrix<double> & JCprot::getd2Pij_dt2(double d) const
{
  _exp = exp(- 20./19. * d);
	for(unsigned int i = 0; i < size_; i++)
  {
		for(unsigned int j = 0; j < size_; j++)
    {
			_p(i,j) = (i==j) ? 20./19. * _exp : - 20./361. * _exp;
		}
	}
	return _p;
}

/******************************************************************************/

