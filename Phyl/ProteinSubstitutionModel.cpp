//
// File: ProteinSubstitutionModel.cpp
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Wed Jan 21 13:59:18 2004
//

#include "ProteinSubstitutionModel.h"

// From the MTL:
#include <mtl/mtl.h>
#include <mtl/mtl2lapack.h>
#include <mtl/dense1D.h>
#include <mtl/utils.h>
using namespace mtl2lapack;

// From NumCalc:
#include <NumCalc/MatrixTools.h>
using namespace MatrixTools;

ProteinSubstitutionModel::ProteinSubstitutionModel(const Alphabet * alpha) :
	AbstractSubstitutionModel(alpha) {}

ProteinSubstitutionModel::~ProteinSubstitutionModel() {}

void ProteinSubstitutionModel::updateMatrices()
{
	// Now computes eigen values and vectors:
	Matrix Pi = diag<Matrix, double>(_freq);
	_generator = _exchangeability * Pi;
	// Normalization:
	double scale = getScale();
	mtl::scale(_generator, 1/scale);
	lapack_matrix<double>::type D(20, 20), R(20, 20), L(20, 20);
	copy(_generator, D);
	mtl::dense1D< complex<double> > wr(20);
	int info;
	info = geev(GEEV_CALC_RIGHT, D, wr, L, R);
	if(info > 0) throw Exception("ERROR!!! Failed to compute eigen values (convergence not reached).");
	//Check eigen values:
	for(unsigned int i = 0; i < 20; i++) {
		if(wr[i].imag() != 0) throw Exception("ERROR!!! Non real eigen value.");
		_eigenValues[i] = wr[i].real();
	}
	//copy(trans(L), _leftEigenVectors); //Not exactly equal to R^-1 !!!
	copy(R, _rightEigenVectors);
	
	//Compute left eigen vectors:
	copy(R, D);
	mtl::dense1D< int > ipivot(20);
	L = getId< lapack_matrix<double>::type >(20);
	info = gesv(D, ipivot, L);
	copy(L, _leftEigenVectors);
}


