//
// File: AbstractSubstitutionModel.cpp
// Created by: Julien Dutheil
// Created on: Tue May 27 10:31:49 2003
//

/*
Copyright ou © ou Copr. CNRS, (16 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant à fournir des classes
pour l'analyse de données phylogénétiques.

Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.
*/

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

Julien.Dutheil@univ-montp2.fr

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

#include "AbstractSubstitutionModel.h"

// From Utils:
#include <Utils/TextTools.h>

// From NumCalc:
#include <NumCalc/MatrixTools.h>
using namespace MatrixOperators;
#include <NumCalc/VectorTools.h>
using namespace VectorFunctions;
using namespace VectorOperators;

// From SeqLib:
#include <Seq/SequenceContainerTools.h>

AbstractSubstitutionModel::AbstractSubstitutionModel(const Alphabet * alpha): alphabet(alpha)
{
	_size = alpha -> getSize();
	_generator.resize(_size, _size);
	_exchangeability.resize(_size, _size);
	_freq.resize(_size);
	_eigenValues.resize(_size);
	_leftEigenVectors.resize(_size, _size);
	_rightEigenVectors.resize(_size, _size);
}

/******************************************************************************/
	
const Alphabet * AbstractSubstitutionModel::getAlphabet() const { return alphabet; }

/******************************************************************************/
	
Vec AbstractSubstitutionModel::getFrequencies() const { return _freq; }

Mat AbstractSubstitutionModel::getExchangeabilityMatrix() const { return _exchangeability; }

Mat AbstractSubstitutionModel::getGenerator() const { return _generator; }

Vec AbstractSubstitutionModel::eigenValues() const { return _eigenValues; }

Mat AbstractSubstitutionModel::horizontalLeftEigenVectors() const { return _leftEigenVectors; }

Mat AbstractSubstitutionModel::verticalRightEigenVectors() const { return _rightEigenVectors; }

Mat AbstractSubstitutionModel::getPij_t(double t) const
{
	return MatrixTools::mult(_rightEigenVectors, exp(_eigenValues*t), _leftEigenVectors);
}

Mat AbstractSubstitutionModel::getdPij_dt(double t) const
{
	return MatrixTools::mult(_rightEigenVectors, _eigenValues * exp(_eigenValues*t), _leftEigenVectors);
}

Mat AbstractSubstitutionModel::getd2Pij_dt2(double t) const
{
	return MatrixTools::mult(_rightEigenVectors, sqr(_eigenValues) * exp(_eigenValues*t), _leftEigenVectors);
}

double AbstractSubstitutionModel::freq(int i) const { return _freq[i]; }

double AbstractSubstitutionModel::Qij(int i, int j) const { return _generator(i, j); }



double AbstractSubstitutionModel::Pij_t(int i, int j, double t) const {	return getPij_t(t)(i,j); }

double AbstractSubstitutionModel::dPij_dt(int i, int j, double t) const { return getdPij_dt(t)(i,j); }

double AbstractSubstitutionModel::d2Pij_dt2(int i, int j, double t) const { return getd2Pij_dt2(t)(i,j); }	

/******************************************************************************/

double AbstractSubstitutionModel::getInitValue(int i, int state) const throw (BadIntException) {
	if(i     < 0 || i     > (int)alphabet -> getSize()) throw BadIntException(i    , "AbstractSubstitutionModel::getInitValue");
	//Old method: do not care about generic characters:
	//if(state < 0 || state > (int)alphabet -> getSize()) throw BadIntException(state, "AbstractSubstitutionModel::getInitValue. Character " + alphabet -> intToChar(state) + " is not allowed in model.");
	//return i == state ? 1 : 0;
	if(state < 0 || !alphabet -> isIntInAlphabet(state)) throw BadIntException(state, "AbstractSubstitutionModel::getInitValue. Character " + alphabet -> intToChar(state) + " is not allowed in model.");
	vector<int> states = alphabet -> getAlias(state);
	for(unsigned int j = 0; j < states.size(); j++) if(i == states[j]) return 1;
	return 0;
}

/******************************************************************************/

void AbstractSubstitutionModel::setFreqFromData(const SequenceContainer & data) {
	map<int, double> freqs = SequenceContainerTools::getFrequencies(data);
	double t = 0;
	for(unsigned int i = 0; i < _size; i++) t += freqs[i];
	for(unsigned int i = 0; i < _size; i++) _freq[i] = freqs[i] / t;
	//Re-compute generator and eigen values:
	updateMatrices();
}

/******************************************************************************/

double AbstractSubstitutionModel::getScale() const {
	return -scalar(MatrixTools::diag<Mat, double>(_generator), _freq);
}

/******************************************************************************/

