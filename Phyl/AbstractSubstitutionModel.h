//
// File: AbstractSubstitutionModel.h
// Created by: Julien Dutheil
// Created on: Tue May 27 10:31:49 2003
//

/*
Copyright ou © ou Copr. Julien Dutheil, (16 Novembre 2004) 

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
Copyright or © or Copr. Julien Dutheil, (November 16, 2004)

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

#ifndef _ABSTRACTSUBSTITUTIONMODEL_H_
#define _ABSTRACTSUBSTITUTIONMODEL_H_

#include "SubstitutionModel.h"

/**
 * @brief Low level implementation of the SubstitutionModel interface.
 */
class AbstractSubstitutionModel : public virtual SubstitutionModel, public virtual AbstractParametrizable
{
	protected:

		/**
		 * @brief The alphabet this model is to use with.
		 */
		const Alphabet * alphabet;

		/**
		 * @brief The size of the alphabet.
		 */
		unsigned int _size;

		/**
		 * @brief The generator matrix \f$Q\f$ of the model.
		 */
		Mat _generator;

		/**
		 * @brief The exchangeability matrix \f$S\f$ of the model.
		 */
		Mat _exchangeability;

		/**
		 * @brief The \f$U^{-1}\f$ matrix made of horizontal left eigen vectors.
		 */
		Mat _leftEigenVectors;

		/**
		 * @brief The \f$U\f$ matrix made of vertical right eigen vectors.
		 */
		Mat _rightEigenVectors;

		/**
		 * @brief The vector of eigen values.
		 */
		Vec _eigenValues;

		/**
		 * @brief The vector of equilibrium frequencies.
		 */
		Vec _freq;

	public:
		AbstractSubstitutionModel(const Alphabet * alpha);
	
		virtual ~AbstractSubstitutionModel() {}
	
	public:
		const Alphabet * getAlphabet() const;

		Vec getFrequencies() const;
		Mat getExchangeabilityMatrix() const;
		Mat getGenerator() const;
		Mat getPij_t(double t) const;
		Mat getdPij_dt(double t) const;
		Mat getd2Pij_dt2(double t) const;
		Vec eigenValues() const;
		Mat horizontalLeftEigenVectors() const;
		Mat verticalRightEigenVectors() const;
		double freq(int i) const;
		double Qij(int i, int j) const;
		double Pij_t    (int i, int j, double t) const;
		double dPij_dt  (int i, int j, double t) const;
		double d2Pij_dt2(int i, int j, double t) const;

		double getInitValue(int i, int state) const throw (BadIntException);
		void setFreqFromData(const SequenceContainer & data);

		void fireParameterChanged(const ParameterList & parameters) {}
		
	protected:
		/**
		 * @brief Compute and diagonalize the \f$Q\f$ matrix and fill the _eigenValues,
		 * _leftEigenVectors and _rightEigenVectors fields.
		 *
		 * This routine uses the MTL interface to Lapack to compute eigen values
		 * and vectors. However, the inverse of the right eigen vectors is used
		 * as left eigen vectors, since the left eigen vectors sent by the lapack
		 * geev routine is not exactly the inverse of the right eigen vectors.
		 * (Dunno why though...)
		 */
		virtual void updateMatrices() = 0;

		/**
		 * @brief Get the scalar product of diagonal elements of the generator
		 * and the frequencies vector.
		 * If the generator is normalized, then scale=1. Otherwise each element
		 * must be multiplied by 1/scale.
		 *
		 * @return Minus the scalar product of diagonal elements and the frequencies vector.
		 */
		double getScale() const;
};


#endif	//_ABSTRACTSUBSTITUTIONMODEL_H_
