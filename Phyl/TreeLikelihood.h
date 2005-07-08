//
// File: TreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 17 17:36:44 2003
//

/*
Copyright ou � ou Copr. CNRS, (16 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant � fournir des classes
pour l'analyse de donn�es phylog�n�tiques.

Ce logiciel est r�gi par la licence CeCILL soumise au droit fran�ais et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffus�e par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilit� au code source et des droits de copie,
de modification et de redistribution accord�s par cette licence, il n'est
offert aux utilisateurs qu'une garantie limit�e.  Pour les m�mes raisons,
seule une responsabilit� restreinte p�se sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les conc�dants successifs.

A cet �gard  l'attention de l'utilisateur est attir�e sur les risques
associ�s au chargement,  � l'utilisation,  � la modification et/ou au
d�veloppement et � la reproduction du logiciel par l'utilisateur �tant 
donn� sa sp�cificit� de logiciel libre, qui peut le rendre complexe � 
manipuler et qui le r�serve donc � des d�veloppeurs et des professionnels
avertis poss�dant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invit�s � charger  et  tester  l'ad�quation  du
logiciel � leurs besoins dans des conditions permettant d'assurer la
s�curit� de leurs syst�mes et ou de leurs donn�es et, plus g�n�ralement, 
� l'utiliser et l'exploiter dans les m�mes conditions de s�curit�. 

Le fait que vous puissiez acc�der � cet en-t�te signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accept� les
termes.
*/

/*
Copyright or � or Copr. CNRS, (November 16, 2004)

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

#ifndef _TREELIKELIHOOD_H_
#define _TREELIKELIHOOD_H_

#include "Tree.h"
#include "Node.h"

// From NumCalc:
#include <NumCalc/ParameterList.h>
#include <NumCalc/Parametrizable.h>
#include <NumCalc/Functions.h>
#include <NumCalc/VectorTools.h>

// From SeqLib:
#include <Seq/Alphabet.h>

/**
 * @brief The TreeLikelihood interface.
 
 * This interface defines the methods needed for computing the likelihood
 * of a phylogenetic tree, given a dataset.
 *
 * Likelihood computing isvery often memory and CPU expensive,
 * hence many algorithms try to store as information as possible to save computations.
 * This interface separates the computation itself (computeLikelihood() method) and the
 * result (getLikelihood() methods).
 */ 
class TreeLikelihood: public virtual DerivableSecondOrder
{
	public:
		virtual ~TreeLikelihood();
	
	public:
		
		/**
		 * @brief Get the likelihood for a site.
		 *
		 * @param site The site index to analyse.
		 * @return The likelihood for site <i>site</i>.
		 */
		virtual double getLikelihoodForASite(unsigned int site) const = 0;

		/**
		 * @brief Get the logarithm of the likelihood for a site.
		 *
		 * @param site The site index to analyse.
		 * @return The logarithm of the likelihood for site <i>site</i>.
		 */
		virtual double getLogLikelihoodForASite(unsigned int site) const = 0;

		/**
		 * @brief Get the likelihood for a site and for a state.
		 *
		 * @param site The site index to analyse.
		 * @param state The state to consider.
		 * @return The likelihood for site <i>site</i> and state <i>state</i>.
		 */
		virtual double getLikelihoodForASiteForAState(unsigned int site, int state) const = 0;

		/**
		 * @brief Get the logarithm of the likelihood for a site and for a state.
		 *
		 * @param site The site index to analyse.
		 * @param state The state to consider.
		 * @return The logarithm of the likelihood for site <i>site</i> and state <i>state</i>.
		 */
		virtual double getLogLikelihoodForASiteForAState(unsigned int site, int state) const = 0;

		/**
		 * @brief Get the likelihood for each site.
		 *
		 * @return A vector with all likelihoods for each site.
		 */
		virtual Vdouble getLikelihoodForEachSite() const = 0;

		/**
		 * @brief Get the logarithm of the likelihood for each site.
		 *
		 * @return A vector with all log likelihoods for each site.
		 */
		virtual Vdouble getLogLikelihoodForEachSite() const = 0;

		/**
		 * @brief Get the likelihood for each site and for each state.
		 *
		 * @return A 2d vector with all likelihoods for each site and for each state.
		 */
		virtual VVdouble getLikelihoodForEachSiteForEachState() const = 0;

		/**
		 * @brief Get the logarithm of the likelihood for each site and for each state.
		 *
		 * @return A 2d vector with all log likelihoods for each site and for each state.
		 */
		virtual VVdouble getLogLikelihoodForEachSiteForEachState() const = 0;
		
		/**
		 * @brief Get the likelihood for the whole dataset.
		 *
		 * @return The likelihood of the dataset.
		 */
		virtual double getLikelihood() const = 0;

		/**
		 * @brief Get the logarithm of the likelihood for the whole dataset.
		 *
		 * @return The logarithm of the likelihood of the dataset.
		 */
		virtual double getLogLikelihood() const = 0;
	
		/**
		 * @brief Get the tree (topology and branch lengths).
		 *
		 * @return The tree of this TreeLikelihood object.
	 	 */
		virtual Tree<Node> * getTree() const = 0;

		/**
		 * @brief Get the number of sites in the dataset.
		 *
		 * @return the number of sites in the dataset.
		 */
		virtual unsigned int getNumberOfSites() const = 0;

		/**
		 * @brief Get the number of states in the alphabet associated to the dataset.
		 *
		 * @return the number of states in the alphabet associated to the dataset.
		 */		
		virtual unsigned int getNumberOfStates() const = 0;
		
		/**
		 * @brief Get the alphabet associated to the dataset.
		 *
		 * @return the alphabet associated to the dataset.
		 */		
		virtual const Alphabet * getAlphabet() const = 0;
		
		/**
		 * @name Retrieve some particular parameters subsets.
		 *
		 * @{
		 */
		
		/**
		 * @brief Get the branch lengths parameters.
		 *
		 * @return A ParameterList with all branch lengths.
		 */
		virtual ParameterList getBranchLengthsParameters() const = 0;
		
		/**
		 * @brief Get the parameters assoicated to substitution model(s).
		 *
		 * @return A ParameterList.
		 */
		virtual ParameterList getSubstitutionModelParameters() const = 0;
		
		/** @} */

		/**
		 * @brief Tell if the derivatives must be computed.
		 *
		 * @param Yes or no.
		 */
		virtual void setComputeDerivatives(bool yn) = 0;

		/**
		 * @brief Tell if the derivatives must be computed.
		 *
		 * @return Yes or no.
		 */
		virtual bool computeDerivatives() const = 0;


};


#endif	//_TREELIKELIHOOD_H_
