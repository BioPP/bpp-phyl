//
// File: AbstractHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thr Dec 23 12:03 2004
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

#ifndef _ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_
#define _ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_

#include "AbstractDiscreteRatesAcrossSitesTreeLikelihood.h"
#include "SubstitutionModel.h"

class AbstractHomogeneousTreeLikelihood: public virtual AbstractDiscreteRatesAcrossSitesTreeLikelihood
{
	protected:
		SubstitutionModel * _model;
		ParameterList _brLenParameters;
		
		mutable map<const Node *, VVVdouble> _pxy;

		mutable map<const Node *, VVVdouble> _dpxy;

		mutable map<const Node *, VVVdouble> _d2pxy;
		
		/**
		 * @brief As previous, but for the global container.
		 *
		 * The size of this vector is equal to the number of sites in the container,
		 * each element corresponds to a site in the container and points to the
		 * corresponding column in the likelihood array of the root node.
		 * If the container contains no repeated site, there will be a strict
		 * equivalence between each site and the likelihood array of the root node.
		 * However, if this is not the case, some pointers may point toward the same
		 * element in the likelihood array.
		 */
		vector<unsigned int> _rootPatternLinks;

		/**
		 * @brief The frequency of each site.
		 */
		vector<unsigned int> _rootWeights;
		
		/**
		 * @brief Pointer toward all nodes in the tree.
		 *
		 * This is used for parameters estimation only, not for likelihood computation.
		 * The order of the nodes in the vector if the order of the named branches and
		 * is initalized once for all in the constructor. It then depends of the tree
		 * topology. This may lead to some problems when we'll act on tree topology...
		 */
		 vector<Node *> _nodes;

		//some values we'll need:
		unsigned int _nbSites,         //the number of sites in the container
		             _nbClasses,       //the number of rate classes
		             _nbStates,        //the number of states in the alphabet
		             _nbNodes;         //the number of nodes in the tree


	public:
		AbstractHomogeneousTreeLikelihood(
			Tree<Node> & tree,
			const SiteContainer & data,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
			bool verbose = true
			)	throw (Exception);

		virtual ~AbstractHomogeneousTreeLikelihood();
		
	public:
		
		/**
		 * @name The TreeLikelihood interface.
		 *
		 * Other methods are implemented in the AbstractTreeLikelihood class.
		 *
		 * @{
		 */
		ParameterList getBranchLengthsParameters() const;
		ParameterList getSubstitutionModelParameters() const;
		/** @} */

		/**
		 * @brief Get the substitution model used for the computation.
		 *
		 * @return A const pointer toward the substitution model of this instance.
		 */
		virtual const SubstitutionModel * getSubstitutionModel() const { return _model; }
		
		/**
		 * @brief Get the substitution model used for the computation.
		 *
		 * @return A pointer toward the substitution model of this instance.
		 */
		virtual SubstitutionModel * getSubstitutionModel() { return _model; }
		
	public: //Specific methods:

		/**
		 * @brief This builds the <i>parameters</i> list from all parametrizable objects,
		 * <i>i.e.</i> substitution model, rate distribution and tree.
		 */
		virtual void initParameters();

		/**
		 * @brief This removes a particular parameter from the list.
		 *
		 * This method may be used to not estimate a parameter after having
		 * fixed its value. The previous method reset all calls of thos one.
		 *
		 * @param name The name of the parameter to ignore.
		 */
		virtual void ignoreParameter(const string & name) throw (ParameterNotFoundException);

		/**
		 * @brief All parameters are stores in a parameter list.
		 *
		 * This function apply these parameters to the substitution model,
		 * to the rate distribution and to the branch lengths.
		 */
		virtual void applyParameters() throw (Exception);	

		virtual void initBranchLengthsParameters();

		void resetLikelihoodArray(VVVdouble & likelihoodArray);

		static void displayLikelihoodArray(const VVVdouble & likelihoodArray);

};

#endif //_ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_

