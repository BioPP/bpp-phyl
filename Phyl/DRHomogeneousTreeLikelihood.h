//
// File: DRHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Fri Oct 17 18:14:51 2003
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

#ifndef _DRHOMOGENEOUSTREELIKELIHOOD_H_
#define _DRHOMOGENEOUSTREELIKELIHOOD_H_

#include "AbstractHomogeneousTreeLikelihood.h"

// From NumCalc:
#include <NumCalc/VectorTools.h>
#include <NumCalc/DiscreteDistribution.h>

// From the STL:
#include <map>

using namespace std;

/**
 * @brief This class implement the computation of likelihood for a tree using the double-recursive
 * method.
 *
 * The substitution model is constant over the tree (homogeneous model).
 * A non uniform distribution of rates among the sites is allowed (ASRV models).</p>
 *
 * The Felsenstein recursive algorithm is used for conputation.
 * As in the HomogeneousTreeLikelihood class, a likelihood tensor is defined:
 * 
 * -Site
 * -Rate class
 * -Ancestral state
 * 
 * However in this class, a node will be attached a set of tensor instead of one single tensor,
 * one tensor for each subtree it defines.
 * These tensors are stored into a map of map, with each node as a primary key and each neighbor
 * of the node (sons + father) a a secondary key.
 *
 * No pattern are used here.
 */
class DRHomogeneousTreeLikelihood : public virtual AbstractHomogeneousTreeLikelihood
{
	protected:
		SiteContainer * _shrunkData;

		/**
		 * @brief This contains all likelihood values used for computation.
		 *
		 * <pre>
		 * x[n][b][i][c][s]
		 *   |---------------> Node n (pointer)
		 *      |------------> Neighbor node of n (pointer)
		 *         |---------> Site i
		 *            |------> Rate class c
		 *               |---> Ancestral state s
		 * </pre> 
		 * We call this the <i>likelihood array</i> for each node.
		 */
		mutable map<const Node *, map<const Node *, VVVdouble> > _likelihoods;

		/**
		 * @brief This contains all likelihood first order derivatives values used for computation.
		 *
		 * <pre>
		 * x[b][i]
		 *   |------------> Neighbor node of n (pointer)
		 *      |---------> Site i
		 * </pre> 
		 * We call this the <i>dLikelihood array</i> for each node.
		 */
		mutable map<const Node *, Vdouble> _dLikelihoods;
	
		/**
		 * @brief This contains all likelihood second order derivatives values used for computation.
		 *
		 * <pre>
		 * x[b][i]
		 *   |------------> Neighbor node of n (pointer)
		 *      |---------> Site i
		 * </pre> 
		 * We call this the <i>d2Likelihood array</i> for each node.
		 */
		mutable map<const Node *, Vdouble> _d2Likelihoods;
	
		mutable map<const Node *, VVdouble> _leavesLikelihoods;

		mutable VVVdouble _rootLikelihoods;
		mutable VVdouble  _rootLikelihoodsS;
		mutable Vdouble   _rootLikelihoodsSR;

		//some values we'll need:
		unsigned int _nbDistinctSites; //the number of distinct sites in the container
		
	public:
		DRHomogeneousTreeLikelihood(
			TreeTemplate<Node> & tree,
			const SiteContainer & data,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
			bool verbose = true)
			throw (Exception);
	
		virtual ~DRHomogeneousTreeLikelihood();
	
	public:

		/**
		 * @name The TreeLikelihood interface.
		 *
		 * Other methods are implemented in the AbstractTreeLikelihood class.
		 *
		 * @{
		 */
		double getLikelihood () const;
		double getLogLikelihood() const;
		double getLikelihoodForASite (unsigned int site) const;
		double getLogLikelihoodForASite(unsigned int site) const;
		/** @} */

		void computeTreeLikelihood();

		
		/**
		 * @name The DiscreteRatesAcrossSites interface implementation:
		 *
		 * @{
		 */
		double getLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
		double getLogLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
		double getLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const;
		double getLogLikelihoodForASiteForARateClassForAState(unsigned int site, unsigned int rateClass, int state) const;
		/** @} */
	
		/**
		 * @brief Implements the Function interface.
		 *
		 * Update the parameter list and call the applyParameters() method.
		 * Then compute the likelihoods at each node (computeLikelihood() method)
		 * and call the getLogLikelihood() method.
		 *
		 * If a subset of the whole parameter list is passed to the function,
		 * only these parameters are updated and the other remain constant (i.e.
		 * equal to their last value).
		 *
		 * @param parameters The parameter list to pass to the function.
		 */
		void setParameters(const ParameterList & parameters) throw (ParameterNotFoundException, ConstraintException);
		double getValue() const throw(Exception);
		
		/**
		 * @name DerivableFirstOrder interface.
		 *
		 * @{
		 */
		double getFirstOrderDerivative(const string & variable) const throw (Exception);
		/** @{ */

		/**
		 * @name DerivableSecondOrder interface.
		 *
		 * @{
		 */
		double getSecondOrderDerivative(const string & variable) const throw (Exception);
		double getSecondOrderDerivative(const string & variable1, const string & variable2) const throw (Exception) { return 0; } // Not implemented for now.
		/** @} */
		
	public:	// Specific methods:
	
		virtual map<const Node *, VVVdouble> getLikelihoodArraysForNode(const Node * node)
	  { return _likelihoods[node]; }

		virtual VVVdouble getTransitionProbabilitiesForNode(const Node * node)
		{ return _pxy[node]; }

		virtual vector<unsigned int> getRootPatternLinks() const
		{ return _rootPatternLinks; }

		virtual VVVdouble getRootLikelihoods() const
		{ return _rootLikelihoods; }

		virtual unsigned int getNumberOfDistinctSites() const
		{ return _nbDistinctSites; }

		virtual VVVdouble computeLikelihoodAtNode(const Node * node) const;
		
		virtual VVdouble getLeafLikelihoods(const Node * node) const
		{ return _leavesLikelihoods[node]; }

		virtual const SiteContainer * getShrunkData() const
		{ return _shrunkData; }

	protected:
		
		/**
		 * @brief This method initializes the leaves according to a sequence container.
		 *
		 * Here the container _shrunkData is used.
		 * Likelihood is set to 1 for the state corresponding to the sequence site,
		 * otherwise it is set to 0.
		 *
		 * All likelihood arrays at each nodes are initialized according to alphabet
		 * size and sequences length, and filled with 1.
		 *
		 * NB: This method is recursive.
		 *
		 * @param node      The node defining the subtree to analyse.
		 * @param sequences The sequence container to use.
		 */
		virtual void initTreeLikelihoods(const Node * node, const SequenceContainer & sequences) throw (Exception);

		/**
		 * Initialize the arrays corresponding to each son node for the node passed as argument.
		 * The method is called for each son node and the result stored in the corresponding array.
		 */
		virtual void computeSubtreeLikelihoodPostfix(const Node * node); //Recursive method.
		/**
		 * This method initilize the remaining lieklihood arrays, corresponding to father nodes.
		 * It must be called after the postfix method because it requires that the arrays for
		 * son nodes to be be computed.
		 */
		virtual void computeSubtreeLikelihoodPrefix(const Node * node); //Recursive method.

		virtual void computeRootLikelihood();

		virtual void computeTreeDLikelihoodAtNode(const Node * node);
		virtual void computeTreeDLikelihoods();
		
		virtual void computeTreeD2LikelihoodAtNode(const Node * node);
		virtual void computeTreeD2Likelihoods();

		void fireParameterChanged(const ParameterList & params);

		void resetLikelihoodArrays(const Node * node);
	
		/**
		 * @brief This method is mainly for debugging purpose.
		 *
		 * @param node The node at which likelihood values must be displayed.
		 */
		virtual void displayLikelihood(const Node * node);

};


#endif	//_DRHOMOGENEOUSTREELIKELIHOOD_H_
