//
// File: RTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Tue May 15 14:30 2012
// From file: RNonHomogeneousTreeLikelihood.h
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _RTREELIKELIHOOD_H_
#define _RTREELIKELIHOOD_H_

#include "AbstractTreeLikelihood.h"
#include "SubstitutionProcess.h"
#include "RTreeLikelihoodData.h"

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp
{
namespace newlik
{

  /**
   * @brief This class implements the single recursion likelihood computation for a tree.
   *
   * This class uses an instance of the RTreeLikelihoodData for conditionnal likelihood storage.
   *
   * This class can also use a simple or recursive site compression.
   * In the simple case, computations for identical sites are not duplicated.
   * In the recursive case, computations for identical sub-sites (<i>site patterns </i>) are also not duplicated:
   * Following N. Galtier (personal communication), we define a Pattern as a distinct site
   * in a sub-dataset corresponding to the dataset with sequences associated to a particular subtree.
   * The likelihood computation is the same for a given site, hence the idea is to save time from
   * performing many times the same coputation.
   * The network between all patterns is defined by the _patternLinks double map, initialized in the
   * initLikelihoodsWithPatterns() method. This initialisation takes more time than the classic
   * initTreeLikelihood one, where all likelihoods for a given site <i>i</i> are at the <i>i</i> coordinate
   * in the likelihood tensor, but is really faster when computing the likelihood (computeLikelihoods() method).
   * Hence, if you have to compute likelihood many times while holding the tree topology unchanged,
   * you should use patterns.
   * This decreases the likelihood computation time, but at a cost: some time is spent to establish the patterns
   * relationships. Whether to use or not patterns depends on what you actllay need:
   * - The more you compute likelihoods without changing the data or topology, the more patterns are interesting
   *   (this divides the cost of computing patterns by the number of computation performed).
   *   Patterns are hence usefull when you have a high number of computation to perform, while optimizing numerical
   *   parameters for instance).
   * - Patterns are more likely to occur whith small alphabet (nucleotides).
   *
   * Important note: The input tree will be considered as rooted, since the likelihood of non-stationary models
   * depends on the position of the root. If the input tree is not rooted, it will be considered as a rotted tree
   * with a root multifurcation.
   */
  class RTreeLikelihood :
    public AbstractTreeLikelihood
  {
  private:

    mutable std::auto_ptr<RTreeLikelihoodData> likelihoodData_;
    double minusLogLik_;

  public:
    /**
     * @brief Build a new Simple Recursive Tree Likelihood object without data.
     *
     * This constructor only initialize the parameters.
     * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
     *
     * @param tree The tree to use.
     * @param modelSet The set of substitution models to use.
     * @param rDist The rate across sites distribution to use.
     * If true, any rooted tree will be unrooted before likelihood computation.
     * @param verbose Should I display some info?
     * @param usePatterns Tell if recursive site compression should be performed.
     * @throw Exception in an error occured.
     */
    RTreeLikelihood(
        SubstitutionProcess* process,
        bool verbose = true,
        bool usePatterns = true)
      throw (Exception);
	
    /**
     * @brief Build a new Simple Recursive Tree Likelihood object and compute the corresponding likelihood.
     *
     * This constructor initializes all parameters, data, and likelihood arrays.
     *
     * @param data Sequences to use.
     * @param modelSet The set of substitution models to use.
     * @param rDist The rate across sites distribution to use.
     * @param verbose Should I display some info?
     * @param usePatterns Tell if recursive site compression should be performed.
     * @throw Exception in an error occured.
     */
    RTreeLikelihood(
        const SiteContainer& data,
        SubstitutionProcess* process,
        bool verbose = true,
        bool usePatterns = true)
      throw (Exception);

    RTreeLikelihood(const RTreeLikelihood& lik);
    
    RTreeLikelihood& operator=(const RTreeLikelihood& lik);

    virtual ~RTreeLikelihood() {} //smart pointers take care of everything.

    RTreeLikelihood* clone() const { return new RTreeLikelihood(*this); }
	
  private:

    /**
     * @brief Method called by constructors.
     */
    void init_(bool usePatterns) throw (Exception);
	
  public:

    /**
     * @name The TreeLikelihood interface.
     *
     * Other methods are implemented in the AbstractHomogeneousTreeLikelihood class.
     *
     * @{
     */
    void setData(const SiteContainer& sites) throw (Exception);
    unsigned int getSiteIndex(unsigned int site) const throw (IndexOutOfBoundsException) { return likelihoodData_->getRootArrayPosition(site); }
		
    double getLogLikelihood() const;
    double getLogLikelihoodForASite(unsigned int site) const;
    double getLikelihoodForASiteForAClass(unsigned int site, unsigned int modelClass) const;
    double getLikelihoodForASiteForAState(unsigned int site, int state) const;
    double getLikelihoodForASiteForAClassForAState(unsigned int site, unsigned int modelClass, int state) const;
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
    void setParameters(const ParameterList& parameters) throw (ParameterNotFoundException, ConstraintException);
    double getValue() const throw(Exception);
		
    /**
     * @name DerivableFirstOrder interface.
     *
     * @{
     */
    double getFirstOrderDerivative(const std::string& variable) const throw (Exception);
    /** @} */

    /**
     * @name DerivableSecondOrder interface.
     *
     * @{
     */
    double getSecondOrderDerivative(const std::string& variable) const throw (Exception);
    double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const throw (Exception) { return 0; } // Not implemented for now.
    /** @} */
	
  public:	// Specific methods:
	
    RTreeLikelihoodData* getLikelihoodData() { return likelihoodData_.get(); }
    const RTreeLikelihoodData* getLikelihoodData() const { return likelihoodData_.get(); }
 
    virtual void computeTreeLikelihood();

    virtual double getDLikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
    virtual double getDLikelihoodForASite(unsigned int site) const;
    virtual double getDLogLikelihoodForASite(unsigned int site) const;
    virtual double getDLogLikelihood() const;
    
    virtual void computeTreeDLikelihood(const std::string& variable);

    virtual double getD2LikelihoodForASiteForARateClass(unsigned int site, unsigned int rateClass) const;
    virtual double getD2LikelihoodForASite(unsigned int site) const;
    virtual double getD2LogLikelihoodForASite(unsigned int site) const;
    virtual double getD2LogLikelihood() const;
		
    virtual void computeTreeD2Likelihood(const std::string& variable);

	
  protected:
			
    /**
     * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
     *
     * @param node The root of the subtree.
     */
    virtual void computeSubtreeLikelihood(const Node * node); //Recursive method.			
    virtual void computeDownSubtreeDLikelihood(const Node *);
    virtual void computeDownSubtreeD2Likelihood(const Node *);
	
    void fireParameterChanged(const ParameterList & params);
	
    /**
     * @brief This method is mainly for debugging purpose.
     *
     * @param node The node at which likelihood values must be displayed.
     */
    virtual void displayLikelihood(const Node* node);

  };


} //end of namespace newlik.
} //end of namespace bpp.

#endif	//_RTREELIKELIHOOD_H_

