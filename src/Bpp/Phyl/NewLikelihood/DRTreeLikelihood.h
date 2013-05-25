//
// File: DRTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Apr 26 20:07 2013
// From file: DRNonHomogeneousTreeLikelihood.h
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

#ifndef _NEWLIK_DRTREELIKELIHOOD_H_
#define _NEWLIK_DRTREELIKELIHOOD_H_

#include "AbstractTreeLikelihood.h"
#include "DRTreeLikelihoodData.h"

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp
{
namespace newlik
{

/**
 * @brief This class implements the likelihood computation for a tree using the double-recursive
 * algorithm, allowing for non-homogeneous models of substitutions.
 *
 * The substitution model is the same over the tree (homogeneous model).
 * A non-uniform distribution of rates among the sites is allowed (ASRV models).</p>
 *
 * This class uses an instance of the DRASDRTreeLikelihoodData for conditionnal likelihood storage.
 *
 * All nodes share the same site patterns.
 *
 * Important note: The input tree will be considered as rooted, since the likelihood of non-stationary models
 * depends on the position of the root. If the input tree is not rooted, it will be considered as a rotted tree
 * with a root multifurcation.
 */
class DRTreeLikelihood:
  public AbstractTreeLikelihood
{
  private:
    mutable std::auto_ptr<DRTreeLikelihoodData> likelihoodData_;
    double minusLogLik_;
   
  public:
    /**
     * @brief Build a new Double Recursive Tree Likelihood object without data.
     *
     * This constructor only initialize the parameters.
     * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
     *
     * @param process The substitution process to use.
     * @param verbose Should I display some info?
     * @throw Exception in an error occured.
     */
    DRTreeLikelihood(
      SubstitutionProcess* process,
      bool verbose = true)
      throw (Exception);
  
    /**
     * @brief Build a new DRNonHomogeneousTreeLikelihood object and compute the corresponding likelihood.
     *
     * This constructor initializes all parameters, data, and likelihood arrays.
     *
     * @param data Sequences to use.
     * @param process The substitution process to use.
     * @param verbose Should I display some info?
     * @throw Exception in an error occured.
     */
    DRTreeLikelihood(
      const SiteContainer& data,
      SubstitutionProcess* process,
      bool verbose = true)
      throw (Exception);

    /**
     * @brief Copy constructor.
     */ 
    DRTreeLikelihood(const DRTreeLikelihood& lik);
    
    DRTreeLikelihood& operator=(const DRTreeLikelihood& lik);

    virtual ~DRTreeLikelihood() {} //smart pointers take care of everything.

    DRTreeLikelihood* clone() const { return new DRTreeLikelihood(*this); }

  private:

    /**
     * @brief Method called by constructors.
     */
    void init_() throw (Exception);

  public:

    /**
     * @name The TreeLikelihood interface.
     *
     * Other methods are implemented in the AbstractTreeLikelihood class.
     *
     * @{
     */
    void setData(const SiteContainer& sites) throw (Exception);
    
    size_t getSiteIndex(size_t site) const throw (IndexOutOfBoundsException) { return likelihoodData_->getRootArrayPosition(site); }
    
    DRTreeLikelihoodData* getLikelihoodData() { return likelihoodData_.get(); }
    const DRTreeLikelihoodData* getLikelihoodData() const { return likelihoodData_.get(); }
  
    double getLogLikelihood() const;
    double getLikelihoodForASite(size_t site) const;
    double getLikelihoodForASiteForAClass(size_t site, size_t modelClass) const;
    double getLikelihoodForASiteForAState(size_t site, int state) const;
    double getLikelihoodForASiteForAClassForAState(size_t site, size_t modelClass, int state) const;
    /** @} */

    /**
     * @brief Implements the Function interface.
     *
     * Update the parameter list and call the fireParameterChanged() method.
     *
     * If a subset of the whole parameter list is passed to the function,
     * only these parameters are updated and the other remain constant (i.e.
     * equal to their last value).
     *
     * @param parameters The parameter list to pass to the function.
     */
    void setParameters(const ParameterList& parameters) throw (ParameterNotFoundException, ConstraintException);
    
    /**
     * @brief Function and NNISearchable interface.
     */
    double getValue() const throw (Exception);
    
    /**
     * @name DerivableFirstOrder interface.
     *
     * @{
     */
    double getFirstOrderDerivative(const std::string& variable) const throw (Exception);
    /** @{ */

    /**
     * @name DerivableSecondOrder interface.
     *
     * @{
     */
    double getSecondOrderDerivative(const std::string& variable) const throw (Exception);
    double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const throw (Exception) { return 0; } // Not implemented for now.
    /** @} */
    
  public:  // Specific methods:

    virtual void computeTreeLikelihood();
    
    virtual void computeLikelihoodAtNode(int nodeId, VVVdouble& likelihoodArray) const
    {
      computeLikelihoodAtNode_(tree_->getNode(nodeId), likelihoodArray);
    }
      
  protected:
    virtual void computeLikelihoodAtNode_(const Node* node, VVVdouble& likelihoodArray) const;

  
    /**
     * Initialize the arrays corresponding to each son node for the node passed as argument.
     * The method is called for each son node and the result stored in the corresponding array.
     */
    virtual void computeSubtreeLikelihoodPostfix(const Node* node); //Recursive method.
    /**
     * This method initilize the remaining likelihood arrays, corresponding to father nodes.
     * It must be called after the postfix method because it requires that the arrays for
     * son nodes to be be computed.
     */
    virtual void computeSubtreeLikelihoodPrefix(const Node* node); //Recursive method.

    virtual void computeRootLikelihood();

    virtual void computeTreeDLikelihoodAtNode(const Node* node);
    virtual void computeTreeDLikelihoods();
    
    virtual void computeTreeD2LikelihoodAtNode(const Node* node);
    virtual void computeTreeD2Likelihoods();

    void fireParameterChanged(const ParameterList& params);

    void resetLikelihoodArrays(const Node* node);
  
    /**
     * @brief This method is mainly for debugging purpose.
     *
     * @param node The node at which likelihood values must be displayed.
     */
    virtual void displayLikelihood(const Node* node);

    /**
     * @brief Compute conditional likelihoods.
     *
     * This method is the "core" likelihood computation function, performing all the product uppon all nodes, the summation for each ancestral state and each rate class.
     * It is designed for inner usage, and a maximum efficiency, so no checking is performed on the input parameters.
     * Use with care!
     * 
     * @param iLik A vector of likelihood arrays, one for each conditional node.
     * @param tProb A vector of transition probabilities, one for each node.
     * @param oLik The likelihood array to store the computed likelihoods.
     * @param nbNodes The number of nodes = the size of the input vectors.
     * @param nbDistinctSites The number of distinct sites (the first dimension of the likelihood array).
     * @param nbClasses The number of rate classes (the second dimension of the likelihood array).
     * @param nbStates The number of states (the third dimension of the likelihood array).
     * @param reset Tell if the output likelihood array must be initalized prior to computation.
     * If true, the resetLikelihoodArray method will be called.
     */
    static void computeLikelihoodFromArrays(
        const std::vector<const VVVdouble*>& iLik,
        const std::vector<const VVVdouble*>& tProb,
        VVVdouble& oLik, size_t nbNodes,
        size_t nbDistinctSites,
        size_t nbClasses,
        size_t nbStates,
        bool reset = true);

    /**
     * @brief Compute conditional likelihoods.
     *
     * This method is the "core" likelihood computation function, performing all the product uppon all nodes, the summation for each ancestral state and each rate class.
     * This function is specific to non-reversible models: the subtree containing the root is specified separately.
     * It is designed for inner usage, and a maximum efficiency, so no checking is performed on the input parameters.
     * Use with care!
     * 
     * @param iLik A vector of likelihood arrays, one for each conditional node.
     * @param tProb A vector of transition probabilities, one for each node.
     * @param iLikR The likelihood array for the subtree containing the root of the tree.
     * @param tProbR The transition probabilities for thr subtree containing the root of the tree.
     * @param oLik The likelihood array to store the computed likelihoods.
     * @param nbNodes The number of nodes = the size of the input vectors.
     * @param nbDistinctSites The number of distinct sites (the first dimension of the likelihood array).
     * @param nbClasses The number of rate classes (the second dimension of the likelihood array).
     * @param nbStates The number of states (the third dimension of the likelihood array).
     * @param reset Tell if the output likelihood array must be initalized prior to computation.
     * If true, the resetLikelihoodArray method will be called.
     */
    static void computeLikelihoodFromArrays(
        const std::vector<const VVVdouble*>& iLik,
        const std::vector<const VVVdouble*>& tProb,
        const VVVdouble* iLikR,
        const VVVdouble* tProbR,
        VVVdouble& oLik,
        size_t nbNodes,
        size_t nbDistinctSites,
        size_t nbClasses,
        size_t nbStates,
        bool reset = true);

};

} //end of namespace newlik.
} //end of namespace bpp.

#endif  //_NEWLIK_DRTREELIKELIHOOD_H_

