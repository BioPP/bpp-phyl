//
// File: DoubleRecursiveTreeLikelihoodCalculation.h
// Created by: Laurent Guéguen
// Created on: lundi 21 juillet 2014, à 09h 19
// From file: DoubleRecursiveNonHomogeneousTreeLikelihood.h
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

#ifndef _DOUBLERECURSIVETREELIKELIHOODCALCULATION_H_
#define _DOUBLERECURSIVETREELIKELIHOODCALCULATION_H_

#include "AbstractTreeLikelihoodCalculation.h"
#include "DoubleRecursiveTreeLikelihoodData.h"

#include <Bpp/Numeric/VectorTools.h>
//#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp
{
namespace newlik
{

/**
 * @brief This class implements the likelihood computation for a tree
 * using the double-recursive algorithm.
 *
 * This class uses an instance of the
 * DoubleRecursiveTreeLikelihoodData for conditionnal likelihood
 * storage.
 *
 * All nodes share the same site patterns.
 *
 * Important note: The input tree will be considered as rooted, since
 * the likelihood of non-stationary models depends on the position of
 * the root. If the input tree is not rooted, it will be considered as
 * a rooted tree with a root multifurcation.
 *
 */

  class DoubleRecursiveTreeLikelihoodCalculation:
    public AbstractTreeLikelihoodCalculation
  {
private:
    mutable std::auto_ptr<DoubleRecursiveTreeLikelihoodData> likelihoodData_;
    int root1_, root2_; // Needed only in case of reparametrization of branch length at root node.
    // TODO: have to be initialized properly! We do not care of that for now. jdutheil on 11/12/12.

    // booleans to say if the Dlikelihoods are null
  
    bool nullDLogLikelihood_;
    bool nullD2LogLikelihood_;


    // Node being currently derivated

    int compNId_;
    
  public:
    /**
     * @brief Build a new DoubleRecursiveTreeLikelihoodCalculation
     * object without data.
     *
     * This constructor only initialize the parameters.
     * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
     *
     * @param process The substitution process to use.
     * @param verbose Should I display some info?
     * @throw Exception in an error occured.
     */
    
    DoubleRecursiveTreeLikelihoodCalculation(
      SubstitutionProcess* process,
      bool verbose = true)
      throw (Exception);
  
    /**
     * @brief Build a new DoubleRecursiveTreeLikelihoodCalculation
     * object and compute the corresponding likelihood.
     *
     * This constructor initializes all parameters, data, and likelihood arrays.
     *
     * @param data Sequences to use.
     * @param process The substitution process to use.
     * @param verbose Should I display some info?
     * @throw Exception in an error occured.
     */
    DoubleRecursiveTreeLikelihoodCalculation(
      const SiteContainer& data,
      SubstitutionProcess* process,
      bool verbose = true)
      throw (Exception);

    /**
     * @brief Copy constructor.
     */ 
    DoubleRecursiveTreeLikelihoodCalculation(const DoubleRecursiveTreeLikelihoodCalculation& lik);
    
    DoubleRecursiveTreeLikelihoodCalculation& operator=(const DoubleRecursiveTreeLikelihoodCalculation& lik);

    virtual ~DoubleRecursiveTreeLikelihoodCalculation() {} //smart pointers take care of everything.

    DoubleRecursiveTreeLikelihoodCalculation* clone() const { return new DoubleRecursiveTreeLikelihoodCalculation(*this); }

  private:

    /**
     * @brief Method called by constructors.
     */
    void init_() throw (Exception);

  public:

    const Alphabet* getAlphabet() const throw (TreeLikelihoodCalculationNotInitializedException)
    {
      if (!initialized_)
        throw new TreeLikelihoodCalculationNotInitializedException("DoubleRecursiveTreeLikelihoodCalculation::getAlphabet().");
      return data_->getAlphabet();
    }

    size_t getNumberOfSites() const throw (TreeLikelihoodCalculationNotInitializedException) {
      if (!initialized_)
        throw new TreeLikelihoodCalculationNotInitializedException("DoubleRecursiveTreeLikelihoodCalculation::getNumberOfSites().");
      return data_->getNumberOfSites();
    }

    bool isInitialized() const { return initialized_; }

    void setData(const SiteContainer& sites) throw (Exception);
  
    const SiteContainer* getData() const throw (TreeLikelihoodCalculationNotInitializedException) {
      if (!initialized_)
        throw new TreeLikelihoodCalculationNotInitializedException("DoubleRecursiveTreeLikelihoodCalculation::getData().");
      return data_.get();
    }
  
    size_t getSiteIndex(size_t site) const throw (IndexOutOfBoundsException) {
      return likelihoodData_->getRootArrayPosition(site);
    }
    
    DoubleRecursiveTreeLikelihoodData* getLikelihoodData() { return likelihoodData_.get(); }
    
    const DoubleRecursiveTreeLikelihoodData* getLikelihoodData() const { return likelihoodData_.get(); }
  
    //  double getLogLikelihood() const;
    double getLikelihoodForASite(size_t site) const;

    double getLikelihoodForASiteForAState(size_t site, int state) const;

    double getLikelihoodForASiteForAClass(size_t site, size_t modelClass) const;

    double getLikelihoodForASiteForAClassForAState(size_t site, size_t modelClass, int state) const;
 
    double getDLogLikelihoodForASite(size_t site) const;

    double getD2LogLikelihoodForASite(size_t site) const;

    void computeTreeLikelihood();

    void computeTreeDLogLikelihood(const std::string& variable);

    void computeTreeD2LogLikelihood(const std::string& variable);

  private:
    void computeLikelihoodAtNode_(const Node* node, VVVdouble& likelihoodArray) const;
  
    /**
     * Initialize the arrays corresponding to each son node for the
     * node passed as argument. The method is called for each son node
     * and the result stored in the corresponding array.
     */

    void computeSubtreeLikelihoodPostfix_(const Node* node); //Recursive method.

    /**
     * This method initilize the remaining likelihood arrays,
     * corresponding to father nodes. It must be called after the
     * postfix method because it requires that the arrays for son
     * nodes to be be computed.
     */

    void computeSubtreeLikelihoodPrefix_(const Node* node); //Recursive method.

    //void computeRootLikelihood();

    void computeTreeDLogLikelihoodAtNode(const Node* node);
//    void computeTreeDLikelihoods();
    
    void computeTreeD2LogLikelihoodAtNode(const Node* node);
    //  void computeTreeD2Likelihoods();

//    void resetLikelihoodArrays(const Node* node);
    void computeRootLikelihood();
  
    /**
     * @brief This method is mainly for debugging purpose.
     *
     * @param node The node at which likelihood values must be displayed.
     */
    virtual void displayLikelihood(const Node* node);

};

} //end of namespace newlik.
} //end of namespace bpp.

#endif  //_DOUBLERECURSIVETREELIKELIHOOD_H_

