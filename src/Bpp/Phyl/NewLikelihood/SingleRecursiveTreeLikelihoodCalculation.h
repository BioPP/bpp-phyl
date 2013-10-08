//
// File: SingleRecursiveTreeLikelihoodCalculation.h
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

#ifndef _SINGLERECURSIVETREELIKELIHOODCALCULATION_H_
#define _SINGLERECURSIVETREELIKELIHOODCALCULATION_H_

#include "AbstractTreeLikelihoodCalculation.h"
#include "SingleRecursiveTreeLikelihoodData.h"

#include <Bpp/Numeric/VectorTools.h>
//#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp
{
namespace newlik
{
/**
 * @brief This class implements the single recursion likelihood computation for a tree.
 *
 * This class uses an instance of the SingleRecursiveTreeLikelihoodData for conditionnal likelihood storage.
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
 * depends on the position of the root. If the input tree is not rooted, it will be considered as a rooted tree
 * with a root multifurcation.
 */
class SingleRecursiveTreeLikelihoodCalculation:
  public AbstractTreeLikelihoodCalculation
{
private:
  mutable std::auto_ptr<SingleRecursiveTreeLikelihoodData> likelihoodData_;
  int root1_, root2_; // Needed only in case of reparametrization of branch length at root node.
  // TODO: have to be initialized properly! We do not care of that for now. jdutheil on 11/12/12.

  // booleans to say if the Dlikelihoods are null
  
  bool nullDLikelihood_;
  bool nullD2Likelihood_;
  
public:
  /**
   * @brief Build a new Simple Recursive Tree Likelihood object without data.
   *
   * This constructor only initialize the parameters.
   * To compute a likelihood, you will need to call the setData() and the computeTreeLikelihood() methods.
   *
   * @param process The substitution process to use.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @throw Exception in an error occured.
   */
  SingleRecursiveTreeLikelihoodCalculation(
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
   * @param process The substitution process to use.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @throw Exception in an error occured.
   */
  SingleRecursiveTreeLikelihoodCalculation(
    const SiteContainer& data,
    SubstitutionProcess* process,
    bool verbose = true,
    bool usePatterns = true)
  throw (Exception);

  SingleRecursiveTreeLikelihoodCalculation(const SingleRecursiveTreeLikelihoodCalculation& lik);

  SingleRecursiveTreeLikelihoodCalculation& operator=(const SingleRecursiveTreeLikelihoodCalculation& lik);

  virtual ~SingleRecursiveTreeLikelihoodCalculation() {} // smart pointers take care of everything.

  SingleRecursiveTreeLikelihoodCalculation* clone() const { return new SingleRecursiveTreeLikelihoodCalculation(*this); }

private:
  /**
   * @brief Method called by constructors.
   */
  void init_(bool usePatterns) throw (Exception);

public:

  const Alphabet* getAlphabet() const throw (TreeLikelihoodCalculationNotInitializedException)
  {
    if (!initialized_)
      throw new TreeLikelihoodCalculationNotInitializedException("SingleRecursiveTreeLikelihoodCalculation::getAlphabet().");
    return data_->getAlphabet();
  }

  size_t getNumberOfSites() const throw (TreeLikelihoodCalculationNotInitializedException) {
    if (!initialized_)
      throw new TreeLikelihoodCalculationNotInitializedException("SingleRecursiveTreeLikelihoodCalculation::getNumberOfSites().");
    return data_->getNumberOfSites();
  }

  bool isInitialized() const { return initialized_; }

  void setData(const SiteContainer& sites) throw (Exception);
  
  const SiteContainer* getData() const throw (TreeLikelihoodCalculationNotInitializedException) {
    if (!initialized_)
      throw new TreeLikelihoodCalculationNotInitializedException("SingleRecursiveTreeLikelihoodCalculation::getData().");
    return data_.get();
  }
  
  size_t getSiteIndex(size_t site) const throw (TreeLikelihoodCalculationNotInitializedException, IndexOutOfBoundsException) {
    if (!initialized_)
      throw new TreeLikelihoodCalculationNotInitializedException("SingleRecursiveTreeLikelihoodCalculation::getSiteIndex().");
    return likelihoodData_->getRootArrayPosition(site);
  }

  SingleRecursiveTreeLikelihoodData* getLikelihoodData() { return likelihoodData_.get(); }
  const SingleRecursiveTreeLikelihoodData* getLikelihoodData() const { return likelihoodData_.get(); }

  double getLikelihoodForASite(size_t site) const;

  double getLikelihoodForASiteForAState(size_t site, int state) const;

  double getLikelihoodForASiteForAClass(size_t site, size_t classIndex) const;

  double getLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state) const;

  double getDLogLikelihood() const;
  
  double getD2LogLikelihood() const;
 
  double getDLikelihoodForASite(size_t site) const;
  
  double getD2LikelihoodForASite(size_t site) const;
  
  void computeTreeLikelihood();
  void computeTreeDLikelihood(const std::string& variable);
  void computeTreeD2Likelihood(const std::string& variable);
 
public:
  // Specific methods:
  double getDLogLikelihoodForASite(size_t site) const;
  double getDLikelihoodForASiteForAClass(size_t site, size_t classIndex) const;
  double getDLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state) const;

  double getD2LogLikelihoodForASite(size_t site) const;
  double getD2LikelihoodForASiteForAClass(size_t site, size_t classIndex) const;
  double getD2LikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state) const;

protected:

  /**
   * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
   *
   * @param node The root of the subtree.
   */
  void computeSubtreeLikelihood_(const Node* node); // Recursive method.
  void computeDownSubtreeDLikelihood_(const Node* node);
  void computeDownSubtreeD2Likelihood_(const Node* node);

  /**
   * @brief This method is mainly for debugging purpose.
   *
   * @param node The node at which likelihood values must be displayed.
   */
  void displayLikelihood(const Node* node);
};
} // end of namespace newlik.
} // end of namespace bpp.

#endif  // _SINGLERECURSIVETREELIKELIHOODCALCULATION_H_

