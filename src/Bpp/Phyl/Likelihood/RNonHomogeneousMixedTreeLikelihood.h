//
// File: RNonHomogeneousMixedLikelihood.h
// Created by: Laurent Gueguen
// Created on: jeudi 11 novembre 2010, à 07h 56
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _RNONHOMOGENEOUSMIXEDTREELIKELIHOOD_H_
#define _RNONHOMOGENEOUSMIXEDTREELIKELIHOOD_H_

#include "RNonHomogeneousTreeLikelihood.h"
#include "../Model/SubstitutionModelSet.h"

#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

using namespace std;
namespace bpp
{
/**
 *@ brief A class to compute the average of several
 *RNonHomogeneousTreeLikelihood defined from a Mixed Substitution
 *Model.
 *
 * This class is made recursively. At each node, we test if an
 * expansion of a mixed model is necessary. This is the case when this
 * model points towards different subtrees under this node, or towards
 * a son of this node and a branch under it. If an expansion is
 * necessary, a vector of RNonHomogeneousMixedLikelihood* is built
 * with all the submodels combinations.
 *
 * Note that this approach is not the most efficient, since a graph
 * based one would avoid some computations, but it seems much more
 * difficult to do it in the extant hierarchy.
 **/

class RNonHomogeneousMixedTreeLikelihood :
  public RNonHomogeneousTreeLikelihood
{
private:

  /**
   * @brief the map of the branch numbers to the vectors of the
   * TreeLikelihoods for the expanded model on this branch.
   *
   */
  map<int, vector<RNonHomogeneousMixedTreeLikelihood*> > mvTreeLikelihoods_;


  /**
   * @brief the map of the branch numbers to the vector of the probabilities of the
   * TreeLikelihoods for the expanded model on this branch.
   *
   */
  
  map<int, vector<double> > mvProbas_;

  /**
   * @brief the vector of the number of the submodels that are taken
   * into account in each mixed model. If the model is not mixed, the
   * associated value is 0. If the model is not expanded at that
   * point, the associated value is -1.
   *
   */
  
  vector<int>   vNumModels_;

  /**
   * @brief the number of the number of the node under which tree the
   * Treelikelihood is computed.
   *
   */
  
  int  upperNode_;

  /**
   * @brief a flag to say if this object is the head of the hierarchy
   *
   **/

  bool main_;
  
  /**
   * @brief Build a new RNonHomogeneousMixeTreeLikelihood object
   * without data.
   *
   * This constructor only initialize the parameters. To compute a
   * likelihood, you will need to call the setData() and the
   * computeTreeLikelihood() methods.
   *
   * @param tree The tree to use.
   * @param modelSet The set of substitution models to use.
   * @param vsubmod the vector of the numbers of the sub models used in
   *  the mixed models.
   * @param upperNode the number of the node under which the treelikelihood
   *  is computed.
   * @param rDist The rate across sites distribution to use.
   *  If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @throw Exception in an error occured.
   */

  RNonHomogeneousMixedTreeLikelihood(const Tree& tree,
                                     SubstitutionModelSet* modelSet,
                                     std::vector<int> &vsubmod,
                                     int upperNode,
                                     DiscreteDistribution* rDist,
                                     bool verbose,
                                     bool usePatterns,
                                     bool main=false);

  /**
   * @brief Build a new RNonHomogeneousMixeTreeLikelihood object
   * with data.
   *
   * This constructor only initialize the parameters. To compute a
   * likelihood, you will need to call the setData() and the
   * computeTreeLikelihood() methods.
   *
   * @param tree The tree to use.
   * @param data Sequences to use.
   * @param modelSet The set of substitution models to use.
   * @param vsubmod the vector of the numbers of the sub models used in
   *  the mixed models.
   * @param upperNode the number of the node under which the treelikelihood
   *  is computed.
   * @param rDist The rate across sites distribution to use.
   *  If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @throw Exception in an error occured.
   */

  RNonHomogeneousMixedTreeLikelihood(const Tree& tree,
                                     const SiteContainer& data,
                                     SubstitutionModelSet* modelSet,
                                     std::vector<int> &vsubmod,
                                     int upperNode,
                                     DiscreteDistribution* rDist,
                                     bool verbose,
                                     bool usePatterns,
                                     bool main=false);


  /**
   * brief method where the recursive structure is built.
   *
   */
  
  void init(const Tree& tree,
            const SiteContainer* pdata,
            SubstitutionModelSet* modelSet,
            std::vector<int> &vsubmod,
            DiscreteDistribution* rDist,
            bool verbose,
            bool usePatterns);

  
public:
  /**
   * @brief Build a new RNonHomogeneousMixeTreeLikelihood object
   * without data.
   *
   * This constructor only initialize the parameters. To compute a
   * likelihood, you will need to call the setData() and the
   * computeTreeLikelihood() methods.
   *
   * @param tree The tree to use.
   * @param modelSet The set of substitution models to use.
   * @param rDist The rate across sites distribution to use.
   * If true, any rooted tree will be unrooted before likelihood computation.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @throw Exception in an error occured.
   */
  RNonHomogeneousMixedTreeLikelihood(
    const Tree& tree,
    SubstitutionModelSet* modelSet,
    DiscreteDistribution* rDist,
    bool verbose = true,
    bool usePatterns = true,
    bool main = true)
  throw (Exception);

  /**
   * @brief Build a new RNonHomogeneousMixedTreeLikelihood object
   * and compute the corresponding likelihood.
   *
   * This constructor initializes all parameters, data, and
   * likelihood arrays.
   *
   * @param tree The tree to use.
   * @param data Sequences to use.
   * @param modelSet The set of substitution models to use.
   * @param rDist The rate across sites distribution to use.
   * @param verbose Should I display some info?
   * @param usePatterns Tell if recursive site compression should be performed.
   * @throw Exception in an error occured.
   */
  RNonHomogeneousMixedTreeLikelihood(const Tree& tree,
                                     const SiteContainer& data,
                                     SubstitutionModelSet* modelSet,
                                     DiscreteDistribution* rDist,
                                     bool verbose = true,
                                     bool usePatterns = true,
                                     bool main = true)
    throw (Exception);

  RNonHomogeneousMixedTreeLikelihood(const RNonHomogeneousMixedTreeLikelihood& lik);

  RNonHomogeneousMixedTreeLikelihood& operator=(const RNonHomogeneousMixedTreeLikelihood& lik);

  virtual ~RNonHomogeneousMixedTreeLikelihood();

  RNonHomogeneousMixedTreeLikelihood* clone() const { return new RNonHomogeneousMixedTreeLikelihood(*this); }

public:
  /**
   * @name The TreeLikelihood interface.
   *
   * Other methods are implemented in the AbstractHomogeneousTreeLikelihood class.
   *
   * @{
   */
  void setData(const SiteContainer& sites) throw (Exception);

public:
  // Specific methods:
  void initialize() throw (Exception);

  virtual void computeTreeDLikelihood(const string& variable);

  virtual void computeTreeD2Likelihood(const string& variable);

  /**
   * @brief return the probability of this objetct in the hierarchy
   *
   */
  
  double getProbability() const;
protected:
  /**
   * @brief Compute the likelihood for a subtree defined by the Tree::Node <i>node</i>.
   *
   * @param node The root of the subtree.
   */
  virtual void computeSubtreeLikelihood(const Node* node); // Recursive method.

  virtual void computeDownSubtreeDLikelihood(const Node*);

  virtual void computeDownSubtreeD2Likelihood(const Node*);

  void fireParameterChanged(const ParameterList& params);

  void computeTransitionProbabilitiesForNode(const Node* node);

};
} // end of namespace bpp.

#endif  // _RNONHOMOGENEOUSMIXEDTREELIKELIHOOD_H_

