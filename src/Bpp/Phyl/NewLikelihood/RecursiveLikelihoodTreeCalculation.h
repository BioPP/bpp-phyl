//
// File: RecursiveLikelihoodTreeCalculation.h
// Created by: Julien Dutheil, Laurent Guéguen
// Created on: Tue May 15 14:30 2012
// From file: RNonHomogeneousLikelihoodTree.h
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

#ifndef _RECURSIVE_LIKELIHOOD_TREE_CALCULATION_H_
#define _RECURSIVE_LIKELIHOOD_TREE_CALCULATION_H_

#include "AbstractLikelihoodTreeCalculation.h"
#include "RecursiveLikelihoodTree.h"

#include <Bpp/Numeric/VectorTools.h>
//#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp
{
/**
 * @brief This class implements the single recursion likelihood computation for a tree.
 *
 * This class uses an instance of the SingleRecursiveLikelihoodTree for conditionnal likelihood storage.
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
 * initLikelihoodTree one, where all likelihoods for a given site <i>i</i> are at the <i>i</i> coordinate
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
    class RecursiveLikelihoodTreeCalculation:
      public AbstractLikelihoodTreeCalculation
    {
    private:
      mutable std::auto_ptr<RecursiveLikelihoodTree> likelihoodData_;
      int root1_, root2_; // Needed only in case of reparametrization of branch length at root node.
      // TODO: have to be initialized properly! We do not care of that for now. jdutheil on 11/12/12.

    public:
      /**
       * @brief Build a new Recursive Tree Likelihood object without data.
       *
       * This constructor only initialize the parameters.
       * To compute a likelihood, you will need to call the setData() and the computeLikelihoodTree() methods.
       *
       * @param process The substitution process to use.
       * @param verbose Should I display some info?
       * @param usePatterns Tell if recursive site compression should be performed.
       * @throw Exception in an error occured.
       */
      RecursiveLikelihoodTreeCalculation(
        const SubstitutionProcess* process,
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
      RecursiveLikelihoodTreeCalculation(
        const SiteContainer& data,
        const SubstitutionProcess* process,
        bool verbose = true,
        bool usePatterns = true)
        throw (Exception);

      /**
       * @brief Copy constructor.
       */ 
      RecursiveLikelihoodTreeCalculation(const RecursiveLikelihoodTreeCalculation& lik);

      RecursiveLikelihoodTreeCalculation& operator=(const RecursiveLikelihoodTreeCalculation& lik);

      virtual ~RecursiveLikelihoodTreeCalculation() {} // smart pointers take care of everything.

      RecursiveLikelihoodTreeCalculation* clone() const { return new RecursiveLikelihoodTreeCalculation(*this); }

    private:
      /**
       * @brief Method called by constructors.
       */
      void init_(bool usePatterns) throw (Exception);

    public:

      /**
       * @brief Methods to retrieve full likelihood.
       * !!! These methods do not check computation is up to date.
       * Call ComputeTreeDXLikelihood for this.
       *
       */
      
      AbstractLikelihoodTree& getLikelihoodData() { return *likelihoodData_.get(); }

      const AbstractLikelihoodTree& getLikelihoodData() const { return *likelihoodData_.get(); }

    protected:
      
      /**
       * @brief check the ComputingNodes to update recursively the
       * Likelihoods flags. If at least one ComputingNode has changed
       * all Above Likelihood flags of the trees are set to false.
       * Below DXLikelihood flags of the changed nodes are set to
       * false, and all their ancestors.
       *
       *
       * Set up2date_ to false if at least one ComputingNode has
       * changed, otherwise does not change up2date_.
       *
       */
      

      void updateLikelihoodFlags_();
      

    public:
      
      void updateLikelihood()
      {
        updateLikelihoodFlags_();
      }
      
      /**
       * @brief Compute derivatives
       *
       */
       
      void computeTreeDLogLikelihood(const std::string& variable);

      void computeTreeD2LogLikelihood(const std::string& variable);

      /*
       * @brief compute likelihood and set up2date_ to true.
       *
       */
      
      void computeTreeLikelihood();

      /**
       *@brief Compute the Likelihood at a given node, using ONLY
       *     downward recursion from the root, so upward recursion
       *     should be up to date.
       *
       *@param node the given node
       *
       */
    
      void computeLikelihoodsAtNode(int nodeId);

    public:

    };
} // end of namespace bpp.

#endif  // _RECURSIVE_LIKELIHOOD_TREE_CALCULATION_H_

