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
  
      bool nullDLikelihood_;
      bool nullD2Likelihood_;


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
        const SubstitutionProcess* process,
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
        const SubstitutionProcess* process,
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
      void setData(const SiteContainer& sites);
  
      DoubleRecursiveTreeLikelihoodData* getLikelihoodData() { return likelihoodData_.get(); }
    
      const DoubleRecursiveTreeLikelihoodData* getLikelihoodData() const { return likelihoodData_.get(); }
  
      double getLikelihoodForASite(size_t site);

      double getLikelihoodForASiteForAState(size_t site, int state);

      double getLikelihoodForASiteForAClass(size_t site, size_t modelClass);

      double getLikelihoodForASiteForAClassForAState(size_t site, size_t modelClass, int state);
 
      double getDLikelihoodForASite(size_t site);

      double getD2LikelihoodForASite(size_t site);

      void computeTreeDLogLikelihood(const std::string& variable);

      void computeTreeD2LogLikelihood(const std::string& variable);

      /**
       *@brief Compute the conditional Likelihood Array at a given node.
       *
       *@param node the given node
       *@param likelihoodArray the class/site/state array storing the
       *       likelihoods
       *@param sonNode in case the computing is called by a son, the
       *       Node of this son (which partial likelihood is not
       *       included in the computation); default : 0 (which means
       *       all sons are considered).
       *
       */
    
      void computeLikelihoodAtNode(const Node* node, VVVdouble& likelihoodArray, const Node* sonNode =0 );

    protected:
      void computeTreeLikelihood();
    
    private:
  
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

      void computeTreeDLikelihoodAtNode(const Node* node);
    
      void computeTreeD2LikelihoodAtNode(const Node* node);

      void computeRootLikelihood();
  
      /**
       * @brief This method is mainly for debugging purpose.
       *
       * @param node The node at which likelihood values must be displayed.
       */
      virtual void displayLikelihood(const Node* node);

    public:
    
      /**
       * @brief Compute the posterior probabilities for each state and
       * each class of each distinct site.
       *
       * @param nodeId The id of the node at which probabilities must be
       * computed.
       * @return A 3-dimensional array, with probabilities for each
       * site, each rate and each state.
       */
    
      VVVdouble getPosteriorProbabilitiesForEachStateForEachClass(int nodeId);

      /**
       * @brief Compute the posterior probabilities for each state for a
       * given node.
       *
       * This method calls the
       * getPosteriorProbabilitiesForEachStateForEachClass function and
       * average the probabilities over all sites and classes,
       * resulting in a one-dimensionnal frequency array, with one
       * frequency per model state.
       *
       * @param nodeId The id of the node at which probabilities must be
       * computed.
       * @return vector of double with state frequencies for the given
       * node.
       */

      Vdouble getPosteriorStateFrequencies(int nodeId);
    };

  } //end of namespace newlik.
} //end of namespace bpp.

#endif  //_DOUBLERECURSIVETREELIKELIHOOD_H_

