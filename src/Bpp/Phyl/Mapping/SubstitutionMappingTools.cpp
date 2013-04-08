//
// File: SubstitutionMappingTools.cpp
// Created by: Julien Dutheil
// Created on: Wed Apr 5 13:04 2006
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#include "SubstitutionMappingTools.h"
#include "../Likelihood/DRTreeLikelihoodTools.h"
#include "../Likelihood/MarginalAncestralStateReconstruction.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/DataTable.h>

using namespace bpp;

// From the STL:
#include <iomanip>

using namespace std;

/******************************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeSubstitutionVectors(
  const DRTreeLikelihood& drtl,
  SubstitutionCount& substitutionCount,
  bool verbose) throw (Exception)
{
  //Preamble:
  if (!drtl.isInitialized()) throw Exception("SubstitutionMappingTools::computeSubstitutionVectors(). Likelihood object is not initialized.");
                                   
  //A few variables we'll need:
  
  const TreeTemplate<Node> tree(drtl.getTree());
  const SiteContainer*    sequences = drtl.getData();
  const DiscreteDistribution* rDist = drtl.getRateDistribution();
    
  size_t nbSites         = sequences->getNumberOfSites();
  size_t nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  size_t nbStates        = sequences->getAlphabet()->getSize();
  size_t nbClasses       = rDist->getNumberOfCategories();
  size_t nbTypes         = substitutionCount.getNumberOfSubstitutionTypes();
  vector<const Node*> nodes    = tree.getNodes();
  const vector<size_t>* rootPatternLinks
                               = &drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes         = nodes.size();
  
  // We create a new ProbabilisticSubstitutionMapping object:
  ProbabilisticSubstitutionMapping* substitutions = new ProbabilisticSubstitutionMapping(tree, &substitutionCount, nbSites);
                                   
  // Store likelihood for each rate for each site:
  VVVdouble lik;
  drtl.computeLikelihoodAtNode(tree.getRootId(), lik);
  Vdouble Lr(nbDistinctSites, 0);
  Vdouble rcProbs = rDist->getProbabilities();
  Vdouble rcRates = rDist->getCategories();
  for (size_t i = 0; i < nbDistinctSites; i++)
  {
    VVdouble* lik_i = &lik[i];
    for (size_t c = 0; c < nbClasses; c++)
    {
      Vdouble* lik_i_c = &(*lik_i)[c];
      double rc = rDist->getProbability(c);
      for (size_t s = 0; s < nbStates; s++)
      {
        Lr[i] += (*lik_i_c)[s] * rc;
      }
    }
  }

  // Compute the number of substitutions for each class and each branch in the tree:
  if (verbose) ApplicationTools::displayTask("Compute joint node-pairs likelihood", true);
  
  for (size_t l = 0; l < nbNodes; l++)
  {
    //For each node,
    const Node* currentNode = nodes[l];

    const Node* father = currentNode->getFather();

    double d = currentNode->getDistanceToFather();
 
    if (verbose) ApplicationTools::displayGauge(l, nbNodes-1, '>');
    VVdouble substitutionsForCurrentNode(nbDistinctSites);
    for (size_t i = 0; i < nbDistinctSites; ++i)
      substitutionsForCurrentNode[i].resize(nbTypes);

    // Now we've got to compute likelihoods in a smart manner... ;)
    VVVdouble likelihoodsFatherConstantPart(nbDistinctSites);
    for (size_t i = 0; i < nbDistinctSites; i++)
    {
      VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
      likelihoodsFatherConstantPart_i->resize(nbClasses);
      for (size_t c = 0; c < nbClasses; c++)
      {
        Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
        likelihoodsFatherConstantPart_i_c->resize(nbStates);
        double rc = rDist->getProbability(c);
        for (size_t s = 0; s < nbStates; s++)
        {
          //(* likelihoodsFatherConstantPart_i_c)[s] = rc * model->freq(s);
          //freq is already accounted in the array
          (* likelihoodsFatherConstantPart_i_c)[s] = rc;
        }
      }
    }
    
    // First, what will remain constant:
    size_t nbSons =  father->getNumberOfSons();
    for (size_t n = 0; n < nbSons; n++)
    {
      const Node* currentSon = father->getSon(n);
      if (currentSon->getId() != currentNode->getId())
      {
        const VVVdouble* likelihoodsFather_son = &drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentSon->getId());

        //Now iterate over all site partitions:
        auto_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentSon->getId()));
        VVVdouble pxy;
        bool first;
        while (mit->hasNext())
        {
          TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
          auto_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
          first = true;
          while (sit->hasNext())
          {
            size_t i = sit->next();
            //We retrieve the transition probabilities for this site partition:
            if (first)
            {
              pxy = drtl.getTransitionProbabilitiesPerRateClass(currentSon->getId(), i);
              first = false;
            }
            const VVdouble* likelihoodsFather_son_i = &(*likelihoodsFather_son)[i];
            VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
            for (size_t c = 0; c < nbClasses; c++)
            {
              const Vdouble* likelihoodsFather_son_i_c = &(*likelihoodsFather_son_i)[c];
              Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
              VVdouble* pxy_c = &pxy[c];
              for (size_t x = 0; x < nbStates; x++)
              {
                Vdouble* pxy_c_x = &(*pxy_c)[x];
                double likelihood = 0.;
                for (size_t y = 0; y < nbStates; y++)
                {
                  likelihood += (*pxy_c_x)[y] * (*likelihoodsFather_son_i_c)[y];
                }
                (*likelihoodsFatherConstantPart_i_c)[x] *= likelihood;
              }
            }
          }
        }      
      }
    }
    if (father->hasFather())
    {
      const Node* currentSon = father->getFather();
      const VVVdouble* likelihoodsFather_son = &drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentSon->getId());
      //Now iterate over all site partitions:
      auto_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(father->getId()));
      VVVdouble pxy;
      bool first;
      while (mit->hasNext())
      {
        TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
        auto_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
        first = true;
        while (sit->hasNext())
        {
          size_t i = sit->next();
          //We retrieve the transition probabilities for this site partition:
          if (first)
          {
            pxy = drtl.getTransitionProbabilitiesPerRateClass(father->getId(), i);
            first = false;
          }
          const VVdouble* likelihoodsFather_son_i = &(*likelihoodsFather_son)[i];
          VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
          for (size_t c = 0; c < nbClasses; c++)
          {
            const Vdouble* likelihoodsFather_son_i_c = &(*likelihoodsFather_son_i)[c];
            Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
            VVdouble* pxy_c = &pxy[c]; 
            for (size_t x = 0; x < nbStates; x++)
            {
              double likelihood = 0.;
              for (size_t y = 0; y < nbStates; y++)
              {
                Vdouble* pxy_c_x = &(*pxy_c)[y];
                likelihood += (*pxy_c_x)[x] * (*likelihoodsFather_son_i_c)[y];
              }
              (*likelihoodsFatherConstantPart_i_c)[x] *= likelihood;
            }
          }
        }
      }      
    }
    else
    {
      //Account for root frequencies:
      for (size_t i = 0; i < nbDistinctSites; i++)
      {
        vector<double> freqs = drtl.getRootFrequencies(i);
        VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
        for (size_t c = 0; c < nbClasses; c++)
        {
          Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
          for (size_t x = 0; x < nbStates; x++)
          {
            (*likelihoodsFatherConstantPart_i_c)[x] *= freqs[x];
          }
        }
      }      
    }

    // Then, we deal with the node of interest.
    // We first average uppon 'y' to save computations, and then uppon 'x'.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)

    //Iterate over all site partitions:
		const VVVdouble* likelihoodsFather_node = &(drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentNode->getId()));
    auto_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentNode->getId()));
    VVVdouble pxy;
    bool first;
    while (mit->hasNext())
    {
      TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
      substitutionCount.setSubstitutionModel(bmd->getModel());
      //compute all nxy first:
      VVVVdouble nxy(nbClasses);
      for (size_t c = 0; c < nbClasses; ++c)
      {
        VVVdouble* nxy_c = &nxy[c];
        double rc = rcRates[c];
        nxy_c->resize(nbTypes);
        for (size_t t = 0; t < nbTypes; ++t)
        {
          VVdouble* nxy_c_t = &(*nxy_c)[t];
          Matrix<double>* nijt = substitutionCount.getAllNumbersOfSubstitutions(d * rc, t + 1);
          nxy_c_t->resize(nbStates);
          for (size_t x = 0; x < nbStates; ++x)
          {
            Vdouble* nxy_c_t_x = &(*nxy_c_t)[x];
            nxy_c_t_x->resize(nbStates);
            for (size_t y = 0; y < nbStates; ++y)
              (*nxy_c_t_x)[y] = (*nijt)(x, y);
          }
          delete nijt;
        }
      }

      //Now loop over sites:
      auto_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
      first = true;
      while (sit->hasNext())
      {
        size_t i = sit->next();
        //We retrieve the transition probabilities and substitution counts for this site partition:
        if (first)
        {
          pxy = drtl.getTransitionProbabilitiesPerRateClass(currentNode->getId(), i);  
          first = false;
        }
        const VVdouble* likelihoodsFather_node_i = &(*likelihoodsFather_node)[i];
        VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
        for (size_t c = 0; c < nbClasses; ++c)
        {
          const Vdouble* likelihoodsFather_node_i_c = &(*likelihoodsFather_node_i)[c];
          Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
          const VVdouble* pxy_c = &pxy[c];
          VVVdouble* nxy_c = &nxy[c];
          for (size_t x = 0; x < nbStates; ++x)
          {
            double* likelihoodsFatherConstantPart_i_c_x = &(*likelihoodsFatherConstantPart_i_c)[x];
            const Vdouble* pxy_c_x = &(*pxy_c)[x];
            for (size_t y = 0; y < nbStates; ++y)
            {
              double likelihood_cxy = (*likelihoodsFatherConstantPart_i_c_x)
                 *(*pxy_c_x)[y]
                 *(*likelihoodsFather_node_i_c)[y];

              for (size_t t = 0; t < nbTypes; ++t) {
                // Now the vector computation:
                substitutionsForCurrentNode[i][t] += likelihood_cxy * (*nxy_c)[t][x][y];
                //                                   <------------>   <--------------->
                // Posterior probability                   |                 | 
                // for site i and rate class c *           |                 |
                // likelihood for this site----------------+                 |
                //                                                           |
                //Substitution function for site i and rate class c----------+
              }
            }          
          }
        }
      }
    }
    
    //Now we just have to copy the substitutions into the result vector:
    for (size_t i = 0; i < nbSites; ++i)
      for (size_t t = 0; t < nbTypes; ++t)
        (*substitutions)(l, i, t) = substitutionsForCurrentNode[(* rootPatternLinks)[i]][t] / Lr[(* rootPatternLinks)[i]];
  }
  if (verbose)
  {
    if (ApplicationTools::message) *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
  return substitutions;
}

/**************************************************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(
  const DRTreeLikelihood& drtl,
  SubstitutionCount& substitutionCount,
  bool verbose) throw (Exception)
{
  //Preamble:
  if (!drtl.isInitialized()) throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(). Likelihood object is not initialized.");
                                   
  //A few variables we'll need:
  const TreeTemplate<Node> tree(drtl.getTree());
  const SiteContainer*    sequences = drtl.getData();
  const DiscreteDistribution* rDist = drtl.getRateDistribution();
    
  size_t nbSites         = sequences->getNumberOfSites();
  size_t nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  size_t nbStates        = sequences->getAlphabet()->getSize();
  size_t nbClasses       = rDist->getNumberOfCategories();
  size_t nbTypes         = substitutionCount.getNumberOfSubstitutionTypes();
  vector<const Node *> nodes   = tree.getNodes();
  const vector<size_t> * rootPatternLinks
                               = &drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes = nodes.size();
  
  // We create a new ProbabilisticSubstitutionMapping object:
  ProbabilisticSubstitutionMapping* substitutions = new ProbabilisticSubstitutionMapping(tree, &substitutionCount, nbSites);
                                   
  Vdouble rcRates = rDist->getCategories();

  // Compute the number of substitutions for each class and each branch in the tree:
  if (verbose) ApplicationTools::displayTask("Compute joint node-pairs likelihood", true);
  
  for (size_t l = 0; l < nbNodes; ++l)
  {
    // For each node,
    const Node* currentNode = nodes[l];

    const Node* father = currentNode->getFather();

    double d = currentNode->getDistanceToFather();
    
    if (verbose) ApplicationTools::displayGauge(l, nbNodes-1, '>');
    VVdouble substitutionsForCurrentNode(nbDistinctSites);
    for (size_t i = 0; i < nbDistinctSites; ++i)
      substitutionsForCurrentNode[i].resize(nbTypes);

    // Now we've got to compute likelihoods in a smart manner... ;)
    VVVdouble likelihoodsFatherConstantPart(nbDistinctSites);
    for (size_t i = 0; i < nbDistinctSites; ++i)
    {
      VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
      likelihoodsFatherConstantPart_i->resize(nbClasses);
      for (size_t c = 0; c < nbClasses; ++c)
      {
        Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
        likelihoodsFatherConstantPart_i_c->resize(nbStates);
        double rc = rDist->getProbability(c);
        for (size_t s = 0; s < nbStates; ++s)
        {
          //(* likelihoodsFatherConstantPart_i_c)[s] = rc * model->freq(s);
          //freq is already accounted in the array
          (* likelihoodsFatherConstantPart_i_c)[s] = rc;
        }
      }
    }
    
    // First, what will remain constant:
    size_t nbSons =  father->getNumberOfSons();
    for (size_t n = 0; n < nbSons; ++n)
    {
      const Node* currentSon = father->getSon(n);    
      if (currentSon->getId() != currentNode->getId())
      {
        const VVVdouble* likelihoodsFather_son = &drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentSon->getId());

        //Now iterate over all site partitions:
        auto_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentSon->getId()));
        VVVdouble pxy;
        bool first;
        while (mit->hasNext())
        {
          TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
          auto_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
          first = true;
          while (sit->hasNext())
          {
            size_t i = sit->next();
            //We retrieve the transition probabilities for this site partition:
            if (first)
            {
              pxy = drtl.getTransitionProbabilitiesPerRateClass(currentSon->getId(), i);
              first = false;
            }
            const VVdouble* likelihoodsFather_son_i = &(*likelihoodsFather_son)[i];
            VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
            for (size_t c = 0; c < nbClasses; ++c)
            {
              const Vdouble* likelihoodsFather_son_i_c = &(*likelihoodsFather_son_i)[c];
              Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
              VVdouble* pxy_c = & pxy[c]; 
              for (size_t x = 0; x < nbStates; ++x)
              {
                Vdouble* pxy_c_x = &(*pxy_c)[x];
                double likelihood = 0.;
                for (size_t y = 0; y < nbStates; ++y)
                {
                  likelihood += (*pxy_c_x)[y] * (*likelihoodsFather_son_i_c)[y];
                }
                (*likelihoodsFatherConstantPart_i_c)[x] *= likelihood;
              }
            }
          }
        }      
      }
    }
    if (father->hasFather())
    {
      const Node* currentSon = father->getFather();
      const VVVdouble* likelihoodsFather_son = &drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentSon->getId());
      //Now iterate over all site partitions:
      auto_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(father->getId()));
      VVVdouble pxy;
      bool first;
      while (mit->hasNext())
      {
        TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
        auto_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
        first = true;
        while (sit->hasNext())
        {
          size_t i = sit->next();
          //We retrieve the transition probabilities for this site partition:
          if (first)
          {
            pxy = drtl.getTransitionProbabilitiesPerRateClass(father->getId(), i);
            first = false;
          }
          const VVdouble* likelihoodsFather_son_i = &(*likelihoodsFather_son)[i];
          VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
          for (size_t c = 0; c < nbClasses; ++c)
          {
            const Vdouble* likelihoodsFather_son_i_c = &(*likelihoodsFather_son_i)[c];
            Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
            VVdouble* pxy_c = &pxy[c]; 
            for (size_t x = 0; x < nbStates; ++x)
            {
              double likelihood = 0.;
              for (size_t y = 0; y < nbStates; ++y)
              {
                Vdouble* pxy_c_x = &(*pxy_c)[y];
                likelihood += (*pxy_c_x)[x] * (*likelihoodsFather_son_i_c)[y];
              }
              (*likelihoodsFatherConstantPart_i_c)[x] *= likelihood;
            }
          }
        }
      }      
    }
    else
    {
      //Account for root frequencies:
      for (size_t i = 0; i < nbDistinctSites; ++i)
      {
        vector<double> freqs = drtl.getRootFrequencies(i);
        VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
        for (size_t c = 0; c < nbClasses; ++c)
        {
          Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
          for (size_t x = 0; x < nbStates; ++x)
          {
            (*likelihoodsFatherConstantPart_i_c)[x] *= freqs[x]; 
          }
        }
      }      
    }

    // Then, we deal with the node of interest.
    // We first average uppon 'y' to save computations, and then uppon 'x'.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)

    //Iterate over all site partitions:
		const VVVdouble* likelihoodsFather_node = &drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentNode->getId());
    auto_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentNode->getId()));
    VVVdouble pxy;
    bool first;
    while (mit->hasNext())
    {
      TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
      substitutionCount.setSubstitutionModel(bmd->getModel());
      //compute all nxy first:
      VVVVdouble nxy(nbClasses);
      for (size_t c = 0; c < nbClasses; ++c)
      {
        double rc = rcRates[c];
        VVVdouble* nxy_c = &nxy[c];
        nxy_c->resize(nbTypes);
        for (size_t t = 0; t < nbTypes; ++t)
        {
          VVdouble* nxy_c_t = &(*nxy_c)[t];
          nxy_c_t->resize(nbStates);
          Matrix<double>* nijt = substitutionCount.getAllNumbersOfSubstitutions(d * rc, t + 1);
          for (size_t x = 0; x < nbStates; ++x)
          {
            Vdouble* nxy_c_t_x = &(*nxy_c_t)[x];
            nxy_c_t_x->resize(nbStates);
            for (size_t y = 0; y < nbStates; ++y)
            {
              (*nxy_c_t_x)[y] = (*nijt)(x, y);
            }
          }
          delete nijt;
        }
      }

      //Now loop over sites:
      auto_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
      first = true;
      while (sit->hasNext())
      {
        size_t i = sit->next();
        //We retrieve the transition probabilities and substitution counts for this site partition:
        if (first)
        {
          pxy = drtl.getTransitionProbabilitiesPerRateClass(currentNode->getId(), i);  
          first = false;
        }
        const VVdouble* likelihoodsFather_node_i = &(*likelihoodsFather_node)[i];
        VVdouble* likelihoodsFatherConstantPart_i = &likelihoodsFatherConstantPart[i];
        RowMatrix<double> pairProbabilities(nbStates, nbStates);
        MatrixTools::fill(pairProbabilities, 0.);
        VVVdouble subsCounts(nbStates);
        for (size_t j = 0; j < nbStates; ++j) {
          subsCounts[j].resize(nbStates);
          for (size_t k = 0; k < nbStates; ++k) {
            subsCounts[j][k].resize(nbTypes);
          }
        }
        for (size_t c = 0; c < nbClasses; ++c)
        {
          const Vdouble* likelihoodsFather_node_i_c = &(*likelihoodsFather_node_i)[c];
          Vdouble* likelihoodsFatherConstantPart_i_c = &(*likelihoodsFatherConstantPart_i)[c];
          const VVdouble* pxy_c = &pxy[c];
          VVVdouble* nxy_c = &nxy[c];
          for (size_t x = 0; x < nbStates; ++x)
          {
            double* likelihoodsFatherConstantPart_i_c_x = &(*likelihoodsFatherConstantPart_i_c)[x];
            const Vdouble* pxy_c_x = &(*pxy_c)[x];
            for (size_t y = 0; y < nbStates; ++y)
            {
              double likelihood_cxy = (*likelihoodsFatherConstantPart_i_c_x)
                 *(*pxy_c_x)[y]
                 *(*likelihoodsFather_node_i_c)[y];
              pairProbabilities(x, y) += likelihood_cxy; // Sum over all rate classes.
              for (size_t t = 0; t < nbTypes; ++t) {
                subsCounts[x][y][t] += likelihood_cxy * (* nxy_c)[t][x][y];
              }
            }
          }
        }
        // Now the vector computation:
        // Here we do not average over all possible pair of ancestral states,
        // We only consider the one with max likelihood:
        vector<size_t> xy = MatrixTools::whichMax(pairProbabilities);
        for (size_t t = 0; t < nbTypes; ++t) {
          substitutionsForCurrentNode[i][t] += subsCounts[xy[0]][xy[1]][t] / pairProbabilities(xy[0], xy[1]);
        }
      }
    }
    //Now we just have to copy the substitutions into the result vector:
    for (size_t i = 0; i < nbSites; i++)
      for (size_t t = 0; t < nbTypes; t++)
        (*substitutions)(l, i, t) = substitutionsForCurrentNode[(*rootPatternLinks)[i]][t];
  }
  if (verbose)
  {
    if (ApplicationTools::message) *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
  return substitutions;
}

/**************************************************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(
  const DRTreeLikelihood& drtl,
  SubstitutionCount& substitutionCount,
  bool verbose) throw (Exception)
{
  //Preamble:
  if (!drtl.isInitialized()) throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(). Likelihood object is not initialized.");
                                   
  //A few variables we'll need:
  
  const TreeTemplate<Node> tree(drtl.getTree());
  const SiteContainer*    sequences = drtl.getData();
  const DiscreteDistribution* rDist = drtl.getRateDistribution();
  const Alphabet*             alpha = sequences->getAlphabet();
    
  size_t nbSites         = sequences->getNumberOfSites();
  size_t nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  size_t nbStates        = alpha->getSize();
  size_t nbTypes         = substitutionCount.getNumberOfSubstitutionTypes();
  vector<const Node*> nodes    = tree.getNodes();
  const vector<size_t>* rootPatternLinks
                               = &drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes = nodes.size();
  
  // We create a new ProbabilisticSubstitutionMapping object:
  ProbabilisticSubstitutionMapping* substitutions = new ProbabilisticSubstitutionMapping(tree, &substitutionCount, nbSites);
  
  // Compute the whole likelihood of the tree according to the specified model:
  
  Vdouble rcRates = rDist->getCategories();

  // Compute the number of substitutions for each class and each branch in the tree:
  if (verbose) ApplicationTools::displayTask("Compute marginal ancestral states");
  MarginalAncestralStateReconstruction masr(&drtl);
  map<int, vector<size_t> > ancestors = masr.getAllAncestralStates();
  if (verbose) ApplicationTools::displayTaskDone();

  // Now we just have to compute the substitution vectors:
  if (verbose) ApplicationTools::displayTask("Compute substitution vectors", true);
  
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* currentNode = nodes[l];
    
    const Node* father = currentNode->getFather();

    double d = currentNode->getDistanceToFather();

    vector<size_t> nodeStates = ancestors[currentNode->getId()]; //These are not 'true' ancestors ;)
    vector<size_t> fatherStates = ancestors[father->getId()];
    
    //For each node,
    if (verbose) ApplicationTools::displayGauge(l, nbNodes-1, '>');
    VVdouble substitutionsForCurrentNode(nbDistinctSites);
    for (size_t i = 0; i < nbDistinctSites; ++i)
      substitutionsForCurrentNode[i].resize(nbTypes);
    
    // Here, we have no likelihood computation to do!

    // Then, we deal with the node of interest.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)
    // Iterate over all site partitions:
    auto_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentNode->getId()));
    while (mit->hasNext())
    {
      TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
      substitutionCount.setSubstitutionModel(bmd->getModel());
      //compute all nxy first:
      VVVdouble nxyt(nbTypes);
      for (size_t t = 0; t < nbTypes; ++t)
      {
        nxyt[t].resize(nbStates);
        Matrix<double>* nxy = substitutionCount.getAllNumbersOfSubstitutions(d, t + 1);
        for (size_t x = 0; x < nbStates; ++x) {
          nxyt[t][x].resize(nbStates);
          for (size_t y = 0; y < nbStates; ++y) {
            nxyt[t][x][y] = (*nxy)(x, y);
          }
        }
        delete nxy;
      }
      //Now loop over sites:
      auto_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
      while (sit->hasNext())
      {
        size_t i = sit->next();
        size_t fatherState = fatherStates[i];
        size_t nodeState   = nodeStates[i];
        if (fatherState >= nbStates || nodeState >= nbStates)
          for (size_t t = 0; t < nbTypes; ++t)
            substitutionsForCurrentNode[i][t] = 0; // To be conservative! Only in case there are generic characters.
        else
          for (size_t t = 0; t < nbTypes; ++t)
            substitutionsForCurrentNode[i][t] = nxyt[t][fatherState][nodeState];
      }
    }
    
    //Now we just have to copy the substitutions into the result vector:
    for (size_t i = 0; i < nbSites; i++)
      for (size_t t = 0; t < nbTypes; t++)
        (*substitutions)(l, i, t) = substitutionsForCurrentNode[(*rootPatternLinks)[i]][t];
  }
  if (verbose)
  {
    if (ApplicationTools::message) *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
  return substitutions;
}

/**************************************************************************************************/

ProbabilisticSubstitutionMapping* SubstitutionMappingTools::computeSubstitutionVectorsMarginal(
  const DRTreeLikelihood& drtl,
  SubstitutionCount& substitutionCount,
  bool verbose) throw (Exception)
{
  //Preamble:
  if (!drtl.isInitialized()) throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsMarginal(). Likelihood object is not initialized.");
                                   
  //A few variables we'll need:
  
  const TreeTemplate<Node>tree(drtl.getTree());
  const SiteContainer*    sequences = drtl.getData();
  const DiscreteDistribution* rDist = drtl.getRateDistribution();
    
  size_t nbSites         = sequences->getNumberOfSites();
  size_t nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  size_t nbStates        = sequences->getAlphabet()->getSize();
  size_t nbClasses       = rDist->getNumberOfCategories();
  size_t nbTypes         = substitutionCount.getNumberOfSubstitutionTypes();
  vector<const Node*> nodes    = tree.getNodes();
  const vector<size_t>* rootPatternLinks
                               = &drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  size_t nbNodes = nodes.size();
  
  // We create a new ProbabilisticSubstitutionMapping object:
  ProbabilisticSubstitutionMapping* substitutions = new ProbabilisticSubstitutionMapping(tree, &substitutionCount, nbSites);
                                     
  // Compute the whole likelihood of the tree according to the specified model:
  
  Vdouble rcProbs = rDist->getProbabilities();
  Vdouble rcRates = rDist->getCategories();

  //II) Compute the number of substitutions for each class and each branch in the tree:
  if (verbose) ApplicationTools::displayTask("Compute marginal node-pairs likelihoods", true);
  
  for (size_t l = 0; l < nbNodes; l++)
  {
    const Node* currentNode = nodes[l];
    
    const Node* father = currentNode->getFather();

    double d = currentNode->getDistanceToFather();
    
    //For each node,
    if (verbose) ApplicationTools::displayGauge(l, nbNodes-1, '>');
    VVdouble substitutionsForCurrentNode(nbDistinctSites);
    for (size_t i = 0; i < nbDistinctSites; ++i)
      substitutionsForCurrentNode[i].resize(nbTypes);

    // Then, we deal with the node of interest.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)
    VVVdouble probsNode   = DRTreeLikelihoodTools::getPosteriorProbabilitiesForEachStateForEachRate(drtl, currentNode->getId());
    VVVdouble probsFather = DRTreeLikelihoodTools::getPosteriorProbabilitiesForEachStateForEachRate(drtl, father->getId());

    //Iterate over all site partitions:
    auto_ptr<TreeLikelihood::ConstBranchModelIterator> mit(drtl.getNewBranchModelIterator(currentNode->getId()));
    while (mit->hasNext())
    {
      TreeLikelihood::ConstBranchModelDescription* bmd = mit->next();
      substitutionCount.setSubstitutionModel(bmd->getModel());
      //compute all nxy first:
      VVVVdouble nxy(nbClasses);
      for (size_t c = 0; c < nbClasses; ++c)
      {
        VVVdouble* nxy_c = &nxy[c];
        double rc = rcRates[c];
        nxy_c->resize(nbTypes);
        for (size_t t = 0; t < nbTypes; ++t)
        {
          VVdouble* nxy_c_t = &(*nxy_c)[t];
          Matrix<double>* nijt = substitutionCount.getAllNumbersOfSubstitutions(d * rc, t + 1);
          nxy_c_t->resize(nbStates);
          for (size_t x = 0; x < nbStates; ++x)
          {
            Vdouble* nxy_c_t_x = &(*nxy_c_t)[x];
            nxy_c_t_x->resize(nbStates);
            for (size_t y = 0; y < nbStates; ++y)
            {
              (*nxy_c_t_x)[y] = (*nijt)(x, y);
            }
          }
          delete nijt;
        }
      }

      //Now loop over sites:
      auto_ptr<TreeLikelihood::SiteIterator> sit(bmd->getNewSiteIterator());
      while (sit->hasNext())
      {
        size_t i = sit->next();
        VVdouble* probsNode_i   = & probsNode[i];
        VVdouble* probsFather_i = & probsFather[i];
        for (size_t c = 0; c < nbClasses; ++c)
        {
          Vdouble* probsNode_i_c   = &(*probsNode_i)[c];
          Vdouble* probsFather_i_c = &(*probsFather_i)[c];
          VVVdouble* nxy_c = &nxy[c];
          for (size_t x = 0; x < nbStates; ++x)
          {
            for (size_t y = 0; y < nbStates; ++y)
            {
              double prob_cxy = (*probsFather_i_c)[x] * (*probsNode_i_c)[y];
              // Now the vector computation:
              for (size_t t = 0; t < nbTypes; ++t) {
                substitutionsForCurrentNode[i][t] += prob_cxy * (*nxy_c)[t][x][y];
                //                                   <------>   <--------------->
                // Posterior probability                 |                | 
                // for site i and rate class c *         |                |
                // likelihood for this site--------------+                |
                //                                                        |
                //Substitution function for site i and rate class c-------+
              }
            }
          }
        }
      }
    }
    
    //Now we just have to copy the substitutions into the result vector:
    for (size_t i = 0; i < nbSites; ++i)
      for (size_t t = 0; t < nbTypes; ++t)
        (*substitutions)(l, i, t) = substitutionsForCurrentNode[(* rootPatternLinks)[i]][t];
  }
  if (verbose)
  {
    if (ApplicationTools::message) *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
  return substitutions;
}

/**************************************************************************************************/

void SubstitutionMappingTools::writeToStream(
  const ProbabilisticSubstitutionMapping& substitutions,
  const SiteContainer& sites,
  size_t type,
  ostream& out)
  throw (IOException) 
{
  if (!out) throw IOException("SubstitutionMappingTools::writeToFile. Can't write to stream.");
  out << "Branches";
  out << "\tMean";
  for (size_t i = 0; i < substitutions.getNumberOfSites(); i++)
  {
    out << "\tSite" << sites.getSite(i).getPosition();
  }
  out << endl;
  
  for (size_t j = 0; j < substitutions.getNumberOfBranches(); j++)
  {
    out << substitutions.getNode(j)->getId() << "\t" << substitutions.getNode(j)->getDistanceToFather();
    for (size_t i = 0; i < substitutions.getNumberOfSites(); i++)
    {
      out << "\t" << substitutions(j, i, type);
    }
    out << endl;
  }
}

/**************************************************************************************************/

void SubstitutionMappingTools::readFromStream(istream& in, ProbabilisticSubstitutionMapping& substitutions, size_t type)
  throw (IOException)
{
  try
  {
    DataTable* data = DataTable::read(in, "\t", true, -1);
    vector<string> ids = data->getColumn(0);
    data->deleteColumn(0);//Remove ids
    data->deleteColumn(0);//Remove means
    //Now parse the table:
    size_t nbSites = data->getNumberOfColumns();
    substitutions.setNumberOfSites(nbSites);
    size_t nbBranches = data->getNumberOfRows();
    for (size_t i = 0; i < nbBranches; i++)
    {
      int id = TextTools::toInt(ids[i]);
      size_t br = substitutions.getNodeIndex(id);
      for (size_t j = 0; j < nbSites; j++)
      {
        substitutions(br, j, type) = TextTools::toDouble((*data)(i, j));
      }
    }
    //Parse the header:
    for (size_t i = 0; i < nbSites; i++)
    {
      string siteTxt = data->getColumnName(i);
      int site = 0;
      if (siteTxt.substr(0,4) == "Site") site = TextTools::to<int>(siteTxt.substr(4));
      else site = TextTools::to<int>(siteTxt);
      substitutions.setSitePosition(i, site);
    }
    
    delete data;
  }
  catch (Exception& e)
  {
    throw IOException(string("Bad input file. ") + e.what());
  }
}

/**************************************************************************************************/
    
vector<double> SubstitutionMappingTools::computeTotalSubstitutionVectorForSite(const SubstitutionMapping& smap, size_t siteIndex)
{
  size_t nbBranches = smap.getNumberOfBranches();
  size_t nbTypes    = smap.getNumberOfSubstitutionTypes();
  Vdouble v(nbBranches);
  for (size_t l = 0; l < nbBranches; ++l) {
    v[l] = 0;
    for (size_t t = 0; t < nbTypes; ++t) {
      v[l] += smap(l, siteIndex, t);
    }
  }
  return v;
}

/**************************************************************************************************/

double SubstitutionMappingTools::computeNormForSite(const SubstitutionMapping& smap, size_t siteIndex)
{
  double sumSquare = 0;
  for (size_t l = 0; l < smap.getNumberOfBranches(); ++l) {
    double sum = 0;
    for (size_t t = 0; t < smap.getNumberOfSubstitutionTypes(); ++t) {
      sum += smap(l, siteIndex, t);
    }
    sumSquare += sum * sum;
  }
  return sqrt(sumSquare);
}

/**************************************************************************************************/
    
vector<double> SubstitutionMappingTools::computeSumForBranch(const SubstitutionMapping& smap, size_t branchIndex)
{
  size_t nbSites = smap.getNumberOfSites();
  size_t nbTypes = smap.getNumberOfSubstitutionTypes();
  Vdouble v(nbTypes, 0);
  for (size_t i = 0; i < nbSites; ++i) {
    for (size_t t = 0; t < nbTypes; ++t) {
      v[t] += smap(branchIndex, i, t);
    }
  }
  return v;
}

/**************************************************************************************************/

vector<double> SubstitutionMappingTools::computeSumForSite(const SubstitutionMapping& smap, size_t siteIndex)
{
  size_t nbBranches = smap.getNumberOfBranches();
  size_t nbTypes = smap.getNumberOfSubstitutionTypes();
  Vdouble v(nbTypes, 0);
  for (size_t i = 0; i < nbBranches; ++i) {
    for (size_t t = 0; t < nbTypes; ++t) {
      v[t] += smap(i, siteIndex, t);
    }
  }
  return v;
}

/**************************************************************************************************/

