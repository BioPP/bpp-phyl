//
// File: SubstitutionMappingTools.cpp
// Created by: Julien Dutheil
// Created on: Wed Apr 5 13:04 2006
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004, 2005, 2006)

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
#include "DRTreeLikelihoodTools.h"
#include "MarginalAncestralStateReconstruction.h"

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>

// From NumCalc:
#include <NumCalc/MatrixTools.h>
#include <NumCalc/DataTable.h>

using namespace bpp;

// From the STL:
#include <iomanip>

using namespace std;

/******************************************************************************/

ProbabilisticSubstitutionMapping * SubstitutionMappingTools::computeSubstitutionVectors(
  const DRTreeLikelihood & drtl,
  SubstitutionCount & substitutionCount,
  bool verbose) throw (Exception)
{
  //Preamble:
  if(!drtl.isInitialized()) throw Exception("SubstitutionMappingTools::computeSubstitutionVectors(). Likelihood object is not initialized.");
                                   
  //A few variables we'll need:
  
  const TreeTemplate<Node> *    tree = dynamic_cast<const TreeTemplate<Node> *>(drtl.getTree());
  const SiteContainer *    sequences = drtl.getData();
  const DiscreteDistribution * rDist = drtl.getRateDistribution();
    
  unsigned int nbSites         = sequences->getNumberOfSites();
  unsigned int nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  unsigned int nbStates        = sequences->getAlphabet()->getSize();
  unsigned int nbClasses       = rDist->getNumberOfCategories();
  vector<const Node *> nodes   = tree->getNodes();
  const vector<unsigned int> * rootPatternLinks
                               = & drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  unsigned int nbNodes         = nodes.size();
  
  // We create a new ProbabilisticSubstitutionMapping object:
  ProbabilisticSubstitutionMapping * substitutions = new ProbabilisticSubstitutionMapping(*tree, nbSites);
                                   
  // Store likelihood for each rate for each site:
  VVVdouble l;
  drtl.computeLikelihoodAtNode(tree->getRootId(), l);
  Vdouble Lr(nbDistinctSites, 0);
  Vdouble rcProbs = rDist->getProbabilities();
  Vdouble rcRates = rDist->getCategories();
  for(unsigned int i = 0; i < nbDistinctSites; i++)
  {
    VVdouble * l_i = & l[i];
    for(unsigned int c = 0; c < nbClasses; c++)
    {
      Vdouble * l_i_c = & (* l_i)[c];
      double rc = rDist->getProbability(c);
      for(unsigned int s = 0; s < nbStates; s++)
      {
        Lr[i] += (* l_i_c)[s] * rc;
      }
    }
  }

  // Compute the number of substitutions for each class and each branch in the tree:
  if(verbose) ApplicationTools::displayTask("Compute joint node-pairs likelihood", true);
  
  for(unsigned int l = 0; l < nbNodes; l++)
  {
    const Node * currentNode = nodes[l];

    const Node * father = currentNode->getFather();

    double d = currentNode->getDistanceToFather();

    substitutionCount.setSubstitutionModel(drtl.getSubstitutionModelForNode(currentNode->getId()));
    
    //For each node,
    if(verbose) ApplicationTools::displayGauge(l, nbNodes-1, '>');
    Vdouble substitutionsForCurrentNode(nbDistinctSites, 0);

    //compute all nxy first:
    VVVdouble nxy = VVVdouble(nbClasses);
    for(unsigned int c = 0; c < nbClasses; c++)
    {
      VVdouble * nxy_c = & nxy[c];
      double rc = rcRates[c];
      Matrix<double> * nijt = substitutionCount.getAllNumbersOfSubstitutions(d * rc);
      nxy_c->resize(nbStates);
      for(unsigned int x = 0; x < nbStates; x++)
      {
        Vdouble * nxy_c_x = & (* nxy_c)[x];
        nxy_c_x->resize(nbStates);
        for(unsigned int y = 0; y < nbStates; y++)
        {
          (* nxy_c_x)[y] = (* nijt)(x, y);
        }
      }
      delete nijt;
    }
    
    // Now we've got to compute likelihoods in a smart manner... ;)

    VVVdouble likelihoodsFatherConstantPart(nbDistinctSites);
    for(unsigned int i = 0; i < nbDistinctSites; i++)
    {
      VVdouble * likelihoodsFatherConstantPart_i = & likelihoodsFatherConstantPart[i];
      likelihoodsFatherConstantPart_i->resize(nbClasses);
      for(unsigned int c = 0; c < nbClasses; c++)
      {
        Vdouble * likelihoodsFatherConstantPart_i_c = & (* likelihoodsFatherConstantPart_i)[c];
        likelihoodsFatherConstantPart_i_c->resize(nbStates);
        double rc = rDist->getProbability(c);
        for(unsigned int s = 0; s < nbStates; s++)
        {
          //(* likelihoodsFatherConstantPart_i_c)[s] = rc * model->freq(s);
          //freq is already accounted in the array
          (* likelihoodsFatherConstantPart_i_c)[s] = rc;
        }
      }
    }
    
    // First, what will remain constant:
    unsigned int nbSons =  father->getNumberOfSons();
    for(unsigned int n = 0; n < nbSons; n++)
    {
      const Node * currentSon = father->getSon(n);
      
      if(currentSon->getId() != currentNode->getId())
      {
        VVVdouble pxy = drtl.getTransitionProbabilitiesForNode(currentSon);
        const VVVdouble * likelihoodsFather_son = & drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentSon->getId());
        for(unsigned int i = 0; i < nbDistinctSites; i++)
        {
          const VVdouble * likelihoodsFather_son_i = & (* likelihoodsFather_son)[i];
          VVdouble * likelihoodsFatherConstantPart_i = & likelihoodsFatherConstantPart[i];
          for(unsigned int c = 0; c < nbClasses; c++)
          {
            const Vdouble * likelihoodsFather_son_i_c = & (* likelihoodsFather_son_i)[c];
            Vdouble * likelihoodsFatherConstantPart_i_c = & (* likelihoodsFatherConstantPart_i)[c];
            VVdouble * pxy_c = & pxy[c]; 
            for(unsigned int x = 0; x < nbStates; x++)
            {
              Vdouble * pxy_c_x = & (* pxy_c)[x];
              double likelihood = 0.;
              for(unsigned int y = 0; y < nbStates; y++)
              {
                likelihood += (* pxy_c_x)[y] * (* likelihoodsFather_son_i_c)[y];
              }
              (* likelihoodsFatherConstantPart_i_c)[x] *= likelihood;
            }
          }
        }      
      }
    }
    if(father->hasFather())
    {
      const Node * currentSon = father->getFather();
      VVVdouble pxy = drtl.getTransitionProbabilitiesForNode(father);
      const VVVdouble * likelihoodsFather_son = & drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentSon->getId());
      for(unsigned int i = 0; i < nbDistinctSites; i++)
      {
        const VVdouble * likelihoodsFather_son_i = & (* likelihoodsFather_son)[i];
        VVdouble * likelihoodsFatherConstantPart_i = & likelihoodsFatherConstantPart[i];
        for(unsigned int c = 0; c < nbClasses; c++)
        {
          const Vdouble * likelihoodsFather_son_i_c = & (* likelihoodsFather_son_i)[c];
          Vdouble * likelihoodsFatherConstantPart_i_c = & (* likelihoodsFatherConstantPart_i)[c];
          VVdouble * pxy_c = & pxy[c]; 
          for(unsigned int x = 0; x < nbStates; x++)
          {
            double likelihood = 0.;
            for(unsigned int y = 0; y < nbStates; y++)
            {
              Vdouble * pxy_c_x = & (* pxy_c)[y];
              likelihood += (* pxy_c_x)[x] * (* likelihoodsFather_son_i_c)[y];
            }
            (* likelihoodsFatherConstantPart_i_c)[x] *= likelihood;
          }
        }
      }      
    }
    else
    {
      //Account for root frequencies:
      vector<double> freqs = drtl.getRootFrequencies();
      for(unsigned int i = 0; i < nbDistinctSites; i++)
      {
        VVdouble * likelihoodsFatherConstantPart_i = & likelihoodsFatherConstantPart[i];
        for(unsigned int c = 0; c < nbClasses; c++)
        {
          Vdouble * likelihoodsFatherConstantPart_i_c = & (* likelihoodsFatherConstantPart_i)[c];
          for(unsigned int x = 0; x < nbStates; x++)
          {
            (* likelihoodsFatherConstantPart_i_c)[x] *= freqs[x]; 
          }
        }
      }      
    }

    // Then, we deal with the node of interest.
    // We first average uppon 'y' to save computations, and then uppon 'x'.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)
    const VVVdouble * pxy = & drtl.getTransitionProbabilitiesForNode(currentNode);
		const VVVdouble * likelihoodsFather_node = & drtl.getLikelihoodData()->getLikelihoodArray(father->getId(), currentNode->getId());
    for(unsigned int i = 0; i < nbDistinctSites; i++)
    {
      const VVdouble * likelihoodsFather_node_i = & (* likelihoodsFather_node)[i];
      VVdouble * likelihoodsFatherConstantPart_i = & likelihoodsFatherConstantPart[i];
      for(unsigned int c = 0; c < nbClasses; c++)
      {
        const Vdouble * likelihoodsFather_node_i_c = & (* likelihoodsFather_node_i)[c];
        Vdouble * likelihoodsFatherConstantPart_i_c = & (* likelihoodsFatherConstantPart_i)[c];
        const VVdouble * pxy_c = & (* pxy)[c];
        VVdouble * nxy_c = & nxy[c];
        for(unsigned int x = 0; x < nbStates; x++)
        {
          double * likelihoodsFatherConstantPart_i_c_x = & (* likelihoodsFatherConstantPart_i_c)[x];
          const Vdouble * pxy_c_x = & (* pxy_c)[x];
          Vdouble * nxy_c_x = & (* nxy_c)[x];
          for(unsigned int y = 0; y < nbStates; y++)
          {
            double likelihood_cxy = (* likelihoodsFatherConstantPart_i_c_x)
              * (* pxy_c_x)[y]
              * (* likelihoodsFather_node_i_c)[y];
            // Now the vector computation:
            substitutionsForCurrentNode[i] += likelihood_cxy * (* nxy_c_x)[y];
            //                                <------------>   <------------>
            // Posterior probability                |                 | 
            // for site i and rate class c *        |                 |
            // likelihood for this site-------------|                 |
            //                                                        |
            //Substitution function for site i and rate class c-------|
          }          
        }
      }
    }
    
    //Now we just have to copy the substitutions into the result vector:
    for(unsigned int i = 0; i < nbSites; i++)
    {
      (*substitutions)(l, i) = substitutionsForCurrentNode[(* rootPatternLinks)[i]] / Lr[(* rootPatternLinks)[i]];
    }
  }
  if(verbose)
  {
    if(ApplicationTools::message) *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
  return substitutions;
}

/**************************************************************************************************/

ProbabilisticSubstitutionMapping * SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(
  const DRTreeLikelihood & drtl,
  SubstitutionCount & substitutionCount,
  bool verbose) throw (Exception)
{
  //Preamble:
  if(!drtl.isInitialized()) throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsNoAveraging(). Likelihood object is not initialized.");
                                   
  //A few variables we'll need:
  const TreeTemplate<Node> *    tree = dynamic_cast<const TreeTemplate<Node> *>(drtl.getTree());
  const SiteContainer *    sequences = drtl.getData();
  const DiscreteDistribution * rDist = drtl.getRateDistribution();
    
  unsigned int nbSites         = sequences->getNumberOfSites();
  unsigned int nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  unsigned int nbStates        = sequences->getAlphabet()->getSize();
  unsigned int nbClasses       = rDist->getNumberOfCategories();
  vector<const Node *> nodes   = tree->getNodes();
  const vector<unsigned int> * rootPatternLinks
                               = & drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  unsigned int nbNodes = nodes.size();
  
  // We create a new ProbabilisticSubstitutionMapping object:
  ProbabilisticSubstitutionMapping * substitutions = new ProbabilisticSubstitutionMapping(*tree, nbSites);
                                   
  Vdouble rcRates = rDist->getCategories();

  // Compute the number of substitutions for each class and each branch in the tree:
  if(verbose) ApplicationTools::displayTask("Compute joint node-pairs likelihood", true);
  
  for(unsigned int l = 0; l < nbNodes; l++)
  {
    // For each node,
    const Node * currentNode = nodes[l];

    const Node * father = currentNode->getFather();

    double d = currentNode->getDistanceToFather();
    
    substitutionCount.setSubstitutionModel(drtl.getSubstitutionModelForNode(currentNode->getId()));
    
    if(verbose) ApplicationTools::displayGauge(l, nbNodes-1, '>');
    Vdouble substitutionsForCurrentNode(nbDistinctSites, 0);

    // Compute all nxy first:
    VVVdouble nxy = VVVdouble(nbClasses);
    for(unsigned int c = 0; c < nbClasses; c++)
    {
      VVdouble * nxy_c = & nxy[c];
      double rc = rcRates[c];
      Matrix<double> * nijt = substitutionCount.getAllNumbersOfSubstitutions(d * rc);
      nxy_c -> resize(nbStates);
      for(unsigned int x = 0; x < nbStates; x++)
      {
        Vdouble * nxy_c_x = & (* nxy_c)[x];
        nxy_c_x->resize(nbStates);
        for(unsigned int y = 0; y < nbStates; y++)
        {
          (* nxy_c_x)[y] = (* nijt)(x, y);
        }
      }
      delete nijt;
    }
    
    map<int, VVVdouble> likelihoodsFather = drtl.getLikelihoodData()->getLikelihoodArrays(father->getId());

    // Now we've got to compute likelihoods in a smart manner... ;)

    VVVdouble likelihoodsFatherConstantPart(nbDistinctSites);
    for(unsigned int i = 0; i < nbDistinctSites; i++)
    {
      VVdouble * likelihoodsFatherConstantPart_i = & likelihoodsFatherConstantPart[i];
      likelihoodsFatherConstantPart_i->resize(nbClasses);
      for(unsigned int c = 0; c < nbClasses; c++)
      {
        Vdouble * likelihoodsFatherConstantPart_i_c = & (* likelihoodsFatherConstantPart_i)[c];
        likelihoodsFatherConstantPart_i_c->resize(nbStates);
        double rc = rDist->getProbability(c);
        for(unsigned int s = 0; s < nbStates; s++)
        {
          //(* likelihoodsFatherConstantPart_i_c)[s] = rc * model->freq(s);
          //freq is already accounted in the array
          (* likelihoodsFatherConstantPart_i_c)[s] = rc;
        }
      }
    }
    
    // First, what will remain constant:
    unsigned int nbSons = father->getNumberOfSons();
    for(unsigned int n = 0; n < nbSons; n++)
    {
      const Node * currentSon = father->getSon(n);
      
      if(currentSon->getId() != currentNode->getId())
      {
        VVVdouble pxy = drtl.getTransitionProbabilitiesForNode(currentSon);
        VVVdouble * likelihoodsFather_son = & likelihoodsFather[currentSon->getId()];
        for(unsigned int i = 0; i < nbDistinctSites; i++)
        {
          VVdouble * likelihoodsFather_son_i = & (* likelihoodsFather_son)[i];
          VVdouble * likelihoodsFatherConstantPart_i = & likelihoodsFatherConstantPart[i];
          for(unsigned int c = 0; c < nbClasses; c++)
          {
            Vdouble * likelihoodsFather_son_i_c = & (* likelihoodsFather_son_i)[c];
            Vdouble * likelihoodsFatherConstantPart_i_c = & (* likelihoodsFatherConstantPart_i)[c];
            VVdouble * pxy_c = & pxy[c]; 
            for(unsigned int x = 0; x < nbStates; x++)
            {
              Vdouble * pxy_c_x = & (* pxy_c)[x];
              double likelihood = 0.;
              for(unsigned int y = 0; y < nbStates; y++)
              {
                likelihood += (* pxy_c_x)[y] * (* likelihoodsFather_son_i_c)[y];
              }
              (* likelihoodsFatherConstantPart_i_c)[x] *= likelihood;
            }
          }
        }      
      }
    }
    if(father -> hasFather())
    {
      const Node * currentSon = father->getFather();
      VVVdouble pxy = drtl.getTransitionProbabilitiesForNode(father);
      VVVdouble * likelihoodsFather_son = & likelihoodsFather[currentSon->getId()];
      for(unsigned int i = 0; i < nbDistinctSites; i++)
      {
        VVdouble * likelihoodsFather_son_i = & (* likelihoodsFather_son)[i];
        VVdouble * likelihoodsFatherConstantPart_i = & likelihoodsFatherConstantPart[i];
        for(unsigned int c = 0; c < nbClasses; c++)
        {
          Vdouble * likelihoodsFather_son_i_c = & (* likelihoodsFather_son_i)[c];
          Vdouble * likelihoodsFatherConstantPart_i_c = & (* likelihoodsFatherConstantPart_i)[c];
          VVdouble * pxy_c = & pxy[c]; 
          for(unsigned int x = 0; x < nbStates; x++)
          {
            Vdouble * pxy_c_x = & (* pxy_c)[x];
            double likelihood = 0.;
            for(unsigned int y = 0; y < nbStates; y++)
            {
              likelihood += (* pxy_c_x)[y] * (* likelihoodsFather_son_i_c)[y];
            }
            (* likelihoodsFatherConstantPart_i_c)[x] *= likelihood;
          }
        }
      }      
    }

    // Then, we deal with the node of interest.
    // We first average uppon 'y' to save computations, and then uppon 'x'.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)
    VVVdouble pxy = drtl.getTransitionProbabilitiesForNode(currentNode);
    VVVdouble * likelihoodsFather_node = & likelihoodsFather[currentNode->getId()];
    for(unsigned int i = 0; i < nbDistinctSites; i++)
    {
      VVdouble * likelihoodsFather_node_i = & (* likelihoodsFather_node)[i];
      VVdouble * likelihoodsFatherConstantPart_i = & likelihoodsFatherConstantPart[i];
      RowMatrix<double> pairProbabilities(nbStates, nbStates);
      MatrixTools::fill(pairProbabilities, 0.);
      RowMatrix<double> subsCounts(nbStates, nbStates);
      MatrixTools::fill(subsCounts, 0.);
      
      for(unsigned int c = 0; c < nbClasses; c++)
      {
        Vdouble * likelihoodsFather_node_i_c = & (* likelihoodsFather_node_i)[c];
        Vdouble * likelihoodsFatherConstantPart_i_c = & (* likelihoodsFatherConstantPart_i)[c];
        VVdouble * pxy_c = & pxy[c];
        VVdouble * nxy_c = & nxy[c];
        for(unsigned int x = 0; x < nbStates; x++)
        {
          double * likelihoodsFatherConstantPart_i_c_x = & (* likelihoodsFatherConstantPart_i_c)[x];
          Vdouble * pxy_c_x = & (* pxy_c)[x];
          Vdouble * nxy_c_x = & (* nxy_c)[x];
          for(unsigned int y = 0; y < nbStates; y++)
          {
            double likelihood_cxy = (* likelihoodsFatherConstantPart_i_c_x)
              * (* pxy_c_x)[y]
              * (* likelihoodsFather_node_i_c)[y];
            pairProbabilities(x, y) += likelihood_cxy; // Sum over all rate classes.
            subsCounts(x, y) += likelihood_cxy * (* nxy_c_x)[y];
            //cout << pairProbabilities(x, y) << "\t" << subsCounts(x, y) << endl;
          }
        }
      }
      // Now the vector computation:
      // Here we do not average over all possible pair of ancestral states,
      // We only consider the one with max likelihood:
      vector<unsigned int> xy = MatrixTools::whichmax(pairProbabilities);
      substitutionsForCurrentNode[i] += subsCounts(xy[0], xy[1]) / pairProbabilities(xy[0], xy[1]);
    }
    
    //Now we just have to copy the substitutions into the result vector:
    for(unsigned int i = 0; i < nbSites; i++)
    {
      (*substitutions)(l, i) = substitutionsForCurrentNode[(* rootPatternLinks)[i]];
    }
  }
  if(verbose)
  {
    if(ApplicationTools::message) *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
  return substitutions;
}

/**************************************************************************************************/

ProbabilisticSubstitutionMapping * SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(
  const DRTreeLikelihood & drtl,
  SubstitutionCount & substitutionCount,
  bool verbose) throw (Exception)
{
  //Preamble:
  if(!drtl.isInitialized()) throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsNoAveragingMarginal(). Likelihood object is not initialized.");
                                   
  //A few variables we'll need:
  
  const TreeTemplate<Node> *    tree = dynamic_cast<const TreeTemplate<Node> *>(drtl.getTree());
  const SiteContainer *    sequences = drtl.getData();
  const DiscreteDistribution * rDist = drtl.getRateDistribution();
  const Alphabet *             alpha = sequences->getAlphabet();
    
  unsigned int nbSites         = sequences->getNumberOfSites();
  unsigned int nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  vector<const Node *> nodes   = tree->getNodes();
  const vector<unsigned int> * rootPatternLinks
                               = &drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  unsigned int nbNodes = nodes.size();
  
  // We create a new ProbabilisticSubstitutionMapping object:
  ProbabilisticSubstitutionMapping * substitutions = new ProbabilisticSubstitutionMapping(*tree, nbSites);
  
  // Compute the whole likelihood of the tree according to the specified model:
  
  Vdouble rcRates = rDist->getCategories();

  // Compute the number of substitutions for each class and each branch in the tree:
  if(verbose) ApplicationTools::displayTask("Compute marginal ancestral states");
  MarginalAncestralStateReconstruction masr(drtl);
  map<int, vector<int> > ancestors = masr.getAllAncestralStates();
  if(verbose) ApplicationTools::displayTaskDone();

  // Now we just have to compute the substitution vectors:
  if(verbose) ApplicationTools::displayTask("Compute substitution vectors", true);
  
  for(unsigned int l = 0; l < nbNodes; l++)
  {
    const Node * currentNode = nodes[l];
    
    const Node * father = currentNode->getFather();

    substitutionCount.setSubstitutionModel(drtl.getSubstitutionModelForNode(currentNode->getId()));
    
    double d = currentNode->getDistanceToFather();

    vector<int> nodeStates = ancestors[currentNode->getId()]; //These are not 'true' ancestors ;)
    vector<int> fatherStates = ancestors[father->getId()];
    
    //For each node,
    if(verbose) ApplicationTools::displayGauge(l, nbNodes-1, '>');
    Vdouble substitutionsForCurrentNode(nbDistinctSites, 0.);

    //compute all nxy first:
    Matrix<double> * nxy = substitutionCount.getAllNumbersOfSubstitutions(d);
    
    // Here, we have no likelihood computation to do!

    // Then, we deal with the node of interest.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)
    for(unsigned int i = 0; i < nbDistinctSites; i++)
    {
      int fatherState = fatherStates[i];
      int nodeState   = nodeStates[i];
      if(fatherState >= (int)(alpha->getSize()) || nodeState >= (int)(alpha->getSize()))
        substitutionsForCurrentNode[i] = 0; // To be conservative! Only in case there are generic characters.
      else
        substitutionsForCurrentNode[i] = (* nxy)((unsigned int)fatherStates[i],(unsigned int)nodeStates[i]);
    }
    
    //Now we just have to copy the substitutions into the result vector:
    for(unsigned int i = 0; i < nbSites; i++)
    {
      (*substitutions)(l,i) = substitutionsForCurrentNode[(* rootPatternLinks)[i]];
    }
  }
  if(verbose)
  {
    if(ApplicationTools::message) *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
  return substitutions;
}

/**************************************************************************************************/

ProbabilisticSubstitutionMapping * SubstitutionMappingTools::computeSubstitutionVectorsMarginal(
  const DRTreeLikelihood & drtl,
  SubstitutionCount & substitutionCount,
  bool verbose) throw (Exception)
{
  //Preamble:
  if(!drtl.isInitialized()) throw Exception("SubstitutionMappingTools::computeSubstitutionVectorsMarginal(). Likelihood object is not initialized.");
                                   
  //A few variables we'll need:
  
  const TreeTemplate<Node> *    tree = dynamic_cast<const TreeTemplate<Node> *>(drtl.getTree());
  const SiteContainer *    sequences = drtl.getData();
  const DiscreteDistribution * rDist = drtl.getRateDistribution();
    
  unsigned int nbSites         = sequences->getNumberOfSites();
  unsigned int nbDistinctSites = drtl.getLikelihoodData()->getNumberOfDistinctSites();
  unsigned int nbStates        = sequences->getAlphabet()->getSize();
  unsigned int nbClasses       = rDist->getNumberOfCategories();
  vector<const Node *> nodes   = tree->getNodes();
  const vector<unsigned int> * rootPatternLinks
                               = &drtl.getLikelihoodData()->getRootArrayPositions();
  nodes.pop_back(); // Remove root node.
  unsigned int nbNodes = nodes.size();
  
  // We create a new ProbabilisticSubstitutionMapping object:
  ProbabilisticSubstitutionMapping * substitutions = new ProbabilisticSubstitutionMapping(*tree, nbSites);
                                     
  // Compute the whole likelihood of the tree according to the specified model:
  
  Vdouble rcProbs = rDist->getProbabilities();
  Vdouble rcRates = rDist->getCategories();

  //II) Compute the number of substitutions for each class and each branch in the tree:
  if(verbose) ApplicationTools::displayTask("Compute marginal node-pairs likelihoods", true);
  
  for(unsigned int l = 0; l < nbNodes; l++)
  {
    const Node * currentNode = nodes[l];
    
    const Node * father = currentNode->getFather();

    substitutionCount.setSubstitutionModel(drtl.getSubstitutionModelForNode(currentNode->getId()));
    
    double d = currentNode->getDistanceToFather();
    
    //For each node,
    if(verbose) ApplicationTools::displayGauge(l, nbNodes-1, '>');
    Vdouble substitutionsForCurrentNode(nbDistinctSites, 0);

    //compute all nxy first:
    VVVdouble nxy = VVVdouble(nbClasses);
    for(unsigned int c = 0; c < nbClasses; c++)
    {
      VVdouble * nxy_c = & nxy[c];
      double rc = rcRates[c];
      Matrix<double> * nijt = substitutionCount.getAllNumbersOfSubstitutions(d * rc);
      nxy_c -> resize(nbStates);
      for(unsigned int x = 0; x < nbStates; x++)
      {
        Vdouble * nxy_c_x = & (* nxy_c)[x];
        nxy_c_x->resize(nbStates);
        for(unsigned int y = 0; y < nbStates; y++)
        {
          (* nxy_c_x)[y] = (* nijt)(x, y);
        }
      }
      delete nijt;
    }

    // Then, we deal with the node of interest.
    // ('y' is the state at 'node' and 'x' the state at 'father'.)
    VVVdouble probsNode   = DRTreeLikelihoodTools::getPosteriorProbabilitiesForEachStateForEachRate(drtl, currentNode->getId());
    VVVdouble probsFather = DRTreeLikelihoodTools::getPosteriorProbabilitiesForEachStateForEachRate(drtl, father->getId());
    for(unsigned int i = 0; i < nbDistinctSites; i++)
    {
      VVdouble * probsNode_i   = & probsNode[i];
      VVdouble * probsFather_i = & probsFather[i];
      for(unsigned int c = 0; c < nbClasses; c++)
      {
        Vdouble * probsNode_i_c   = & (* probsNode_i)[c];
        Vdouble * probsFather_i_c = & (* probsFather_i)[c];
        VVdouble * nxy_c = & nxy[c];
        for(unsigned int x = 0; x < nbStates; x++)
        {
          Vdouble * nxy_c_x = & (* nxy_c)[x];
          for(unsigned int y = 0; y < nbStates; y++)
          {
            double prob_cxy = (* probsFather_i_c)[x] * (* probsNode_i_c)[y];
            // Now the vector computation:
            substitutionsForCurrentNode[i] += prob_cxy * (* nxy_c_x)[y];
            //                                <------>   <------------>
            // Posterior probability              |                | 
            // for site i and rate class c *      |                |
            // likelihood for this site-----------|                |
            //                                                     |
            //Substitution function for site i and rate class c----|
          }
        }
      }
    }
    
    //Now we just have to copy the substitutions into the result vector:
    for(unsigned int i = 0; i < nbSites; i++)
    {
      (*substitutions)(l,i) = substitutionsForCurrentNode[(* rootPatternLinks)[i]];
    }
  }
  if(verbose)
  {
    if(ApplicationTools::message) *ApplicationTools::message << " ";
    ApplicationTools::displayTaskDone();
  }
  return substitutions;
}

/**************************************************************************************************/

void SubstitutionMappingTools::writeToStream(
  const ProbabilisticSubstitutionMapping & substitutions,
  const SiteContainer & sites,
  ostream & out)
  throw (IOException) 
{
  if(!out) throw IOException("SubstitutionMappingTools::writeToFile. Can't write to stream.");
  out << "Branches";
  out << "\tMean";
  for(unsigned int i = 0; i < substitutions.getNumberOfSites(); i++)
  {
    out << "\tSite" << sites.getSite(i) -> getPosition();
  }
  out << endl;
  
  for(unsigned int j = 0; j < substitutions.getNumberOfBranches(); j++)
  {
    out << substitutions.getNode(j)->getId() << "\t" << substitutions.getNode(j)->getDistanceToFather();
    for(unsigned int i = 0; i < substitutions.getNumberOfSites(); i++)
    {
      out << "\t" << substitutions(j, i);
    }
    out << endl;
  }
}

/**************************************************************************************************/

void SubstitutionMappingTools::readFromStream(istream & in, ProbabilisticSubstitutionMapping & substitutions)
  throw (IOException)
{
  try {
    DataTable * data = DataTable::read(in, "\t", true, -1);
    vector<string> ids = data->getColumn(0);
    data->deleteColumn(0);//Remove ids
    data->deleteColumn(0);//Remove means
    //Now parse the table:
    unsigned int nbSites = data->getNumberOfColumns();
    substitutions.setNumberOfSites(nbSites);
    unsigned int nbBranches = data->getNumberOfRows();
    for(unsigned int i = 0; i < nbBranches; i++)
    {
      int id = TextTools::toInt(ids[i]);
      unsigned int br = substitutions.getNodeIndex(id);
      for(unsigned int j = 0; j < nbSites; j++)
      {
        substitutions(br, j) = TextTools::toDouble((*data)(i, j));
      }
    }
    //Parse the header:
    for(unsigned int i = 0; i < nbSites; i++)
    {
      string siteTxt = data->getColumnName(i);
      int site = 0;
      if(siteTxt.substr(0,4) == "Site") site = TextTools::to<int>(siteTxt.substr(4));
      else site = TextTools::to<int>(siteTxt);
      substitutions.setSitePosition(i, site);
    }
    
    delete data;
  }
  catch(Exception & e)
  {
    throw IOException(string("Bad input file. ") + e.what());
  }
}

/**************************************************************************************************/

