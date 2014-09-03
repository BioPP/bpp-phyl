// File: DoubleRecursiveTreeLikelihoodCalculation.cpp
// Created by: Julien Dutheil
// Created on: Fri Apr 26 20:07 2013
// From file: DRNonHomogeneousTreeLikelihood.cpp
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

#include "DoubleRecursiveTreeLikelihoodCalculation.h"
#include "../PatternTools.h"
#include "ComputingNode.h"

using namespace std;
using namespace bpp;
using namespace newlik;

/******************************************************************************/

DoubleRecursiveTreeLikelihoodCalculation::DoubleRecursiveTreeLikelihoodCalculation(
  SubstitutionProcess* process,
  bool verbose)
throw (Exception) :
  AbstractTreeLikelihoodCalculation(process, verbose),
  likelihoodData_(0),
  root1_(-1),
  root2_(-1),
  nullDLikelihood_(true),
  nullD2Likelihood_(true),
  compNId_(-1)
{
  init_();
}

/******************************************************************************/

DoubleRecursiveTreeLikelihoodCalculation::DoubleRecursiveTreeLikelihoodCalculation(
  const SiteContainer& data,
  SubstitutionProcess* process,
  bool verbose)
throw (Exception) :
  AbstractTreeLikelihoodCalculation(process, verbose),
  likelihoodData_(0),
  root1_(-1),
  root2_(-1),
  nullDLikelihood_(true),
  nullD2Likelihood_(true),
  compNId_(-1)
{
  init_();
  setData(data);
}

/******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::init_() throw (Exception)
{
  likelihoodData_.reset(new DoubleRecursiveTreeLikelihoodData(
                          process_->getNumberOfClasses()));
}

/******************************************************************************/

DoubleRecursiveTreeLikelihoodCalculation::DoubleRecursiveTreeLikelihoodCalculation(const DoubleRecursiveTreeLikelihoodCalculation& lik) :
  AbstractTreeLikelihoodCalculation(lik),
  likelihoodData_(0),
  root1_(lik.root1_),
  root2_(lik.root2_),
  nullDLikelihood_(lik.nullDLikelihood_),
  nullD2Likelihood_(lik.nullD2Likelihood_),
  compNId_(lik.compNId_)
{
  likelihoodData_.reset(lik.likelihoodData_->clone());
}

/******************************************************************************/

DoubleRecursiveTreeLikelihoodCalculation& DoubleRecursiveTreeLikelihoodCalculation::operator=(const DoubleRecursiveTreeLikelihoodCalculation& lik)
{
  AbstractTreeLikelihoodCalculation::operator=(lik);
  likelihoodData_.reset(lik.likelihoodData_->clone());
  root1_=lik.root1_;
  root2_=lik.root2_;
  
  nullDLikelihood_=lik.nullDLikelihood_;
  nullD2Likelihood_=lik.nullD2Likelihood_;
  compNId_=lik.compNId_;
  
  return *this;
}

/******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::setData(const SiteContainer& sites) throw (Exception)
{
  //try {
    const TreeTemplate<Node>& tt = dynamic_cast<const TreeTemplate<Node>&>(process_->getTree());
    data_.reset(PatternTools::getSequenceSubset(sites, *tt.getRootNode()));
  // } catch (exception& e) {
  //   throw Exception("DEBUG. DoubleRecursiveTreeLikelihoodCalculation::setData. The SubstitutionProcess does not use a TreeTemplate object.");
  // }
  if (verbose_)
    ApplicationTools::displayTask("Initializing data structure");
  likelihoodData_->initLikelihoods(*data_, *process_); //We assume here that all models have the same number of states, and that they have the same 'init' method,
                                                       //Which is a reasonable assumption as long as they share the same alphabet.
  if (verbose_)
    ApplicationTools::displayTaskDone();

  nbSites_         = likelihoodData_->getNumberOfSites();
  nbDistinctSites_ = likelihoodData_->getNumberOfDistinctSites();
  nbStates_        = likelihoodData_->getNumberOfStates();

  if (verbose_)
    ApplicationTools::displayResult("Number of distinct sites", TextTools::toString(nbDistinctSites_));

  initialized_ = true;
  compNId_ = -1;
  
  //Recompute likelihood:
  computeTreeLikelihood();
}

/******************************************************************************/

double DoubleRecursiveTreeLikelihoodCalculation::getLikelihoodForASite(size_t site) const
{
  double l = 0;
  VVVdouble* lla = &likelihoodData_->getRootLikelihoodArray();
  size_t posR=likelihoodData_->getRootArrayPosition(site);

  
  for (size_t c = 0; c < nbClasses_; ++c)
  {
    for (size_t j = 0; j < nbStates_; ++j)
    {
      l += (*lla)[c][posR][j] * process_->getProbabilityForModel(c) * process_->getRootFrequencies()[j];
    }
  }
  
  if (l < 0) l = 0; //May happen because of numerical errors.
  return l;
}

/******************************************************************************/

double DoubleRecursiveTreeLikelihoodCalculation::getLikelihoodForASiteForAClass(size_t site, size_t classIndex) const
{
  size_t posR=likelihoodData_->getRootArrayPosition(site);
  Vdouble* la = &likelihoodData_->getRootLikelihoodArray()[classIndex][posR];

  double l = 0;
  for (size_t i = 0; i < nbStates_; ++i)
  {
    l += (*la)[i] * process_->getRootFrequencies()[i];
  }
  return l;
}

/******************************************************************************/

double DoubleRecursiveTreeLikelihoodCalculation::getLikelihoodForASiteForAState(size_t site, int state) const
{
  double l = 0;
  VVVdouble* lla = &likelihoodData_->getRootLikelihoodArray();
  size_t posR=likelihoodData_->getRootArrayPosition(site);

  for (size_t c = 0; c < nbClasses_; ++c)
  {
    l += (*lla)[c][posR][state] * process_->getProbabilityForModel(c);
  }
  return l;
}

/******************************************************************************/

double DoubleRecursiveTreeLikelihoodCalculation::getLikelihoodForASiteForAClassForAState(size_t site, size_t classIndex, int state) const
{
  return likelihoodData_->getRootLikelihoodArray()[classIndex][likelihoodData_->getRootArrayPosition(site)][state];
}



/******************************************************************************
 *                           Likelihood computation                           *
 ******************************************************************************/

/******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::computeTreeLikelihood()
{
  computeSubtreeLikelihoodPostfix_(process_->getTree().getRootNode());
  computeSubtreeLikelihoodPrefix_(process_->getTree().getRootNode());
  computeRootLikelihood();
}

/******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::computeSubtreeLikelihoodPostfix_(const Node* node)
{
  if (node->getNumberOfSons() == 0)
    return;

  // Set all likelihood arrays to 1 for a start:
  likelihoodData_->getNodeData(node->getId()).resetLikelihoodArrays();

  
  map<int, VVVdouble>* likelihoods_node = &likelihoodData_->getLikelihoodArrays(node->getId());
  
  size_t nbNodes = node->getNumberOfSons();

  for (size_t l = 0; l < nbNodes; l++)
  {
    // For each son node...
    
    const Node* son = node->getSon(l);
    VVVdouble* likelihoods_node_son = &(*likelihoods_node)[son->getId()];

    if (son->isLeaf())
    {
      VVdouble* likelihoods_leaf = &likelihoodData_->getLeafLikelihoods(son->getId());

      // For each rate classe,
      for (size_t c = 0; c < nbClasses_; c++)
      {
        VVdouble* likelihoods_node_son_c = &(*likelihoods_node_son)[c];
        
        for (size_t i = 0; i < nbDistinctSites_; i++)
        {
          // For each site in the sequence,
          Vdouble* likelihoods_leaf_i = &(*likelihoods_leaf)[i];
          Vdouble* likelihoods_node_son_c_i = &(*likelihoods_node_son_c)[i];
          for (size_t x = 0; x < nbStates_; x++)
          {
            // For each initial state,
            (*likelihoods_node_son_c_i)[x] = (*likelihoods_leaf_i)[x];
          }
        }
      }
    }
    else
    {
      computeSubtreeLikelihoodPostfix_(son); // Recursive method:

      size_t nbSons = son->getNumberOfSons();

      map<int, VVVdouble>* likelihoods_son = &likelihoodData_->getLikelihoodArrays(son->getId());

      vector<const VVVdouble*> iLik(nbSons);
      
      for (size_t n = 0; n < nbSons; n++)
      {
        const Node* sonSon = son->getSon(n);
        iLik[n] = &(*likelihoods_son)[sonSon->getId()];
      }
        
      process_->multiplyPartialLikelihoods(likelihoods_node_son, iLik, son->getId(), ComputingNode::D0);
    }
  }
}


/******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::computeSubtreeLikelihoodPrefix_(const Node* node)
{
  if (!node->hasFather())
  {
    // 'node' is the root of the tree.
    // Just call the method on each son node:
    size_t nbSons = node->getNumberOfSons();
    for (size_t n = 0; n < nbSons; n++)
      computeSubtreeLikelihoodPrefix_(node->getSon(n));

    return;
  }
  else
  {
    const Node* father = node->getFather();
    map<int, VVVdouble>* likelihoods_node = &likelihoodData_->getLikelihoodArrays(node->getId());
    map<int, VVVdouble>* likelihoods_father = &likelihoodData_->getLikelihoodArrays(father->getId());
    VVVdouble* likelihoods_node_father = &(*likelihoods_node)[father->getId()];
    if (node->isLeaf())
    {
      TreeLikelihoodData::resetLikelihoodArray(likelihoodData_->getLikelihoodArrays(node->getId())[father->getId()]);
    }

    if (father->isLeaf())
    {
      // If the tree is rooted by a leaf
      VVdouble* likelihoods_leaf = &likelihoodData_->getLeafLikelihoods(father->getId());
      for (size_t c = 0; c < nbClasses_; c++)
      {
        // For each rate classe,
        VVdouble* likelihoods_node_father_c = &(*likelihoods_node_father)[c];

        for (size_t i = 0; i < nbDistinctSites_; i++)
        {
          // For each site in the sequence,
          Vdouble* likelihoods_leaf_i = &(*likelihoods_leaf)[i];
          Vdouble* likelihoods_node_father_c_i = &(*likelihoods_node_father_c)[i];
          for (size_t x = 0; x < nbStates_; x++)
          {
            // For each initial state,
            (*likelihoods_node_father_c_i)[x] = (*likelihoods_leaf_i)[x];
          }
        }
      }
    }
    else
    {
      vector<const Node*> nodes;
      // Add brothers:
      size_t nbFatherSons = father->getNumberOfSons();
      vector<const VVVdouble*> iLik(nbFatherSons);

      for (size_t n = 0; n < nbFatherSons; n++)
      {
        const Node* son = father->getSon(n);
        if (son->getId() != node->getId())
          iLik[n] = &(*likelihoods_father)[son->getId()];
        else
          iLik[n] = 0;
      }

      process_->multiplyPartialLikelihoods(likelihoods_node_father, iLik, father->getId(), ComputingNode::D0);

      if (father->hasFather())
      {
        const Node* fatherFather = father->getFather();
        process_->multiplyPartialLikelihoods(likelihoods_node_father, &(*likelihoods_father)[fatherFather->getId()], father->getId(), ComputingNode::D0);
      }
    }

    if (!father->hasFather())
    {
      const vector<double>& rootFreqs=process_->getRootFrequencies();
      
      // We have to account for the root frequencies:
      for (size_t c = 0; c < nbClasses_; c++)
      { 
        VVdouble* likelihoods_node_father_c = &(*likelihoods_node_father)[c];
        for (size_t i = 0; i < nbDistinctSites_; i++)
        {
          Vdouble* likelihoods_node_father_c_i = &(*likelihoods_node_father_c)[i];
          for (size_t x = 0; x < nbStates_; x++)
          {
            (*likelihoods_node_father_c_i)[x] *= rootFreqs[x];
          }
        }
      }
    }

    // Call the method on each son node:
    size_t nbNodeSons = node->getNumberOfSons();
    for (size_t i = 0; i < nbNodeSons; i++)
    {
      computeSubtreeLikelihoodPrefix_(node->getSon(i)); // Recursive method.
    }
  }
}

/******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::computeRootLikelihood()
{
  const Node* root = process_->getTree().getRootNode();
  VVVdouble* rootLikelihoods = &likelihoodData_->getRootLikelihoodArray();
  // Set all likelihoods to 1 for a start:
  if (root->isLeaf())
  {
    VVdouble* leavesLikelihoods_root = &likelihoodData_->getLeafLikelihoods(root->getId());
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* rootLikelihoods_c = &(*rootLikelihoods)[c];
      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        Vdouble* leavesLikelihoods_root_i = &(*leavesLikelihoods_root)[i];
        Vdouble* rootLikelihoods_c_i = &(*rootLikelihoods_c)[i];
        for (size_t x = 0; x < nbStates_; x++)
        {
          (*rootLikelihoods_c_i)[x] = (*leavesLikelihoods_root_i)[x];
        }
      }
    }
  }
  else
  {
    TreeLikelihoodData::resetLikelihoodArray(*rootLikelihoods);
  }

  map<int, VVVdouble>* likelihoods_root = &likelihoodData_->getLikelihoodArrays(root->getId());
  size_t nbNodes = root->getNumberOfSons();
  vector<const VVVdouble*> iLik(nbNodes);
  for (size_t n = 0; n < nbNodes; n++)
  {
    const Node* son = root->getSon(n);
    iLik[n] = &(*likelihoods_root)[son->getId()];
  }

  process_->multiplyPartialLikelihoods(rootLikelihoods, iLik, root->getId(), ComputingNode::D0);

  VVdouble* rootLikelihoodsS  = &likelihoodData_->getRootStateLikelihoodArray();
  Vdouble* rootLikelihoodsSC = &likelihoodData_->getRootStateClassLikelihoodArray();

  const vector<double>& rootFreqs=process_->getRootFrequencies();

  for (size_t c = 0; c < nbClasses_; c++)
  {
    // For each class
    VVdouble* rootLikelihoods_c = &(*rootLikelihoods)[c];
    Vdouble* rootLikelihoodsS_c = &(*rootLikelihoodsS)[c];
    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      // For each site in the sequence,
      Vdouble* rootLikelihoods_c_i = &(*rootLikelihoods_c)[i];
      double* rootLikelihoodsS_c_i = &(*rootLikelihoodsS_c)[i];
      (*rootLikelihoodsSC)[i] = 0;
      (*rootLikelihoodsS_c_i) = 0;
      for (size_t x = 0; x < nbStates_; x++)
      {
        // For each initial state,
        (*rootLikelihoodsS_c_i) += rootFreqs[x] * (*rootLikelihoods_c_i)[x];
      }
      (*rootLikelihoodsSC)[i] += process_->getProbabilityForModel(c) * (*rootLikelihoodsS_c_i);
    }

    // Final checking (for numerical errors):
    for (size_t i = 0; i < nbDistinctSites_; i++)
      if ((*rootLikelihoodsSC)[i] < 0)
        (*rootLikelihoodsSC)[i] = 0.;
  }
}


/******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::computeLikelihoodAtNode_(const Node* node, VVVdouble& likelihoodArray, const Node* sonNode) const
{
  int nodeId = node->getId();
  likelihoodArray.resize(nbClasses_);
  map<int, VVVdouble>* likelihoods_node = &likelihoodData_->getLikelihoodArrays(nodeId);

  // Initialize likelihood array:
  if (node->isLeaf())
  {
    VVdouble* leavesLikelihoods_node = &likelihoodData_->getLeafLikelihoods(nodeId);
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* likelihoodArray_c = &likelihoodArray[c];
      likelihoodArray_c->resize(nbDistinctSites_);

      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        Vdouble* leavesLikelihoods_node_i = &(*leavesLikelihoods_node)[i];
        Vdouble* likelihoodArray_c_i = &(*likelihoodArray_c)[i];
        likelihoodArray_c_i->resize(nbStates_);
        for (size_t x = 0; x < nbStates_; x++)
        {
          (*likelihoodArray_c_i)[x] = (*leavesLikelihoods_node_i)[x];
        }
      }
    }
  }
  else
  {
    // Otherwise:
    // Set all likelihoods to 1 for a start:
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* likelihoodArray_c = &likelihoodArray[c];
      likelihoodArray_c->resize(nbDistinctSites_);
      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        Vdouble* likelihoodArray_c_i = &(*likelihoodArray_c)[i];
        likelihoodArray_c_i->resize(nbStates_);
        for (size_t x = 0; x < nbStates_; x++)
        {
          (*likelihoodArray_c_i)[x] = 1.;
        }
      }
    }
  }

  size_t nbNodes = node->getNumberOfSons();

  vector<const VVVdouble*> iLik(nbNodes);
  
  for (size_t n = 0; n < nbNodes; n++)
  {
    const Node* son = node->getSon(n);

    iLik[n] = (son != sonNode)?&(*likelihoods_node)[son->getId()]:0;
  }

  process_->multiplyPartialLikelihoods(&likelihoodArray, iLik, node->getId(), ComputingNode::D0);

  const vector<double>& rootFreqs=process_->getRootFrequencies();

  if (node->hasFather())
  {
    const Node* father = node->getFather();

    process_->multiplyPartialLikelihoods(&likelihoodArray, &(*likelihoods_node)[father->getId()], node->getId(), ComputingNode::D0);
  }
  else
  {
    // We have to account for the equilibrium frequencies:
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* likelihoodArray_c = &likelihoodArray[c];
      for (size_t i = 0; i < nbDistinctSites_; i++)
      {
        Vdouble* likelihoodArray_c_i = &(*likelihoodArray_c)[i];
        for (size_t x = 0; x < nbStates_; x++)
        {
          (*likelihoodArray_c_i)[x] *= rootFreqs[x];
        }
      }
    }
  }
}


/******************************************************************************
 *                           First Order Derivatives                          *
 ******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::computeTreeDLikelihoodAtNode(const Node* node)
{
  const Node* father = node->getFather();
  VVVdouble* likelihoods_father_node = &likelihoodData_->getLikelihoodArray(father->getId(), node->getId());
  
  VVdouble* dLikelihoods_node = &likelihoodData_->getDLikelihoodArray(node->getId());

  VVVdouble larray;
  computeLikelihoodAtNode_(father, larray, node);

  process_->multiplyPartialLikelihoods(&larray, likelihoods_father_node, node->getId(), ComputingNode::D1);

  for (size_t c = 0; c < nbClasses_; c++)
  {
    VVdouble* larray_c = &larray[c];
    Vdouble* dLikelihoods_node_c=&(*dLikelihoods_node)[c];

    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      double* dLikelihoods_node_c_i=&(*dLikelihoods_node_c)[i];
      Vdouble* larray_c_i = &(*larray_c)[i];
      (*dLikelihoods_node_c_i) = 0;
      for (size_t x = 0; x < nbStates_; x++)
        (*dLikelihoods_node_c_i) += (*larray_c_i)[x];
    }
  }
}

/******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::computeTreeDLogLikelihood(const std::string& variable)
{
  int brId;
  try {
    brId = atoi(variable.substr(5).c_str());
  }
  catch (std::exception const& e)
  {
    nullDLikelihood_=true;
    return;
  }

  compNId_=brId;

  nullDLikelihood_=false;
  const Node* branch = process_->getTree().getNode(brId);

  computeTreeDLikelihoodAtNode(branch);
}

double DoubleRecursiveTreeLikelihoodCalculation::getDLikelihoodForASite(size_t site) const
{
  if ((nullDLikelihood_) || (compNId_==-1))
    return 0;
  
  // Derivative of the sum is the sum of derivatives:
  double dl = 0;
  VVdouble* ldla = &likelihoodData_->getDLikelihoodArray(compNId_);
  size_t posR=likelihoodData_->getRootArrayPosition(site);

  for (size_t c = 0; c < nbClasses_; c++)
  {
    dl += (*ldla)[c][posR] * process_->getProbabilityForModel(c);
  }
  return  dl;
}

/******************************************************************************
 *                           Second Order Derivatives                         *
 ******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::computeTreeD2LikelihoodAtNode(const Node* node)
{
  const Node* father = node->getFather();
  VVVdouble* likelihoods_father_node = &likelihoodData_->getLikelihoodArray(father->getId(), node->getId());
  
  VVdouble* d2Likelihoods_node = &likelihoodData_->getD2LikelihoodArray(node->getId());

  VVVdouble larray;
  
  computeLikelihoodAtNode_(father, larray, node);

  process_->multiplyPartialLikelihoods(&larray, likelihoods_father_node, node->getId(), ComputingNode::D2);

  for (size_t c = 0; c < nbClasses_; c++)
  {
    VVdouble* larray_c = &larray[c];
    Vdouble* d2Likelihoods_node_c=&(*d2Likelihoods_node)[c];

    for (size_t i = 0; i < nbDistinctSites_; i++)
    {
      double* d2Likelihoods_node_c_i=&(*d2Likelihoods_node_c)[i];
      Vdouble* larray_c_i = &(*larray_c)[i];
      (*d2Likelihoods_node_c_i) = 0;
      for (size_t x = 0; x < nbStates_; x++)
        (*d2Likelihoods_node_c_i) += (*larray_c_i)[x];
    }
  }
}

/******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::computeTreeD2LogLikelihood(const std::string& variable)
{
  int brId;
  try {
    brId = atoi(variable.substr(5).c_str());
  }
  catch (std::exception const& e)
  {
    nullD2Likelihood_=true;
    return;
  }

  compNId_=brId;
  nullD2Likelihood_=false;
  const Node* branch = process_->getTree().getNode(brId);

  computeTreeD2LikelihoodAtNode(branch);
}

double DoubleRecursiveTreeLikelihoodCalculation::getD2LikelihoodForASite(size_t site) const
{
  if ((nullDLikelihood_) || (compNId_==-1))
    return 0;
  
  // Derivative of the sum is the sum of derivatives:
  double d2l = 0;
  VVdouble* ldla = &likelihoodData_->getD2LikelihoodArray(compNId_);
  size_t posR=likelihoodData_->getRootArrayPosition(site);

  for (size_t c = 0; c < nbClasses_; c++)
  {
    d2l += (*ldla)[c][posR] * process_->getProbabilityForModel(c);
  }
  return d2l;
}

/******************************************************************************/

void DoubleRecursiveTreeLikelihoodCalculation::displayLikelihood(const Node* node)
{
  cout << "Likelihoods at node " << node->getId() << ": " << endl;
  for (size_t n = 0; n < node->getNumberOfSons(); n++)
  {
    const Node* subNode = node->getSon(n);
    cout << "Array for sub-node " << subNode->getId() << endl;
    displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId(), subNode->getId()));
  }
  if (node->hasFather())
  {
    const Node* father = node->getFather();
    cout << "Array for father node " << father->getId() << endl;
    displayLikelihoodArray(likelihoodData_->getLikelihoodArray(node->getId(), father->getId()));
  }
  cout << "                                         ***" << endl;
}

/*******************************************************************************/

