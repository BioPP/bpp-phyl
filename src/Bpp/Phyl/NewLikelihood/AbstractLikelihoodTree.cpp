//
// File: AbstractLikelihoodTree.cpp
// Created by: Laurent Guéguen
// Created on: mardi 23 juin 2015, à 09h 43
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

#include "AbstractLikelihoodTree.h"

#include "SubstitutionProcess.h"

using namespace bpp;
using namespace std;

AbstractLikelihoodTree::AbstractLikelihoodTree(const SubstitutionProcess& process) :
  rootPatternLinks_(), rootWeights_(), alphabet_(0),
  shrunkData_(0), nbSites_(0), nbStates_(0),
  nbClasses_(process.getNumberOfClasses()),
  vProbClass_(process.getClassProbabilities()),
  nbDistinctSites_(0)
{
}

AbstractLikelihoodTree::AbstractLikelihoodTree(const AbstractLikelihoodTree& atd) :
  rootPatternLinks_(atd.rootPatternLinks_),
  rootWeights_(atd.rootWeights_),
  alphabet_(atd.alphabet_),
  shrunkData_(0),
  nbSites_(atd.nbSites_), nbStates_(atd.nbStates_),
  nbClasses_(atd.nbClasses_),
  vProbClass_(atd.vProbClass_),
  nbDistinctSites_(atd.nbDistinctSites_)
{
  if (atd.shrunkData_.get())
    shrunkData_.reset(dynamic_cast<SiteContainer*>(atd.shrunkData_->clone()));
}

AbstractLikelihoodTree& AbstractLikelihoodTree::operator=(const AbstractLikelihoodTree& atd)
{
  rootPatternLinks_ = atd.rootPatternLinks_;
  rootWeights_      = atd.rootWeights_;
  alphabet_         = atd.alphabet_;
  nbSites_          = atd.nbSites_;
  nbStates_         = atd.nbStates_;
  nbClasses_        = atd.nbClasses_;
  vProbClass_       = atd.vProbClass_;
  
  nbDistinctSites_  = atd.nbDistinctSites_;
  if (atd.shrunkData_.get())
    shrunkData_.reset(dynamic_cast<SiteContainer*>(atd.shrunkData_->clone()));
  else
    shrunkData_.reset();

  return *this;
}

AbstractLikelihoodTree::~AbstractLikelihoodTree()
{
}


/******************************************************************************/

VVVdouble AbstractLikelihoodTree::getPosteriorProbabilitiesForEachStateForEachClass(int nodeId)
{
  VVVdouble vRes(nbClasses_);

  for (size_t i=0; i<nbClasses_; i++)
  {
    getNodeData(nodeId, i).getPosteriorProbabilitiesForEachState(vRes[i]);
    vRes[i]*=vProbClass_[i];
  }
  
  return vRes;
}

/******************************************************************************/

Vdouble AbstractLikelihoodTree::getPosteriorStateFrequencies(int nodeId)
{
  VVVdouble probs = getPosteriorProbabilitiesForEachStateForEachClass(nodeId);
  Vdouble freqs(getNumberOfStates());
  double sumw = 0, w;
  
  for (size_t j = 0; j < getNumberOfClasses(); j++)
  {
    for (size_t i = 0; i < getNumberOfDistinctSites(); i++)
    {
      w = getWeight(i);
      sumw += w;
      for (size_t k = 0; k < getNumberOfStates(); k++)
      {
        freqs[k] += probs[j][i][k] * w;
      }
    }
  }

  for (size_t k = 0; k < getNumberOfStates(); k++)
  {
    freqs[k] /= sumw;
  }
  return freqs;  
}



