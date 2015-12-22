//
// File: HmmPhyloEmissionProbabilities.cpp
// Created by: Laurent Guéguen
// Created on: mardi 24 septembre 2013, à 10h 09
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

#include "HmmPhyloEmissionProbabilities.h"

using namespace bpp;
using namespace std;

HmmPhyloEmissionProbabilities::HmmPhyloEmissionProbabilities(const HmmPhyloAlphabet* alphabet) :
  AbstractParametrizable(""),
  phylAlph_(alphabet),
  emProb_(alphabet->getNumberOfSites()),
  dEmProb_(),
  d2EmProb_(),
  upToDate_(false),
  nbSites_(alphabet->getNumberOfSites())
{
  for (size_t i=0;i<emProb_.size();i++)
    emProb_[i].resize(getNumberOfStates());
}


void HmmPhyloEmissionProbabilities::setHmmStateAlphabet(const HmmStateAlphabet* stateAlphabet) throw (HmmUnvalidAlphabetException)
{
  if (stateAlphabet==NULL)
    throw HmmUnvalidAlphabetException("Null alphabet in HmmPhyloEmissionProbabilities::setHmmStateAlphabet");
  if (dynamic_cast<const HmmPhyloAlphabet*>(stateAlphabet)==NULL)
    throw HmmUnvalidAlphabetException("Non PhyloLikelihood alphabet in HmmPhyloEmissionProbabilities::setHmmStateAlphabet");
  
  phylAlph_=dynamic_cast<const HmmPhyloAlphabet*>(stateAlphabet);
  nbSites_=phylAlph_->getNumberOfSites();
  emProb_.resize(nbSites_);
  for (size_t i=0;i<emProb_.size();i++)
    emProb_[i].resize(getNumberOfStates());
  upToDate_=false;
}


void HmmPhyloEmissionProbabilities::computeEmissionProbabilities_() const
{
  phylAlph_->updateLikelihood();
  phylAlph_->computeLikelihood();
  
  for (size_t i=0;i<nbSites_;i++)
    for (size_t j=0;j<getNumberOfStates();j++)
      emProb_[i][j]= phylAlph_->getPhyloLikelihood(j).getLikelihoodForASite(i);

  upToDate_=true;
}

void HmmPhyloEmissionProbabilities::computeDEmissionProbabilities(std::string& variable) const
{
  phylAlph_->computeDLogLikelihood(variable);
  
  if (dEmProb_.size()!=nbSites_)
  {
    dEmProb_.resize(nbSites_);
    for (size_t i=0; i<nbSites_;i++)
      dEmProb_[i].resize(getNumberOfStates());
  }

  for (size_t i=0;i<nbSites_;i++)
    for (size_t j=0;j<getNumberOfStates();j++)
    {
      const AlignedPhyloLikelihood& apl=phylAlph_->getPhyloLikelihood(j);
      
      dEmProb_[i][j]= apl.getDLogLikelihoodForASite(i) * apl.getLikelihoodForASite(i);
    }
}
  
void HmmPhyloEmissionProbabilities::computeD2EmissionProbabilities(std::string& variable) const
{
  phylAlph_->computeD2LogLikelihood(variable);
  
  if (d2EmProb_.size()!=nbSites_)
  {
    d2EmProb_.resize(nbSites_);
    for (size_t i=0; i<nbSites_;i++)
      d2EmProb_[i].resize(getNumberOfStates());
  }

  for (size_t i=0;i<nbSites_;i++)
    for (size_t j=0;j<getNumberOfStates();j++)
    {
      const AlignedPhyloLikelihood& apl=phylAlph_->getPhyloLikelihood(j);
      double x= apl.getDLogLikelihoodForASite(i);
      
      d2EmProb_[i][j]= (apl.getD2LogLikelihoodForASite(i) + x*x) * apl.getLikelihoodForASite(i);
    }
}


  

