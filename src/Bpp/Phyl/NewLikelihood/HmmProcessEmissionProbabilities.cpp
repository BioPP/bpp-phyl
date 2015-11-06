//
// File: HmmEmissionProbabilities.cpp
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

#include "HmmProcessEmissionProbabilities.h"

using namespace bpp;
using namespace std;

HmmProcessEmissionProbabilities::HmmProcessEmissionProbabilities(const HmmProcessAlphabet* alphabet,
                                                                 const MultiProcessSequencePhyloLikelihood* multiPL) :
  AbstractParametrizable(""),
  procAlph_(alphabet),
  multiPL_(multiPL),
  emProb_(multiPL->getNumberOfSites()),
  dEmProb_(),
  d2EmProb_(),
  upToDate_(false)
{
  for (size_t i=0;i<emProb_.size();i++)
    emProb_[i].resize(multiPL_->getNumberOfSubstitutionProcess());

  addParameters_(multiPL_->getSequenceEvolution().getSubstitutionProcessParameters(true));
}


void HmmProcessEmissionProbabilities::setHmmStateAlphabet(const HmmStateAlphabet* stateAlphabet) throw (HmmUnvalidAlphabetException)
{
  if (stateAlphabet==NULL)
    throw HmmUnvalidAlphabetException("Null alphabet in HmmProcessEmissionProbabilities::setHmmStateAlphabet");
  if (dynamic_cast<const HmmProcessAlphabet*>(stateAlphabet)==NULL)
    throw HmmUnvalidAlphabetException("Non Process alphabet in HmmProcessEmissionProbabilities::setHmmStateAlphabet");
  
  procAlph_=dynamic_cast<const HmmProcessAlphabet*>(stateAlphabet);
  emProb_.resize(multiPL_->getNumberOfSites());
  upToDate_=false;
}


void HmmProcessEmissionProbabilities::updateEmissionProbabilities_() const
{
  multiPL_->computeLikelihood();
  
  for (size_t i=0;i<emProb_.size();i++)
    for (size_t j=0;j<emProb_[j].size();j++)
      emProb_[i][j]= multiPL_->getLikelihoodForASiteForAProcess(i, j);
  
  upToDate_=true;
}

void HmmProcessEmissionProbabilities::computeDEmissionProbabilities(std::string& variable) const
{
  for (size_t i=0;i<multiPL_->getNumberOfSubstitutionProcess();i++)
    multiPL_->computeDLogLikelihoodForAProcess(variable,i);

  if (dEmProb_.size()!=multiPL_->getNumberOfSites())
  {
    dEmProb_.resize(multiPL_->getNumberOfSites());
    for (size_t i=0; i<multiPL_->getNumberOfSites();i++)
      dEmProb_[i].resize(multiPL_->getNumberOfSubstitutionProcess());
  }

  for (size_t i=0;i<dEmProb_.size();i++)
    for (size_t j=0;j<dEmProb_[j].size();j++)
      dEmProb_[i][j]= multiPL_->getDLogLikelihoodForASiteForAProcess(i, j) * multiPL_->getLikelihoodForASiteForAProcess(i, j);
}
  
void HmmProcessEmissionProbabilities::computeD2EmissionProbabilities(std::string& variable) const
{
  for (size_t i=0;i<multiPL_->getNumberOfSubstitutionProcess();i++)
    multiPL_->computeD2LogLikelihoodForAProcess(variable,i);

  if (d2EmProb_.size()!=multiPL_->getNumberOfSites())
  {
    d2EmProb_.resize(multiPL_->getNumberOfSites());
    for (size_t i=0; i<multiPL_->getNumberOfSites();i++)
      d2EmProb_[i].resize(multiPL_->getNumberOfSubstitutionProcess());
  }

  for (size_t i=0;i<d2EmProb_.size();i++)
    for (size_t j=0;j<d2EmProb_[j].size();j++){
      double x= multiPL_->getDLogLikelihoodForASiteForAProcess(i, j);
      d2EmProb_[i][j]= (multiPL_->getD2LogLikelihoodForASiteForAProcess(i, j) + x*x)*multiPL_->getLikelihoodForASiteForAProcess(i, j);
    }
}

  

