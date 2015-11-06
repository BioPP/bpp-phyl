//
// File: MultiProcessSequenceEvolution.cpp
// Created by: Laurent Guéguen
// Created on: mardi 28 avril 2015, à 15h 18
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

#include "MultiProcessSequenceEvolution.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

MultiProcessSequenceEvolution::MultiProcessSequenceEvolution(
  SubstitutionProcessCollection* processColl,
  vector<size_t> nProc,
  const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  processColl_(processColl),
  nProc_(nProc)
{
  // initialize parameters:

  for (size_t i=0; i<nProc_.size(); i++)
    includeParameters_(processColl_->getSubstitutionProcessParameters(nProc_[i],true));

//  setNamespace(prefix);
}

/******************************************************************************/

void MultiProcessSequenceEvolution::fireParameterChanged(const ParameterList& parameters)
{
  AbstractParameterAliasable::fireParameterChanged(parameters);
  processColl_->matchParametersValues(parameters);
}

/******************************************************************************/

ParameterList MultiProcessSequenceEvolution::getSubstitutionProcessParameters(bool independent) const
{
  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcessParameters(nProc_[i], independent));

  return pl;
}

ParameterList MultiProcessSequenceEvolution::getSubstitutionModelParameters(bool independent) const
{ 
  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcess(nProc_[i]).getSubstitutionModelParameters(independent));

  return pl;
}

ParameterList MultiProcessSequenceEvolution::getRateDistributionParameters(bool independent) const
{
  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcess(nProc_[i]).getRateDistributionParameters(independent));

  return pl;
}

ParameterList MultiProcessSequenceEvolution::getRootFrequenciesParameters(bool independent) const
{
  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcess(nProc_[i]).getRootFrequenciesParameters(independent));

  return pl;
}

ParameterList MultiProcessSequenceEvolution::getBranchLengthParameters(bool independent) const
{ 
  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcess(nProc_[i]).getBranchLengthParameters(independent));

  return pl;
}

ParameterList MultiProcessSequenceEvolution::getNonDerivableParameters() const
{
  // patch, to be fixed properly later
  return getIndependentParameters();

  ParameterList pl;
  
  for (size_t i=0; i<nProc_.size(); i++)
    pl.includeParameters(processColl_->getSubstitutionProcess(nProc_[i]).getNonDerivableParameters());

  return pl;
}

/******************************************************************************/


void MultiProcessSequenceEvolution::setParameters(const ParameterList& parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  setParametersValues(parameters);
}


/*************************************************************************/

bool MultiProcessSequenceEvolution::isCompatibleWith(const SiteContainer& data) const
{
  for (size_t i=0; i<nProc_.size(); i++)
    if ( !processColl_->getSubstitutionProcess(nProc_[i]).isCompatibleWith(data))
      return false;
  
  return true;
}

/*****************************************************************************/

bool MultiProcessSequenceEvolution::hasDerivableParameter(const std::string& name) const
{
  for (size_t i=0; i<nProc_.size(); i++)
    if (processColl_->getSubstitutionProcess(nProc_[i]).hasDerivableParameter(name))
      return true;
  
  return false;
}

