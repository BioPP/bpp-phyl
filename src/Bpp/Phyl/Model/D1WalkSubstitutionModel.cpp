//
// File: D1WalkSubstitutionModel.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: mercredi 13 juillet 2016, ÃÂ  08h 55
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/


#include "D1WalkSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

#include <cmath>
#include <map>

using namespace std;

/******************************************************************************/

D1WalkSubstitutionModel::D1WalkSubstitutionModel(std::shared_ptr<const IntegerAlphabet> alpha, unsigned short method) :
  AbstractParameterAliasable("D1Walk."),
  AbstractReversibleSubstitutionModel(alpha, std::shared_ptr<StateMapInterface>(new CanonicalStateMap(alpha, false)), "D1Walk."),
  freqSet_(0)
{
  freqSet_ = std::make_shared<FullFrequencySet>(getStateMap(), true, method);
  freqSet_->setNamespace("D1Walk.");
  computeFrequencies(false);

  addParameters_(freqSet_->getParameters());
  // Exchangeability Matrix:
  for (unsigned int i = 0; i < size_-1; i++)
    exchangeability_(i, i+1) = 1;

  for (unsigned int i = 1; i < size_; i++)
    exchangeability_(i, i-1) = 1;
  
  updateMatrices_();
}

/******************************************************************************/

void D1WalkSubstitutionModel::updateMatrices_()
{
  // Frequencies:
  freq_ = freqSet_->getFrequencies();

  AbstractReversibleSubstitutionModel::updateMatrices_();
}


/******************************************************************************/

void D1WalkSubstitutionModel::setFreq(std::map<int, double>& freqs)
{
  for (auto i : freqs)
  {
    freq_[(size_t)i.first] = i.second;
  }

  freqSet_->setFrequencies(freq_);
  // Update parameters and re-compute generator and eigen values:
  matchParametersValues(freqSet_->getParameters());
}

/******************************************************************************/
