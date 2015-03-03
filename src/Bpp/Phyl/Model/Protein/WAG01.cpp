//
// File: WAG01.cpp
// Created by: Laurent Gueguen
// Created on: mardi 28 septembre 2010, à 14h 43
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

#include "WAG01.h"

//From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

WAG01::WAG01(const ProteicAlphabet* alpha) :
  AbstractParameterAliasable("WAG01."),
  AbstractReversibleProteinSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "WAG01."),
  freqSet_(0)
{
  #include "__WAG01ExchangeabilityCode"
  #include "__WAG01FrequenciesCode"
  freqSet_ = new FixedProteinFrequenciesSet(alpha, freq_);
  updateMatrices();  
}

WAG01::WAG01(const ProteicAlphabet* alpha, ProteinFrequenciesSet* freqSet, bool initFreqs) :
  AbstractParameterAliasable("WAG01+F."),
  AbstractReversibleProteinSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "WAG01+F."),
  freqSet_(freqSet)
{
  #include "__WAG01ExchangeabilityCode"
  #include "__WAG01FrequenciesCode"
  freqSet_->setNamespace("WAG01+F."+freqSet_->getNamespace());
  if (initFreqs) freqSet_->setFrequencies(freq_);
  else freq_ = freqSet_->getFrequencies();
  addParameters_(freqSet_->getParameters());
  updateMatrices();  
}

/******************************************************************************/

void WAG01::setFreqFromData(const SequenceContainer& data, double pseudoCount)
{
  map<int, int> counts;
  SequenceContainerTools::getCounts(data, counts);
  double t = 0;
  for (int i = 0; i < static_cast<int>(size_); i++)
  {
    t += (counts[i] + pseudoCount);
  }
  for (size_t i = 0; i < size_; ++i) freq_[i] = (static_cast<double>(counts[static_cast<int>(i)]) + pseudoCount) / t;
  freqSet_->setFrequencies(freq_);
  //Update parameters and re-compute generator and eigen values:
  matchParametersValues(freqSet_->getParameters());
}

/******************************************************************************/

