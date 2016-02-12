//
// File: FrequenciesSet.cpp
// Created by: Bastien Boussau
//             Julien Dutheil
// Created on: Tue Aug 21 2007
//

/*
   Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "FrequenciesSet.h"

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/Simplex.h>

// From bpp-phyl
#include "../SubstitutionModel.h"
#include "../AbstractBiblioSubstitutionModel.h"


using namespace bpp;

#include <cmath>
using namespace std;

IntervalConstraint FrequenciesSet::FREQUENCE_CONSTRAINT_MILLI(NumConstants::MILLI(), 1 - NumConstants::MILLI(), false, false);
IntervalConstraint FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL(NumConstants::SMALL(), 1 - NumConstants::SMALL(), false, false);

// ///////////////////////////////////////
// AbstractFrequenciesSet


void AbstractFrequenciesSet::setFrequenciesFromAlphabetStatesFrequencies(const map<int, double>& frequencies)
{
  size_t s = stateMap_->getNumberOfModelStates();
  vector<double> freq(s);
  double x = 0;
  for (size_t i = 0; i < s; ++i)
  {
    map<int, double>::const_iterator it = frequencies.find(stateMap_->getAlphabetStateAsInt(i));
    if (it != frequencies.end())
      freq[i] = it->second;
    else
      freq[i] = 0;
    x += freq[i];
  }

  if (x!=0)
    for (size_t i = 0; i < s; ++i)
      freq[i] /= x;
  else
    for (size_t i = 0; i < s; ++i)
      freq[i] = 1.0/(double)s;

  setFrequencies(freq);
}

const std::map<int, double> AbstractFrequenciesSet::getAlphabetStatesFrequencies() const
{
  map<int, double> fmap;
  for (size_t i = 0; i < stateMap_->getNumberOfModelStates(); ++i) {
    fmap[stateMap_->getAlphabetStateAsInt(i)] += freq_[i];
  }
  return fmap;
}

// ////////////////////////////
// FullFrequenciesSet

FullFrequenciesSet::FullFrequenciesSet(StateMap* stateMap, bool allowNullFreqs, unsigned short method, const string& name) :
  AbstractFrequenciesSet(stateMap, "Full.", name),
  sFreq_(stateMap->getNumberOfModelStates(), method, allowNullFreqs, "Full.")
{
  vector<double> vd;
  double r = 1. / static_cast<double>(stateMap->getNumberOfModelStates());

  for (size_t i = 0; i < stateMap->getNumberOfModelStates(); i++)
    vd.push_back(r);

  sFreq_.setFrequencies(vd);
  addParameters_(sFreq_.getParameters());
  updateFreq_();
}

FullFrequenciesSet::FullFrequenciesSet(StateMap* stateMap, const vector<double>& initFreqs, bool allowNullFreqs, unsigned short method, const string& name) :
  AbstractFrequenciesSet(stateMap, "Full.", name),
  sFreq_(stateMap->getNumberOfModelStates(), method, allowNullFreqs, "Full.")
{
  sFreq_.setFrequencies(initFreqs);
  addParameters_(sFreq_.getParameters());
  updateFreq_();
}

void FullFrequenciesSet::setFrequencies(const vector<double>& frequencies) 
{
  sFreq_.setFrequencies(frequencies);
  setParametersValues(sFreq_.getParameters()); 

  updateFreq_();
}

void FullFrequenciesSet::setNamespace(const std::string& nameSpace)
{
  sFreq_.setNamespace(nameSpace);
  AbstractFrequenciesSet::setNamespace(nameSpace);
}

void FullFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
  AbstractFrequenciesSet::fireParameterChanged(parameters);
  sFreq_.matchParametersValues(parameters);
  updateFreq_();
}

void FullFrequenciesSet::updateFreq_()
{
  for (size_t i = 0; i < getAlphabet()->getSize(); i++)
    getFreq_(i)=sFreq_.prob(i);
}

// ///////////////////////////////////////////
// / FixedFrequenciesSet

FixedFrequenciesSet::FixedFrequenciesSet(StateMap* stateMap, const vector<double>& initFreqs, const string& name) throw (Exception) :
  AbstractFrequenciesSet(stateMap, "Fixed.", name)
{
  if (stateMap->getNumberOfModelStates() != initFreqs.size())
    throw Exception("FixedFrequenciesSet::constructor. size of init vector does not match the number of states in the model.");
  setFrequencies(initFreqs);
}

FixedFrequenciesSet::FixedFrequenciesSet(StateMap* stateMap, const string& name) :
  AbstractFrequenciesSet(stateMap, "Fixed.", name)
{
  size_t n = stateMap->getNumberOfModelStates();
  for (size_t i = 0; i < n; ++i)
  {
    getFreq_(i) = 1. / static_cast<double>(n);
  }
}

void FixedFrequenciesSet::setFrequencies(const vector<double>& frequencies) 
{
  if (frequencies.size() != getNumberOfFrequencies())
    throw DimensionException("FixedFrequenciesSet::setFrequencies", frequencies.size(), getNumberOfFrequencies());
  double sum = 0.0;
  for (size_t i = 0; i < frequencies.size(); i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.00001)
    throw Exception("FixedFrequenciesSet::setFrequencies. Frequencies sum must equal 1 (sum = " + TextTools::toString(sum) + ").");
  setFrequencies_(frequencies);
}

MarkovModulatedFrequenciesSet::MarkovModulatedFrequenciesSet(FrequenciesSet* freqSet, const std::vector<double>& rateFreqs) :
  AbstractFrequenciesSet(new MarkovModulatedStateMap(freqSet->getStateMap(), static_cast<unsigned int>(rateFreqs.size())), "MarkovModulated.", "MarkovModulated." + freqSet->getName()),
  freqSet_(freqSet),
  rateFreqs_(rateFreqs)
{
  freqSet_->setNamespace(std::string("MarkovModulated.") + freqSet_->getNamespace());
  addParameters_(freqSet_->getParameters());
  setFrequencies_(VectorTools::kroneckerMult(rateFreqs, freqSet_->getFrequencies()));
}



//////////////////////////////////
/// From Model


FromModelFrequenciesSet::FromModelFrequenciesSet(const FromModelFrequenciesSet& fmfs):
  AbstractFrequenciesSet(fmfs),
  model_(fmfs.model_->clone())
{}

FromModelFrequenciesSet& FromModelFrequenciesSet::operator=(const FromModelFrequenciesSet& fmfs) 
{
  AbstractFrequenciesSet::operator=(fmfs);
  model_ = fmfs.model_->clone();
  return *this;
}

FromModelFrequenciesSet::~FromModelFrequenciesSet()
{
  delete model_;
}

FromModelFrequenciesSet::FromModelFrequenciesSet(SubstitutionModel* model) :
  AbstractFrequenciesSet(model->getStateMap().clone(), "FromModel."+(model?model->getNamespace():""), "FromModel"),
  model_(model)
{
  model_->setNamespace(getNamespace());
  addParameters_(model_->getParameters());
  setFrequencies_(model_->getFrequencies());
}


void FromModelFrequenciesSet::setNamespace(const std::string& name)
{
  AbstractParameterAliasable::setNamespace(name);
  model_->setNamespace(name);
}


void FromModelFrequenciesSet::setFrequencies(const std::vector<double>& frequencies)
{
  std::map<int, double> freq;
  for (size_t i = 0; i < getNumberOfFrequencies(); ++i) {
    freq[getStateMap().getAlphabetStateAsInt(i)] += frequencies[i];
  }
  model_->setFreq(freq);
  matchParametersValues(model_->getParameters());
}

void FromModelFrequenciesSet::fireParameterChanged(const ParameterList& pl)
{
  AbstractFrequenciesSet::fireParameterChanged(pl);
  model_->matchParametersValues(pl);
  setFrequencies_(model_->getFrequencies());
}
