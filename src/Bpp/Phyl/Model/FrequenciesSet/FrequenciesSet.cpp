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

using namespace bpp;

#include <cmath>
using namespace std;

IntervalConstraint FrequenciesSet::FREQUENCE_CONSTRAINT_MILLI(NumConstants::MILLI(), 1 - NumConstants::MILLI(), false, false);
IntervalConstraint FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL(NumConstants::SMALL(), 1 - NumConstants::SMALL(), false, false);

// ///////////////////////////////////////
// AbstractFrequenciesSet


void AbstractFrequenciesSet::setFrequenciesFromMap(const map<int, double>& frequencies)
{
  size_t s = getAlphabet()->getSize();
  vector<double> freq(s);
  double x = 0;
  for (size_t i = 0; i < s; i++)
  {
    map<int, double>::const_iterator it = frequencies.find(static_cast<int>(i));
    if (it != frequencies.end())
      freq[i] = it->second;
    else
      freq[i] = 0;
    x += freq[i];
  }
  for (size_t i = 0; i < s; i++)
  {
    freq[i] /= x;
  }
  setFrequencies(freq);
}

// ////////////////////////////
// FullFrequenciesSet

FullFrequenciesSet::FullFrequenciesSet(const Alphabet* alphabet, bool allowNullFreqs, const string& name) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Full.", name)
{
  size_t size = alphabet->getSize();

  for (size_t i = 0; i < alphabet->getSize() - 1; i++)
  {
    addParameter_(new Parameter(
      "Full.theta" + TextTools::toString(i + 1),
      1. / static_cast<double>(size - i),
      allowNullFreqs ?
      &Parameter::PROP_CONSTRAINT_IN :
      &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
    getFreq_(i) = 1. / static_cast<double>(size);
  }
  size_t i = alphabet->getSize() - 1;
  getFreq_(i) = 1. / static_cast<double>(size);
}

FullFrequenciesSet::FullFrequenciesSet(const Alphabet* alphabet, const vector<double>& initFreqs, bool allowNullFreqs, const string& name) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Full.", name)
{
  if (initFreqs.size() != alphabet->getSize())
    throw Exception("FullFrequenciesSet(constructor). There must be " + TextTools::toString(alphabet->getSize()) + " frequencies.");
  double sum = VectorTools::sum(initFreqs);
  if (fabs(1. - sum) > NumConstants::SMALL())
  {
    throw Exception("Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  }

  double y = 1;
  for (size_t i = 0; i < alphabet->getSize() - 1; i++)
  {
    addParameter_(new Parameter(
      "Full.theta" + TextTools::toString(i + 1),
      initFreqs[i] / y,
      allowNullFreqs ?
      &Parameter::PROP_CONSTRAINT_IN :
      &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
    getFreq_(i) = initFreqs[i];
    y -= initFreqs[i];
  }
  size_t i = alphabet->getSize() - 1;
  getFreq_(i) = initFreqs[i];
}

void FullFrequenciesSet::setFrequencies(const vector<double>& frequencies) 
{
  const Alphabet* alphabet = getAlphabet();
  if (frequencies.size() != getNumberOfFrequencies())
    throw DimensionException("FullFrequenciesSet::setFrequencies. Invalid number of frequencies.", frequencies.size(), getNumberOfFrequencies());

  if (fabs(1. - VectorTools::sum(frequencies)) >= NumConstants::SMALL())
    throw Exception("FullFrequenciesSet::setFrequencies. Frequencies do not sum to 1 : " + TextTools::toString(VectorTools::sum(frequencies)));

  setFrequencies_(frequencies);

  double y = 1;
  for (size_t i = 0; i < alphabet->getSize() - 1; i++)
  {
    getParameter_("theta" + TextTools::toString(i + 1)).setValue(frequencies[i] / y);
    y -= frequencies[i];
  }
}

void FullFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
  const Alphabet* alphabet = getAlphabet();
  double y = 1;
  size_t i;
  for (i = 0; i < alphabet->getSize() - 1; i++)
  {
    getFreq_(i) = getParameter_("theta" + TextTools::toString(i + 1)).getValue() * y;
    y *= 1 - getParameter_("theta" + TextTools::toString(i + 1)).getValue();
  }

  i = alphabet->getSize() - 1;
  getFreq_(i) = y;
}


// ///////////////////////////////////////////
// / FixedFrequenciesSet

FixedFrequenciesSet::FixedFrequenciesSet(const Alphabet* alphabet, const vector<double>& initFreqs, const string& name) :
  AbstractFrequenciesSet(initFreqs.size(), alphabet, "Fixed.", name)
{
  setFrequencies(initFreqs);
}

FixedFrequenciesSet::FixedFrequenciesSet(const Alphabet* alphabet, size_t nFreqs, const string& name) :
  AbstractFrequenciesSet(nFreqs, alphabet, "Fixed.", name)
{
  size_t size = nFreqs;
  for (size_t i = 0; i < nFreqs; ++i)
  {
    getFreq_(i) = 1. / static_cast<double>(size);
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
  AbstractFrequenciesSet(freqSet->getAlphabet()->getSize() * rateFreqs.size(), freqSet->getAlphabet(), "MarkovModulated.", "MarkovModulated." + freqSet->getName()),
  freqSet_(freqSet),
  rateFreqs_(rateFreqs)
{
  freqSet_->setNamespace(std::string("MarkovModulated.") + freqSet_->getNamespace());
  addParameters_(freqSet_->getParameters());
  setFrequencies_(VectorTools::kroneckerMult(rateFreqs, freqSet_->getFrequencies()));
}

