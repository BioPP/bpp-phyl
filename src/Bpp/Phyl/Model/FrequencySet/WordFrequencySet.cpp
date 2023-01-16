//
// File: WordFrequencySet.cpp
// Authors:
//   Laurent Gueguen
// Created: lundi 2 avril 2012, ÃÂ  14h 02
//

/*
  Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)
  
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


#include "WordFrequencySet.h"

using namespace bpp;

#include <cmath>
using namespace std;

size_t AbstractWordFrequencySet::getSizeFromVector(const std::vector<std::unique_ptr<FrequencySetInterface> >& freqVector)
{
  size_t s = 1;
  size_t l = freqVector.size();

  for (size_t i = 0; i < l; i++)
  {
    s *= freqVector[i]->getAlphabet()->getSize();
  }

  return s;
}

AbstractWordFrequencySet::AbstractWordFrequencySet(std::shared_ptr<const StateMapInterface> stateMap, const string& prefix, const string& name) :
  AbstractFrequencySet(stateMap, prefix, name)
{}

size_t AbstractWordFrequencySet::getLength() const
{
  return getWordAlphabet()->getLength();
}

AbstractWordFrequencySet::~AbstractWordFrequencySet()
{}

// ///////////////////////////////////////////////////////////////////
// // WordFromIndependentFrequencySet


WordFromIndependentFrequencySet::WordFromIndependentFrequencySet(
    shared_ptr<const WordAlphabet> pWA,
    vector<unique_ptr<FrequencySetInterface>>& freqVector,
    const string& prefix,
    const string& name) :
  AbstractWordFrequencySet(
    make_shared<CanonicalStateMap>(pWA, false),
    prefix,
    name),
  vFreq_(),
  vNestedPrefix_()
{
  size_t sf = getSizeFromVector(freqVector);
  if (pWA->getSize() !=  sf)
    throw Exception("WordFromIndependentFrequencySet: Size of the frequencies does not match size of the alphabet : " + TextTools::toString(sf) + " vs " + TextTools::toString(pWA->getSize()));

  size_t l = freqVector.size();

  for (size_t i = 0; i < l; ++i)
  {
    vFreq_.push_back(move(freqVector[i]));
    vNestedPrefix_.push_back(freqVector[i]->getNamespace());
    vFreq_[i]->setNamespace(prefix + TextTools::toString(i + 1) + "_" + vNestedPrefix_[i]);
    addParameters_(vFreq_[i]->getParameters());
  }

  updateFrequencies();
}

WordFromIndependentFrequencySet::WordFromIndependentFrequencySet(
    shared_ptr<const CodonAlphabet> pWA,
    vector<unique_ptr<FrequencySetInterface>>& freqVector,
    const string& prefix,
    const string& name) :
  AbstractWordFrequencySet(
    make_shared<CanonicalStateMap>(pWA, false),
    prefix,
    name),
  vFreq_(),
  vNestedPrefix_()
{
  size_t sf = getSizeFromVector(freqVector);
  if (pWA->getSize() !=  sf)
    throw Exception("WordFromIndependentFrequencySet: Size of the frequencies does not match size of the alphabet : " + TextTools::toString(sf) + " vs " + TextTools::toString(pWA->getSize()));

  size_t l = freqVector.size();

  for (size_t i = 0; i < l; ++i)
  {
    vFreq_.push_back(move(freqVector[i]));
    vNestedPrefix_.push_back(freqVector[i]->getNamespace());
    vFreq_[i]->setNamespace(prefix + TextTools::toString(i + 1) + "_" + vNestedPrefix_[i]);
    addParameters_(vFreq_[i]->getParameters());
  }

  updateFrequencies();
}

WordFromIndependentFrequencySet::WordFromIndependentFrequencySet(const WordFromIndependentFrequencySet& iwfs) :
  AbstractWordFrequencySet(iwfs),
  vFreq_(iwfs.vFreq_.size()),
  vNestedPrefix_(iwfs.vNestedPrefix_)
{
  for (unsigned i = 0; i < iwfs.vFreq_.size(); i++)
  {
    vFreq_[i].reset(iwfs.vFreq_[i]->clone());
  }
  updateFrequencies();
}

WordFromIndependentFrequencySet::~WordFromIndependentFrequencySet()
{}

WordFromIndependentFrequencySet& WordFromIndependentFrequencySet::operator=(const WordFromIndependentFrequencySet& iwfs)
{
  AbstractWordFrequencySet::operator=(iwfs);
  vNestedPrefix_ = iwfs.vNestedPrefix_;

  // Clean current frequencies first:
  vFreq_.resize(iwfs.vFreq_.size());
  for (unsigned i = 0; i < vFreq_.size(); i++)
  {
    vFreq_[i].reset(iwfs.vFreq_[i]->clone());
  }
  updateFrequencies();

  return *this;
}

void WordFromIndependentFrequencySet::fireParameterChanged(const ParameterList& pl)
{
  size_t l = vFreq_.size();

  bool f = 0;
  for (size_t i = 0; i < l; i++)
  {
    f |= vFreq_[i]->matchParametersValues(pl);
  }

  if (f)
    updateFrequencies();
}

void WordFromIndependentFrequencySet::updateFrequencies()
{
  size_t l = vFreq_.size();
  size_t s = getWordAlphabet()->getSize();
  vector< vector<double> > f(l);

  size_t i, p, t, i2;

  for (i = 0; i < l; i++)
  {
    f[i] = vFreq_[i]->getFrequencies();
  }

  for (i = 0; i < s; i++)
  {
    i2 = i;
    getFreq_(i) = 1;
    for (p = l; p > 0; p--)
    {
      t = vFreq_[p - 1]->getAlphabet()->getSize();
      getFreq_(i) *= f[p - 1][i2 % t];
      i2 /= t;
    }
  }
}

void WordFromIndependentFrequencySet::setFrequencies(const vector<double>& frequencies)
{
  if (frequencies.size() != getWordAlphabet()->getSize())
    throw DimensionException("WordFromIndependentFrequencySet::setFrequencies", frequencies.size(), getWordAlphabet()->getSize());
  double sum = 0.0;
  size_t size = frequencies.size();
  for (size_t i = 0; i < size; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.000001)
    throw Exception("WordFromIndependentFrequencySet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");

  size_t d, i, j, k, s, l = vFreq_.size();
  vector<double> freq;

  d = size;
  for (i = 0; i < l; i++)
  {
    s = vFreq_[i]->getAlphabet()->getSize();
    freq.resize(s);
    d /= s;
    for (j = 0; j < s; j++)
    {
      freq[j] = 0;
    }
    for (k = 0; k < size; k++)
    {
      freq[(k / d) % s] += frequencies[k];
    }
    vFreq_[i]->setFrequencies(freq);
  }

  for (i = 0; i < l; i++)
  {
    matchParametersValues(vFreq_[i]->getParameters());
  }

  updateFrequencies();
}


size_t WordFromIndependentFrequencySet::getLength() const
{
  return vFreq_.size();
}

void WordFromIndependentFrequencySet::setNamespace(const std::string& prefix)
{
  AbstractFrequencySet::setNamespace(prefix);
  for (size_t i = 0; i < vFreq_.size(); i++)
  {
    vFreq_[i]->setNamespace(prefix + TextTools::toString(i + 1) + "_" + vNestedPrefix_[i]);
  }
}

std::string WordFromIndependentFrequencySet::getDescription() const
{
  string s = getName() + " : " + vFreq_[0]->getName();
  for (size_t i = 1; i < vFreq_.size(); i++)
  {
    s += " * " + vFreq_[i]->getName();
  }
  return s;
}

// ///////////////////////////////////////////////////////////////////
// // WordFromUniqueFrequencySet


WordFromUniqueFrequencySet::WordFromUniqueFrequencySet(
    shared_ptr<const WordAlphabet> pWA,
    unique_ptr<FrequencySetInterface> pabsfreq,
    const string& prefix,
    const string& name) :
  AbstractWordFrequencySet(
    make_shared<CanonicalStateMap>(pWA, false),
    prefix,
    name),
  pFreq_(move(pabsfreq)),
  NestedPrefix_(pabsfreq->getNamespace()),
  length_(pWA->getLength())
{
  size_t i;

  string st = "";
  for (i = 0; i < length_; ++i)
  {
    st += TextTools::toString(i + 1);
  }

  pFreq_->setNamespace(prefix + st + "_" + NestedPrefix_);
  addParameters_(pFreq_->getParameters());

  updateFrequencies();
}

WordFromUniqueFrequencySet::WordFromUniqueFrequencySet(
    shared_ptr<const CodonAlphabet> pWA,
    unique_ptr<FrequencySetInterface> pabsfreq,
    const string& prefix,
    const string& name) :
  AbstractWordFrequencySet(
    make_shared<CanonicalStateMap>(pWA, false),
    prefix,
    name),
  pFreq_(move(pabsfreq)),
  NestedPrefix_(pabsfreq->getNamespace()),
  length_(pWA->getLength())
{
  size_t i;

  string st = "";
  for (i = 0; i < length_; ++i)
  {
    st += TextTools::toString(i + 1);
  }

  pFreq_->setNamespace(prefix + st + "_" + NestedPrefix_);
  addParameters_(pFreq_->getParameters());

  updateFrequencies();
}

WordFromUniqueFrequencySet::WordFromUniqueFrequencySet(const WordFromUniqueFrequencySet& iwfs) :
  AbstractWordFrequencySet(iwfs),
  pFreq_(iwfs.pFreq_->clone()),
  NestedPrefix_(iwfs.NestedPrefix_),
  length_(iwfs.length_)
{
  updateFrequencies();
}


WordFromUniqueFrequencySet& WordFromUniqueFrequencySet::operator=(const WordFromUniqueFrequencySet& iwfs)
{
  AbstractWordFrequencySet::operator=(iwfs);
  pFreq_.reset(iwfs.pFreq_->clone());
  NestedPrefix_ = iwfs.NestedPrefix_;
  length_ = iwfs.length_;

  updateFrequencies();
  return *this;
}

WordFromUniqueFrequencySet::~WordFromUniqueFrequencySet()
{
  pFreq_ = 0;
}

void WordFromUniqueFrequencySet::fireParameterChanged(const ParameterList& pl)
{
  if (pFreq_->matchParametersValues(pl))
    updateFrequencies();
}

void WordFromUniqueFrequencySet::updateFrequencies()
{
  size_t s = getWordAlphabet()->getSize();
  vector<double> f;
  size_t letsi = pFreq_->getAlphabet()->getSize();

  size_t i, p, i2;

  f = pFreq_->getFrequencies();

  for (i = 0; i < s; i++)
  {
    i2 = i;
    getFreq_(i2) = 1;
    for (p = length_; p > 0; p--)
    {
      getFreq_(i) *= f[i2 % letsi];
      i2 /= letsi;
    }
  }
}

void WordFromUniqueFrequencySet::setFrequencies(const vector<double>& frequencies)
{
  if (frequencies.size() != getWordAlphabet()->getSize())
    throw DimensionException("WordFromUniqueFrequencySet::setFrequencies", frequencies.size(), getWordAlphabet()->getSize());
  double sum = 0.0;
  size_t size = frequencies.size();
  for (size_t i = 0; i < size; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.000001)
    throw Exception("WordFromUniqueFrequencySet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");

  size_t d, i, j, k;
  vector<double> freq;

  size_t letsi = pFreq_->getAlphabet()->getSize();
  freq.resize(letsi);

  for (j = 0; j < letsi; j++)
  {
    freq[j] = 0;
  }

  d = size;
  for (i = 0; i < length_; i++)
  {
    d /= letsi;
    for (k = 0; k < size; k++)
    {
      freq[(k / d) % letsi] += frequencies[k];
    }
  }
  for (j = 0; j < letsi; j++)
  {
    freq[j] /= static_cast<double>(length_);
  }

  pFreq_->setFrequencies(freq);
  matchParametersValues(pFreq_->getParameters());
  updateFrequencies();
}


void WordFromUniqueFrequencySet::setNamespace(const string& prefix)
{
  AbstractFrequencySet::setNamespace(prefix);
  string st = "";
  for (unsigned i = 0; i < length_; i++)
  {
    st += TextTools::toString(i + 1);
  }
  pFreq_->setNamespace(prefix + st + "_" + NestedPrefix_);
}


string WordFromUniqueFrequencySet::getDescription() const
{
  return getName() + " : " + pFreq_->getName() + " * " + TextTools::toString(length_);
}
