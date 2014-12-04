//
// File: WordFrequenciesSet.cpp
// Created by: Laurent Gueguen
// Created on: lundi 2 avril 2012, à 14h 02
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

#include "WordFrequenciesSet.h"


using namespace bpp;

#include <cmath>
using namespace std;

size_t AbstractWordFrequenciesSet::getSizeFromVector(const std::vector<FrequenciesSet*>& freqVector)
{
  size_t s = 1;
  size_t l = freqVector.size();

  for (size_t i = 0; i < l; i++)
  {
    s *= freqVector[i]->getAlphabet()->getSize();
  }

  return s;
}

AbstractWordFrequenciesSet::AbstractWordFrequenciesSet(StateMap* stateMap, const string& prefix, const string& name) :
  AbstractFrequenciesSet(stateMap, prefix, name)
{}

size_t AbstractWordFrequenciesSet::getLength() const
{
  return dynamic_cast<const WordAlphabet*>(getAlphabet())->getLength();
}

AbstractWordFrequenciesSet::~AbstractWordFrequenciesSet()
{}

// ///////////////////////////////////////////////////////////////////
// // WordFromIndependentFrequenciesSet


WordFromIndependentFrequenciesSet::WordFromIndependentFrequenciesSet(
    const WordAlphabet* pWA,
    const std::vector<FrequenciesSet*>& freqVector,
    const string& prefix, const string& name) :
  AbstractWordFrequenciesSet(new CanonicalStateMap(pWA, false), prefix, name),
  vFreq_(),
  vNestedPrefix_()
{
  size_t sf = getSizeFromVector(freqVector);
  if (pWA->getSize() !=  sf)
    throw Exception("WordFromIndependentFrequenciesSet: Size of the frequencies does not match size of the alphabet : " + TextTools::toString(sf) + " vs " + TextTools::toString(pWA->getSize()));

  size_t l = freqVector.size();

  for (size_t i = 0; i < l; i++)
  {
    vFreq_.push_back(freqVector[i]);
    vNestedPrefix_.push_back(freqVector[i]->getNamespace());
    vFreq_[i]->setNamespace(prefix + TextTools::toString(i + 1) + "_" + vNestedPrefix_[i]);
    addParameters_(vFreq_[i]->getParameters());
  }

  updateFrequencies();
}

WordFromIndependentFrequenciesSet::WordFromIndependentFrequenciesSet(const WordFromIndependentFrequenciesSet& iwfs) :
  AbstractWordFrequenciesSet(iwfs),
  vFreq_(iwfs.vFreq_.size()),
  vNestedPrefix_(iwfs.vNestedPrefix_)
{
  for (unsigned i = 0; i < iwfs.vFreq_.size(); i++)
  {
    vFreq_[i] =  iwfs.vFreq_[i]->clone();
  }
  updateFrequencies();
}

WordFromIndependentFrequenciesSet::~WordFromIndependentFrequenciesSet()
{
  for (unsigned i = 0; i < vFreq_.size(); i++)
  {
    delete vFreq_[i];
  }
}

WordFromIndependentFrequenciesSet& WordFromIndependentFrequenciesSet::operator=(const WordFromIndependentFrequenciesSet& iwfs)
{
  AbstractWordFrequenciesSet::operator=(iwfs);
  vNestedPrefix_ = iwfs.vNestedPrefix_;

  //Clean current frequencies first:
  for (unsigned i = 0; i < vFreq_.size(); i++)
  {
    delete vFreq_[i];
  }

  vFreq_.resize(iwfs.vFreq_.size());
  for (unsigned i = 0; i < vFreq_.size(); i++)
  {
    vFreq_[i] = iwfs.vFreq_[i]->clone();
  }
  updateFrequencies();

  return *this;
}

void WordFromIndependentFrequenciesSet::fireParameterChanged(const ParameterList& pl)
{
  AbstractFrequenciesSet::fireParameterChanged(pl);
  size_t l = vFreq_.size();

  bool f = 0;
  for (size_t i = 0; i < l; i++)
  {
    f |= vFreq_[i]->matchParametersValues(pl);
  }

  if (f)
    updateFrequencies();
}

void WordFromIndependentFrequenciesSet::updateFrequencies()
{
  size_t l = vFreq_.size();
  size_t s = getAlphabet()->getSize();
  vector< vector<double> >f(l);

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

void WordFromIndependentFrequenciesSet::setFrequencies(const vector<double>& frequencies) 
{
  if (frequencies.size() != getAlphabet()->getSize())
    throw DimensionException("WordFromIndependentFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  size_t size = frequencies.size();
  for (size_t i = 0; i < size; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.000001)
    throw Exception("WordFromIndependentFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");

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


size_t WordFromIndependentFrequenciesSet::getLength() const
{
  return vFreq_.size();
}

void WordFromIndependentFrequenciesSet::setNamespace(const std::string& prefix)
{
  AbstractFrequenciesSet::setNamespace(prefix);
  for (size_t i = 0; i < vFreq_.size(); i++)
  {
    vFreq_[i]->setNamespace(prefix + TextTools::toString(i + 1) + "_" + vNestedPrefix_[i]);
  }
}

std::string WordFromIndependentFrequenciesSet::getDescription() const
{
  string s = getName() +" : " + vFreq_[0]->getName();
  for (size_t i = 1; i < vFreq_.size(); i++)
  {
    s += " * " + vFreq_[i]->getName();
  }
  return s;
}

// ///////////////////////////////////////////////////////////////////
// // WordFromUniqueFrequenciesSet


WordFromUniqueFrequenciesSet::WordFromUniqueFrequenciesSet(
    const WordAlphabet* pWA,
    FrequenciesSet* pabsfreq,
    const string& prefix,
    const string& name) :
  AbstractWordFrequenciesSet(new CanonicalStateMap(pWA, false), prefix, name),
  pFreq_(pabsfreq),
  NestedPrefix_(pabsfreq->getNamespace()),
  length_(pWA->getLength())
{
  size_t i;

  string st = "";
  for (i = 0; i < length_; i++)
  {
    st += TextTools::toString(i + 1);
  }

  pFreq_->setNamespace(prefix+ st + "_" + NestedPrefix_);
  addParameters_(pFreq_->getParameters());

  updateFrequencies();
}

WordFromUniqueFrequenciesSet::WordFromUniqueFrequenciesSet(const WordFromUniqueFrequenciesSet& iwfs) :
  AbstractWordFrequenciesSet(iwfs),
  pFreq_(iwfs.pFreq_->clone()),
  NestedPrefix_(iwfs.NestedPrefix_),
  length_(iwfs.length_)
{
  updateFrequencies();
}


WordFromUniqueFrequenciesSet& WordFromUniqueFrequenciesSet::operator=(const WordFromUniqueFrequenciesSet& iwfs)
{
  AbstractWordFrequenciesSet::operator=(iwfs);
  delete pFreq_;
  pFreq_ = iwfs.pFreq_->clone();
  NestedPrefix_ = iwfs.NestedPrefix_;
  length_ = iwfs.length_;

  updateFrequencies();
  return *this;
}

WordFromUniqueFrequenciesSet::~WordFromUniqueFrequenciesSet()
{
  if (pFreq_)
    delete pFreq_;
  pFreq_ = 0;
}

void WordFromUniqueFrequenciesSet::fireParameterChanged(const ParameterList& pl)
{
  AbstractFrequenciesSet::fireParameterChanged(pl);
  if (pFreq_->matchParametersValues(pl))
    updateFrequencies();
}

void WordFromUniqueFrequenciesSet::updateFrequencies()
{
  size_t s = getAlphabet()->getSize();
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

void WordFromUniqueFrequenciesSet::setFrequencies(const vector<double>& frequencies) 
{
  if (frequencies.size() != getAlphabet()->getSize())
    throw DimensionException("WordFromUniqueFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  size_t size = frequencies.size();
  for (size_t i = 0; i < size; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.000001)
    throw Exception("WordFromUniqueFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");

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


void WordFromUniqueFrequenciesSet::setNamespace(const string& prefix)
{
  AbstractFrequenciesSet::setNamespace(prefix);
  string st = "";
  for (unsigned i = 0; i < length_; i++)
  {
    st += TextTools::toString(i + 1);
  }
  pFreq_->setNamespace(prefix + st + "_" + NestedPrefix_);
}


string WordFromUniqueFrequenciesSet::getDescription() const
{
  return getName() + " : " + pFreq_->getName() + " * " + TextTools::toString(length_);
}


