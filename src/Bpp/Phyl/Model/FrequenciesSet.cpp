//
// File: FrequenciesSet.cpp
// Created by: Bastien Boussau
//             Julien Dutheil
// Created on: Tue Aug 21 2007
//

/*
   Copyright or (c) or Copr. CNRS, (November 16, 2004)

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

#include <Seq/CodonAlphabet.h>
#include <Seq/AlphabetTools.h>

using namespace bpp;

#include <cmath>
using namespace std;

// From NumCalc:
#include <NumCalc/NumConstants.h>
#include <NumCalc/VectorTools.h>

// ///////////////////////////////////////
// AbstractFrequenciesSet


void AbstractFrequenciesSet::setFrequenciesFromMap(const map<int, double>& frequencies)
{
   unsigned int s = getAlphabet()->getSize();
   vector<double> freq(s);
   double x = 0;
  for (unsigned int i = 0; i < s; i++)
  {
   map<int, double>::const_iterator it = frequencies.find(i);
    if (it != frequencies.end())
      freq[i] = it->second;
    else
      freq[i] = 0;
    x += freq[i];
  }
  for (unsigned int i = 0; i < s; i++)
  {
    freq[i] /= x;
  }
  setFrequencies(freq);
}

// ////////////////////////////
// FullFrequenciesSet

FullFrequenciesSet::FullFrequenciesSet(const Alphabet* alphabet, bool allowNullFreqs) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Full.")
{
   unsigned int size = alphabet->getSize();

  for (unsigned int i = 0; i < alphabet->getSize() - 1; i++)
  {
    Parameter p(
        "Full.theta" + TextTools::toString(i + 1),
        1. / (size - i),
        allowNullFreqs ?
        dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
        dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
    addParameter_(p);
    getFreq_(i) = 1. / size;
  }
  unsigned int i = alphabet->getSize() - 1;
  getFreq_(i) = 1. / size;
}

FullFrequenciesSet::FullFrequenciesSet(const Alphabet* alphabet, const vector<double>& initFreqs, bool allowNullFreqs) throw (Exception) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Full.")
{
  if (initFreqs.size() != alphabet->getSize())
    throw Exception("FullFrequenciesSet(constructor). There must be " + TextTools::toString(alphabet->getSize()) + " frequencies.");
  double sum = VectorTools::sum(initFreqs);
  if (fabs(1. - sum) > NumConstants::SMALL)
  {
    throw Exception("Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  }

  double y = 1;
  for (unsigned int i = 0; i < alphabet->getSize() - 1; i++)
  {
    Parameter p(
        "Full.theta" + TextTools::toString(i + 1),
        initFreqs[i] / y,
        allowNullFreqs ?
        dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
        dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
    addParameter_(p);
    getFreq_(i) = initFreqs[i];
    y -= initFreqs[i];
  }
  unsigned int i = alphabet->getSize() - 1;
  getFreq_(i) = initFreqs[i];
}

void FullFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
   const Alphabet* alphabet = getAlphabet();
  if (frequencies.size() != getNumberOfFrequencies())
    throw DimensionException("FullFrequenciesSet::setFrequencies. Invalid number of frequencies.", frequencies.size(), getNumberOfFrequencies());

  if (fabs(1. - VectorTools::sum(frequencies)) >= NumConstants::SMALL)
    throw Exception("FullFrequenciesSet::setFrequencies. Frequencies do not sum to 1 : " + TextTools::toString(VectorTools::sum(frequencies)));

  setFrequencies_(frequencies);

  double y = 1;
  for (unsigned int i = 0; i < alphabet->getSize() - 1; i++)
  {
    getParameter_("theta" + TextTools::toString(i + 1)).setValue(frequencies[i] / y);
    y -= frequencies[i];
  }
}

void FullFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
   const Alphabet* alphabet = getAlphabet();
   double y = 1;
   unsigned int i;
  for (i = 0; i < alphabet->getSize() - 1; i++)
  {
    getFreq_(i) = getParameter_("theta" + TextTools::toString(i + 1)).getValue() * y;
    y *= 1 - getParameter_("theta" + TextTools::toString(i + 1)).getValue();
  }

  i = alphabet->getSize() - 1;
  getFreq_(i) = y;
}


// ////////////////////////////
// FullCodonFrequenciesSet

FullCodonFrequenciesSet::FullCodonFrequenciesSet(const CodonAlphabet* alphabet, bool allowNullFreqs) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Full.")
{
   unsigned int size = alphabet->getSize() - alphabet->numberOfStopCodons();
   unsigned int j = 0;

  for (unsigned int i = 0; i < alphabet->getSize() - 1; i++)
  {
    if (alphabet->isStop(i))
    {
      getFreq_(i) = 0;
    }
    else
    {
      Parameter p(
          "Full.theta" + TextTools::toString(i + 1),
          1. / (size - j),
          allowNullFreqs ?
          dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
          dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
      addParameter_(p);
      getFreq_(i) = 1. / size;
      j++;
    }
  }
  unsigned int i = alphabet->getSize() - 1;
  getFreq_(i) = (alphabet->isStop(i)) ? 0 : 1. / size;
}


FullCodonFrequenciesSet::FullCodonFrequenciesSet(const CodonAlphabet* alphabet, const vector<double>& initFreqs, bool allowNullFreqs) throw (Exception) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Full.")
{
  if (initFreqs.size() != alphabet->getSize())
    throw Exception("FullCodonFrequenciesSet(constructor). There must be " + TextTools::toString(alphabet->getSize()) + " frequencies.");
  double sum = 0.0;

  for (unsigned int i = 0; i < initFreqs.size(); i++)
  {
    if (!alphabet->isStop(i))
    {
      sum += initFreqs[i];
    }
  }

  if (fabs(1. - sum) > NumConstants::SMALL)
  {
    throw Exception("Non stop frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  }

  double y = 1;
  for (unsigned int i = 0; i < alphabet->getSize() - 1; i++)
  {
    if (alphabet->isStop(i))
    {
      getFreq_(i) = 0;
    }
    else
    {
      Parameter p(
          "Full.theta" + TextTools::toString(i + 1),
          initFreqs[i] / y,
          allowNullFreqs ?
          dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
          dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
      addParameter_(p);
      getFreq_(i) = initFreqs[i];
      y -= initFreqs[i];
    }
  }
  unsigned int i = alphabet->getSize() - 1;
  getFreq_(i) = (alphabet->isStop(i)) ? 0 : initFreqs[i];
}

void FullCodonFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if (frequencies.size() != getAlphabet()->getSize()) throw DimensionException("FullFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  const CodonAlphabet* alphabet = getAlphabet();

  double sum = 0.0;
  unsigned int i;
  for (i = 0; i < frequencies.size(); i++)
  {
    if (!(alphabet->isStop(i)))
      sum += frequencies[i];
  }
  if (fabs(1. - sum) > NumConstants::SMALL)
  {
    throw Exception("FullFrequenciesSet::setFrequencies. Non stop frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  }

  double y = 1;
  for (i = 0; i < alphabet->getSize() - 1; i++)
  {
    if (alphabet->isStop(i))
    {
      getFreq_(i) = 0;
    }
    else
    {
      getParameter_("theta" + TextTools::toString(i + 1)).setValue(frequencies[i] / y);
      y -= frequencies[i];
      getFreq_(i) = frequencies[i];
    }
  }
  i = alphabet->getSize() - 1;
  getFreq_(i) = (alphabet->isStop(i)) ? 0 : frequencies[i];
}

void FullCodonFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
   const CodonAlphabet* alphabet = getAlphabet();
   double y = 1;
   unsigned int i;
  for (i = 0; i < alphabet->getSize() - 1; i++)
  {
    if (!(alphabet->isStop(i))){
      getFreq_(i) = getParameter_("theta" + TextTools::toString(i + 1)).getValue() * y;
      y *= 1 - getParameter_("theta" + TextTools::toString(i + 1)).getValue();
    }
  }

  i = alphabet->getSize() - 1;
  getFreq_(i) = (alphabet->isStop(i)) ? 0 : y;
}


// ///////////////////////////////////////
// FullNucleotideFrequenciesSet


FullNucleotideFrequenciesSet::FullNucleotideFrequenciesSet(const NucleicAlphabet* alphabet, bool allowNullFreqs) :
  AbstractFrequenciesSet(4, alphabet, "Full.")
{
  Parameter thetaP(
    "Full.theta", 0.5,
    allowNullFreqs ?
    dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
    dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
  addParameter_(thetaP);
  Parameter theta1P(
      "Full.theta1", 0.5,
      allowNullFreqs ?
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
  addParameter_(theta1P);
  Parameter theta2P("Full.theta2", 0.5,
      allowNullFreqs ?
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
  addParameter_(theta2P);
  getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
}

FullNucleotideFrequenciesSet::FullNucleotideFrequenciesSet(const NucleicAlphabet* alphabet, double theta, double theta1, double theta2, bool allowNullFreqs) :
  AbstractFrequenciesSet(4, alphabet, "Full.")
{
  Parameter thetaP(
    "Full.theta",
    theta,
    allowNullFreqs ?
    dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
    dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
  addParameter_(thetaP);
  Parameter theta1P(
      "Full.theta1",
      theta1,
      allowNullFreqs ?
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
  addParameter_(theta1P);
  Parameter theta2P(
      "Full.theta2",
      theta2,
      allowNullFreqs ?
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
  addParameter_(theta2P);
  getFreq_(0) = theta1 * (1. - theta);
  getFreq_(1) = (1 - theta2) * theta;
  getFreq_(2) = theta2 * theta;
  getFreq_(3) = (1 - theta1) * (1. - theta);
}

void FullNucleotideFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if (frequencies.size() != 4) throw DimensionException(" FullNucleotideFrequenciesSet::setFrequencies", frequencies.size(), 4);
  double sum = 0.0;
  for (unsigned int i = 0; i < 4; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > NumConstants::SMALL)
    throw Exception("FullNucleotideFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  double theta = frequencies[1] + frequencies[2];
  getParameter_(0).setValue(theta);
  getParameter_(1).setValue(frequencies[0] / (1 - theta));
  getParameter_(2).setValue(frequencies[2] / theta);

  setFrequencies_(frequencies);
}

void FullNucleotideFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
  double theta  = getParameter_(0).getValue();
  double theta1 = getParameter_(1).getValue();
  double theta2 = getParameter_(2).getValue();
  getFreq_(0) = theta1 * (1. - theta);
  getFreq_(1) = (1 - theta2) * theta;
  getFreq_(2) = theta2 * theta;
  getFreq_(3) = (1 - theta1) * (1. - theta);
}

// /////////////////////////////////////////
// GCFrequenciesSet

void GCFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if (frequencies.size() != 4) throw DimensionException("GCFrequenciesSet::setFrequencies", frequencies.size(), 4);
  double sum = 0.0;
  for (unsigned int i = 0; i < 4; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > NumConstants::SMALL)
    throw Exception("GCFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  double theta = frequencies[1] + frequencies[2];
  // We set everything in one shot here:
  getParameter_(0).setValue(theta);
  getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
  getFreq_(1) = getFreq_(2) = theta / 2.;
}

void GCFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
   double theta = getParameter_(0).getValue();
  getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
  getFreq_(1) = getFreq_(2) = theta / 2.;
}


// ///////////////////////////////////////////
// / FixedFrequenciesSet

FixedFrequenciesSet::FixedFrequenciesSet(const Alphabet* alphabet, const vector<double>& initFreqs) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.")
{
  setFrequencies(initFreqs);
}

FixedFrequenciesSet::FixedFrequenciesSet(const Alphabet* alphabet) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.")
{
   unsigned int size = alphabet->getSize();

  for (unsigned int i = 0; i < alphabet->getSize(); i++)
  {
    getFreq_(i) = 1. / size;
  }
}

void FixedFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if (frequencies.size() != getAlphabet()->getSize()) throw DimensionException("FixedFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  for (unsigned int i = 0; i < frequencies.size(); i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.00001)
    throw Exception("FixedFrequenciesSet::setFrequencies. Frequencies sum must equal 1 (sum = " + TextTools::toString(sum) + ").");
  setFrequencies_(frequencies);
}

// ///////////////////////////////////////////
// / FixedCodonFrequenciesSet

FixedCodonFrequenciesSet::FixedCodonFrequenciesSet(const CodonAlphabet* alphabet, const vector<double>& initFreqs) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.")
{
  setFrequencies(initFreqs);
}

FixedCodonFrequenciesSet::FixedCodonFrequenciesSet(const CodonAlphabet* alphabet) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.")
{
   unsigned int size = alphabet->getSize() - alphabet->numberOfStopCodons();

  for (unsigned int i = 0; i < alphabet->getSize(); i++)
  {
    getFreq_(i) = (alphabet->isStop(i)) ? 0 : 1. / size;
  }
}

void FixedCodonFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if (frequencies.size() != getAlphabet()->getSize()) throw DimensionException("FixedFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  for (unsigned int i = 0; i < frequencies.size(); i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.000001)
    throw Exception("FixedCodonFrequenciesSet::setFrequencies. Frequencies sum must equal 1 (sum = " + TextTools::toString(sum) + ").");
  setFrequencies_(frequencies);
}


// /////////////////////////////////////////////
// /////////////////////////////////////////////
// WordFrequenciesSet

unsigned int WordFrequenciesSet::getSizeFromVector(const std::vector<FrequenciesSet*>& freqVector)
{
  unsigned int s = 1;
  unsigned int l = freqVector.size();

  for (unsigned int i = 0; i < l; i++)
  {
    s *= freqVector[i]->getAlphabet()->getSize();
  }

  return s;
}

WordFrequenciesSet::WordFrequenciesSet(unsigned int size, const Alphabet* palph) :
  AbstractFrequenciesSet(size, palph, "Word.")
{}

unsigned int WordFrequenciesSet::getLength() const
{
  return dynamic_cast<const WordAlphabet*>(getAlphabet())->getLength();
}

WordFrequenciesSet::~WordFrequenciesSet()
{}


// ///////////////////////////////////////////////////////////////////
// // WordFromIndependentFrequenciesSet


WordFromIndependentFrequenciesSet::WordFromIndependentFrequenciesSet(
  const WordAlphabet* pWA, const std::vector<FrequenciesSet*>& freqvector) throw (Exception) :
  WordFrequenciesSet(pWA->getSize(), pWA),
  vFreq_(),
  vNestedPrefix_()
{
  unsigned int sf = getSizeFromVector(freqvector);
  if (pWA->getSize() !=  sf)
    throw Exception("WordFromIndependentFrequenciesSet: Size of the frequencies does not match size of the alphabet : " + TextTools::toString(sf) + " vs " + TextTools::toString(pWA->getSize()));

  unsigned int l = freqvector.size();

  for (unsigned int i = 0; i < l; i++)
  {
    vFreq_.push_back(freqvector[i]);
    vNestedPrefix_.push_back(freqvector[i]->getNamespace());
    vFreq_[i]->setNamespace("Word." + TextTools::toString(i + 1) + "_" + vNestedPrefix_[i]);
    addParameters_(vFreq_[i]->getParameters());
  }

  updateFrequencies();
}

WordFromIndependentFrequenciesSet::WordFromIndependentFrequenciesSet(const WordFromIndependentFrequenciesSet& iwfs) :
  WordFrequenciesSet(iwfs),
  vFreq_(iwfs.vFreq_.size()),
  vNestedPrefix_(iwfs.vNestedPrefix_)
{
  for (unsigned i = 0; i < iwfs.vFreq_.size(); i++)
    vFreq_[i] =  iwfs.vFreq_[i]->clone();
  updateFrequencies();
}

WordFromIndependentFrequenciesSet::~WordFromIndependentFrequenciesSet()
{
  for (unsigned i = 0; i < vFreq_.size(); i++)
    delete vFreq_[i];
}

WordFromIndependentFrequenciesSet& WordFromIndependentFrequenciesSet::operator=(const WordFromIndependentFrequenciesSet& iwfs)
{
  WordFrequenciesSet::operator=(iwfs);
  vNestedPrefix_ = iwfs.vNestedPrefix_;

  //Clean current frequencies first:
  for (unsigned i = 0; i < vFreq_.size(); i++)
    delete vFreq_[i];

  vFreq_.resize(iwfs.vFreq_.size());
  for (unsigned i = 0; i < vFreq_.size(); i++)
    vFreq_[i] = iwfs.vFreq_[i]->clone();
  updateFrequencies();

  return *this;
}

void WordFromIndependentFrequenciesSet::fireParameterChanged(const ParameterList& pl)
{
  unsigned int l = vFreq_.size();

  bool f = 0;
  for (unsigned int i = 0; i < l; i++)
  {
    f |= vFreq_[i]->matchParametersValues(pl);
  }

  if (f)
    updateFrequencies();
}

void WordFromIndependentFrequenciesSet::updateFrequencies()
{
  unsigned int l = vFreq_.size();
  unsigned int s = getAlphabet()->getSize();
  vector<double> f[l];

  unsigned int i, p, t, i2;

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

void WordFromIndependentFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if (frequencies.size() != getAlphabet()->getSize())
    throw DimensionException("WordFromIndependentFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  unsigned int size = frequencies.size();
  for (unsigned int i = 0; i < size; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.000001)
    throw Exception("WordFromIndependentFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");

  unsigned int d,i,j,s,l = vFreq_.size();
  int k;
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
    for (k = 0; k < (int)size; k++)
    {
      freq[(k / d) % s] += frequencies[k];
    }
    vFreq_[i]->setFrequencies(freq);
  }

  // updating freq_

  updateFrequencies();
}


unsigned int WordFromIndependentFrequenciesSet::getLength() const
{
  return vFreq_.size();
}

void WordFromIndependentFrequenciesSet::setNamespace(const std::string& prefix)
{
   AbstractFrequenciesSet::setNamespace(prefix);
  for (unsigned int i = 0; i < vFreq_.size(); i++)
  {
    vFreq_[i]->setNamespace(prefix + TextTools::toString(i + 1) + "_" + vNestedPrefix_[i]);
  }
}

std::string WordFromIndependentFrequenciesSet::getName() const
{
  string s = "Word From Independent Frequencies : ";
  for (unsigned int i = 0; i < vFreq_.size(); i++)
  {
    s += vFreq_[i]->getName();
  }
  return s;
}


// ///////////////////////////////////////////////////////////////////
// // WordFromUniqueFrequenciesSet


WordFromUniqueFrequenciesSet::WordFromUniqueFrequenciesSet(const WordAlphabet* pWA, FrequenciesSet* pabsfreq) :
  WordFrequenciesSet(pWA->getSize(), pWA),
  pFreq_(pabsfreq),
  NestedPrefix_(pabsfreq->getNamespace()),
  length_(pWA->getLength())
{
  unsigned int i;

  string st = "";
  for (i = 0; i < length_; i++)
  {
    st += TextTools::toString(i + 1);
  }

  pFreq_->setNamespace("Word." + st + "_" + NestedPrefix_);
  addParameters_(pFreq_->getParameters());

  updateFrequencies();
}

WordFromUniqueFrequenciesSet::WordFromUniqueFrequenciesSet(const WordFromUniqueFrequenciesSet& iwfs) :
  WordFrequenciesSet(iwfs),
  pFreq_(iwfs.pFreq_->clone()),
  NestedPrefix_(iwfs.NestedPrefix_),
  length_(iwfs.length_)
{
  updateFrequencies();
}


WordFromUniqueFrequenciesSet& WordFromUniqueFrequenciesSet::operator=(const WordFromUniqueFrequenciesSet& iwfs)
{
  WordFrequenciesSet::operator=(iwfs);
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
  if (pFreq_->matchParametersValues(pl))
    updateFrequencies();
}

void WordFromUniqueFrequenciesSet::updateFrequencies()
{
  unsigned int s = getAlphabet()->getSize();
  vector<double> f;
  int letsi = pFreq_->getAlphabet()->getSize();

  unsigned int i, p, i2;

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

void WordFromUniqueFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if (frequencies.size() != getAlphabet()->getSize())
    throw DimensionException("WordFromUniqueFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  unsigned int size = frequencies.size();
  for (unsigned int i = 0; i < size; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.000001)
    throw Exception("WordFromUniqueFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");

  unsigned int d,i,j;
  int k;
  vector<double> freq;

  unsigned int letsi = pFreq_->getAlphabet()->getSize();
  freq.resize(letsi);

  for (j = 0; j < letsi; j++)
  {
    freq[j] = 0;
  }

  d = size;
  for (i = 0; i < length_; i++)
  {
    d /= letsi;
    for (k = 0; k < (int)size; k++)
    {
      freq[(k / d) % letsi] += frequencies[k];
    }
  }
  for (j = 0; j < letsi; j++)
  {
    freq[j] /= length_;
  }
  pFreq_->setFrequencies(freq);

  // updating freq_

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


string WordFromUniqueFrequenciesSet::getName() const
{
  return "Word From Unique Frequency : " + pFreq_->getName() + " * " + TextTools::toString(length_);
}


/*********************************************************************/

FrequenciesSet* FrequenciesSet::getFrequenciesSetForCodons(short option, const GeneticCode& gc)
{
  FrequenciesSet* codonFreqs;
  const CodonAlphabet* pCA = dynamic_cast<const CodonAlphabet*>(gc.getSourceAlphabet());

  if (option == F0)
    codonFreqs = new FixedFrequenciesSet(pCA);
  else if (option == F1X4)
    codonFreqs = new WordFromUniqueFrequenciesSet(pCA, new FullNucleotideFrequenciesSet(pCA->getNucleicAlphabet()));
  else if (option == F3X4)
  {
    vector<FrequenciesSet*> v_AFS(3);
    v_AFS[0] = new FullNucleotideFrequenciesSet(pCA->getNucleicAlphabet());
    v_AFS[1] = new FullNucleotideFrequenciesSet(pCA->getNucleicAlphabet());
    v_AFS[2] = new FullNucleotideFrequenciesSet(pCA->getNucleicAlphabet());
    codonFreqs = new WordFromIndependentFrequenciesSet(pCA,v_AFS);
  }
  else if (option == F61)
    codonFreqs = new FullCodonFrequenciesSet(pCA);
  else throw Exception("FrequenciesSet::getFrequencySetForCodons(). Unvalid codon frequency set argument.");
  return codonFreqs;
}

/******************************************************************************/

const short FrequenciesSet::F0   = 0;
const short FrequenciesSet::F1X4 = 1;
const short FrequenciesSet::F3X4 = 2;
const short FrequenciesSet::F61  = 3;

/******************************************************************************/
