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
        "Full.theta_" + TextTools::toString(i + 1),
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
        "Full.theta_" + TextTools::toString(i + 1),
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
    getParameter_("theta_" + TextTools::toString(i + 1)).setValue(frequencies[i] / y);
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
    getFreq_(i) = getParameter_("theta_" + TextTools::toString(i + 1)).getValue() * y;
    y *= 1 - getParameter_("theta_" + TextTools::toString(i + 1)).getValue();
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
          "Full.theta_" + TextTools::toString(i + 1),
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
          "Full.theta_" + TextTools::toString(i + 1),
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
      getParameter_("theta_" + TextTools::toString(i + 1)).setValue(frequencies[i] / y);
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
      getFreq_(i) = getParameter_("theta_" + TextTools::toString(i + 1)).getValue() * y;
      y *= 1 - getParameter_("theta_" + TextTools::toString(i + 1)).getValue();
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
      "Full.theta_1", 0.5,
      allowNullFreqs ?
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
  addParameter_(theta1P);
  Parameter theta2P("Full.theta_2", 0.5,
      allowNullFreqs ?
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
  addParameter_(theta2P);
  getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
}

FullNucleotideFrequenciesSet::FullNucleotideFrequenciesSet(const NucleicAlphabet* alphabet, double theta, double theta_1, double theta_2, bool allowNullFreqs) :
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
      "Full.theta_1",
      theta_1,
      allowNullFreqs ?
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
  addParameter_(theta1P);
  Parameter theta2P(
      "Full.theta_2",
      theta_2,
      allowNullFreqs ?
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_IN) :
      dynamic_cast<Constraint*>(&Parameter::PROP_CONSTRAINT_EX));
  addParameter_(theta2P);
  getFreq_(0) = theta_1 * (1. - theta);
  getFreq_(1) = (1 - theta_2) * theta;
  getFreq_(2) = theta_2 * theta;
  getFreq_(3) = (1 - theta_1) * (1. - theta);
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
  double theta_1 = getParameter_(1).getValue();
  double theta_2 = getParameter_(2).getValue();
  getFreq_(0) = theta_1 * (1. - theta);
  getFreq_(1) = (1 - theta_2) * theta;
  getFreq_(2) = theta_2 * theta;
  getFreq_(3) = (1 - theta_1) * (1. - theta);
}

///////////////////////////////////////////
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


/////////////////////////////////////////////
/// FixedFrequenciesSet

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

/////////////////////////////////////////////
/// FixedCodonFrequenciesSet

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
// IndependentWordFrequenciesSet

IndependentWordFrequenciesSet::IndependentWordFrequenciesSet(
    const std::vector<FrequenciesSet*>& freqvector) :
  AbstractFrequenciesSet(
      getSizeFromVector(freqvector),
      extractAlph_(freqvector), "IndependentWord."),
  vFreq_(),
  vNestedPrefix_(),
  uniqueAbsFreq_(false)
{
  int i, j, k, t;
  int l = freqvector.size();

  i = 0;
  j = 1;
  while (!uniqueAbsFreq_ && i < (l - 1))
  {
    if (freqvector[i] == freqvector[j])
      uniqueAbsFreq_ = true;
    else
    {
      j++;
      if (j == l)
      {
        i++;
        j = i + 1;
      }
    }
  }

  if (!uniqueAbsFreq_)
  {
    for (i = 0; i < l; i++)
    {
      vFreq_.push_back(freqvector[i]);
      vNestedPrefix_.push_back(freqvector[i]->getNamespace());
      vFreq_[i]->setNamespace("IndependentWord." + TextTools::toString(i + 1) + "_" + vNestedPrefix_[i]);
      addParameters_(vFreq_[i]->getParameters());
    }
  }
  else
  {
    string st = "";
    for (i = 0; i < l; i++)
    {
      vFreq_.push_back(freqvector[0]);
      vNestedPrefix_.push_back(freqvector[0]->getNamespace());
      st += TextTools::toString(i + 1);
    }
    vFreq_[0]->setNamespace("IndependentWord." + st + "_" + vNestedPrefix_[0]);
    addParameters_(vFreq_[0]->getParameters());
  }

  vector<double> f[l];
  for (i = 0; i < l; i++)
  {
    f[i] = vFreq_[i]->getFrequencies();
  }

  int s = getAlphabet()->getSize();

  for (i = 0; i < s; i++)
  {
    j = i;
    getFreq_(j) = 1;
    for (k = l - 1; k >= 0; k--)
    {
      t = vFreq_[k]->getAlphabet()->getSize();
      getFreq_(i) *= f[k][j % t];
      j /= t;
    }
  }
}

IndependentWordFrequenciesSet::IndependentWordFrequenciesSet(
    FrequenciesSet* pabsfreq, int num) :
  AbstractFrequenciesSet(
      static_cast<int>(pow(static_cast<double>(pabsfreq->getAlphabet()->getSize()), num)),
      new WordAlphabet(pabsfreq->getAlphabet(), num), "IndependentWord."),
  vFreq_(),
  vNestedPrefix_(),
  uniqueAbsFreq_(true)
{
  int i, j, k, t;

  string st = "";
  for (i = 0; i < num; i++)
  {
    vFreq_.push_back(pabsfreq);
    vNestedPrefix_.push_back(pabsfreq->getNamespace());
    st += TextTools::toString(i + 1);
  }
  vFreq_[0]->setNamespace("IndependentWord." + st + "_" + vNestedPrefix_[0]);
  addParameters_(vFreq_[0]->getParameters());

  vector<double> f[num];
  for (i = 0; i < num; i++)
  {
    f[i] = vFreq_[i]->getFrequencies();
  }

  int s = getAlphabet()->getSize();

  for (i = 0; i < s; i++)
  {
    j = i;
    getFreq_(j) = 1;
    for (k = num - 1; k >= 0; k--)
    {
      t = vFreq_[k]->getAlphabet()->getSize();
      getFreq_(i) *= f[k][j % t];
      j /= t;
    }
  }
}

Alphabet* IndependentWordFrequenciesSet::extractAlph_(const std::vector<FrequenciesSet*>& freqVector)
{
  unsigned int i, l = freqVector.size();
  vector<const Alphabet*> pval;

  for (i = 0; i < l; i++)
  {
    pval.push_back(freqVector[i]->getAlphabet());
  }

  return new WordAlphabet(pval);
}

unsigned int IndependentWordFrequenciesSet::getSizeFromVector(const std::vector<FrequenciesSet*>& freqVector)
{
  unsigned int s = 1;
  unsigned int l = freqVector.size();

  for (unsigned int i = 0; i < l; i++)
  {
    s *= freqVector[i]->getAlphabet()->getSize();
  }

  return s;
}

IndependentWordFrequenciesSet::IndependentWordFrequenciesSet(const IndependentWordFrequenciesSet& iwfs) :
  AbstractFrequenciesSet(
      iwfs.getNumberOfFrequencies(),
      new WordAlphabet(*dynamic_cast<const WordAlphabet*>(iwfs.getAlphabet())), iwfs.getNamespace()),
  vFreq_(iwfs.vFreq_),
  vNestedPrefix_(iwfs.vNestedPrefix_),
  uniqueAbsFreq_(iwfs.uniqueAbsFreq_)
{
  unsigned int l = iwfs.vFreq_.size();

  FrequenciesSet* pAFS = (uniqueAbsFreq_ ? iwfs.vFreq_[0]->clone() : 0);

  for (unsigned i = 0; i < l; i++)
  {
    vFreq_.push_back(uniqueAbsFreq_ ? pAFS : iwfs.vFreq_[i]->clone());
  }
}

IndependentWordFrequenciesSet& IndependentWordFrequenciesSet::operator=(const IndependentWordFrequenciesSet& iwfs)
{
  AbstractFrequenciesSet::operator=(iwfs);
  vFreq_ = iwfs.vFreq_;
  vNestedPrefix_ = iwfs.vNestedPrefix_;
  uniqueAbsFreq_ = iwfs.uniqueAbsFreq_;
  
  unsigned int l = iwfs.vFreq_.size();

  FrequenciesSet* pAFS = (uniqueAbsFreq_ ? iwfs.vFreq_[0]->clone() : 0);

  for (unsigned i = 0; i < l; i++)
  {
    vFreq_[i] = (uniqueAbsFreq_ ? pAFS : iwfs.vFreq_[i]->clone());
  }
  return *this;
}


void IndependentWordFrequenciesSet::fireParameterChanged(const ParameterList& pl)
{
  unsigned int l = vFreq_.size();

  bool f = 0;
  if (uniqueAbsFreq_)
    f = vFreq_[0]->matchParametersValues(pl);
  else
    for (unsigned int i = 0; i < l; i++)
    {
      f |= vFreq_[i]->matchParametersValues(pl);
    }

  if (f)
    updateFrequencies();
}

void IndependentWordFrequenciesSet::updateFrequencies()
{
  unsigned int l = vFreq_.size();
  unsigned int s = getAlphabet()->getSize();
  vector<double> f[l];

  unsigned int i, p, t, i2;
// [Julien, 30/10/09] : We should not need that, has frequencies are automatically updated when a parameter is changed...
//  if (unique_AbsFreq)
//    vFreq_[0]->updateFrequencies();
//  else
//    for (i = 0; i < l; i++)
//      vFreq_[i]->updateFrequencies();

  for (i = 0; i < l; i++)
  {
    f[i] = vFreq_[i]->getFrequencies();
  }

  for (i = 0; i < s; i++)
  {
    i2 = i;
    getFreq_(i2) = 1;
    for (p = l; p > 0; p--)
    {
      t = vFreq_[p - 1]->getAlphabet()->getSize();
      getFreq_(i) *= f[p - 1][i2 % t];
      i2 /= t;
    }
  }
}

void IndependentWordFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if (frequencies.size() != getAlphabet()->getSize()) throw DimensionException("IndependentWordFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  unsigned int size = frequencies.size();
  for (unsigned int i = 0; i < size; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.000001)
    throw Exception("IndependentWordFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");

  unsigned int d,i,j,s,l = vFreq_.size();
  int k;
  vector<double> freq;

  if (uniqueAbsFreq_)
  {
    s = vFreq_[0]->getAlphabet()->getSize();
    freq.resize(s);
    for (j = 0; j < s; j++)
    {
      freq[j] = 0;
    }
    d = size;
    for (i = 0; i < l; i++)
    {
      d /= s;
      for (k = 0; k < (int)size; k++)
      {
        freq[(k / d) % s] += frequencies[k];
      }
    }
    for (j = 0; j < s; j++)
    {
      freq[j] /= l;
    }
    vFreq_[0]->setFrequencies(freq);
  }
  else
  {
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
  }

  // updating freq_

  vector<double> f[l];
  for (k = 0; k < (int)l; k++)
  {
    f[k] = vFreq_[k]->getFrequencies();
  }

  for (i = 0; i < size; i++)
  {
    j = i;
    getFreq_(j) = 1;
    for (k = l - 1; k >= 0; k--)
    {
      s = vFreq_[k]->getAlphabet()->getSize();
      getFreq_(i) *= f[k][j % s];
      j /= s;
    }
  }
}

IndependentWordFrequenciesSet::~IndependentWordFrequenciesSet()
{
  if (getAlphabet())
    delete getAlphabet();
}

void IndependentWordFrequenciesSet::setNamespace(const string& prefix)
{
  AbstractFrequenciesSet::setNamespace(prefix);
  if (uniqueAbsFreq_)
  {
    string st = "";
    for (unsigned i = 0; i < vFreq_.size(); i++)
      st += TextTools::toString(i + 1);
    vFreq_[0]->setNamespace(prefix + st + "_" + vNestedPrefix_[0]);
  }
  else
  {
    for (unsigned int i = 0; i < vFreq_.size(); i++)
    {
      vFreq_[i]->setNamespace(prefix + TextTools::toString(i + 1) + "_" + vNestedPrefix_[i]);
    }
  }
}

string IndependentWordFrequenciesSet::getName() const
{
  string s = "IndependentWord : ";
  for (unsigned int i = 0; i < vFreq_.size(); i++)
  {
    s += vFreq_[i]->getName();
  }
  return s;
}


FrequenciesSet* FrequenciesSet::getFrequencySetForCodons(short option, const GeneticCode& gc)
{
  FrequenciesSet* codonFreqs;
  if (option == F0)
    codonFreqs = new FixedFrequenciesSet(dynamic_cast<const CodonAlphabet*>(gc.getSourceAlphabet()));
  else if (option == F1X4)
    codonFreqs = new IndependentWordFrequenciesSet(new FullNucleotideFrequenciesSet(dynamic_cast<const CodonAlphabet*>(gc.getSourceAlphabet())->getNucleicAlphabet()), 3);
  else if (option == F3X4)
  {
    vector<FrequenciesSet*> v_AFS(3);
    v_AFS[0] = new FullNucleotideFrequenciesSet(dynamic_cast<const CodonAlphabet*>(gc.getSourceAlphabet())->getNucleicAlphabet());
    v_AFS[1] = new FullNucleotideFrequenciesSet(dynamic_cast<const CodonAlphabet*>(gc.getSourceAlphabet())->getNucleicAlphabet());
    v_AFS[2] = new FullNucleotideFrequenciesSet(dynamic_cast<const CodonAlphabet*>(gc.getSourceAlphabet())->getNucleicAlphabet());
    codonFreqs = new IndependentWordFrequenciesSet(v_AFS);
  }
  else if (option == F61)
    codonFreqs = new FullCodonFrequenciesSet(dynamic_cast<const CodonAlphabet*>(gc.getSourceAlphabet()));
  else throw Exception("FrequenciesSet::getFrequencySetForCodons(). Unvalid codon frequency set argument.");
  return codonFreqs;
}

/******************************************************************************/

const short FrequenciesSet::F0   = 0;
const short FrequenciesSet::F1X4 = 1;
const short FrequenciesSet::F3X4 = 2;
const short FrequenciesSet::F61  = 3;

/******************************************************************************/

