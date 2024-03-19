// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/Simplex.h>

#include "FrequencySet.h"

// From bpp-phyl
#include "../SubstitutionModel.h"
#include "../AbstractBiblioSubstitutionModel.h"


using namespace bpp;

#include <cmath>
using namespace std;

shared_ptr<IntervalConstraint> FrequencySetInterface::FREQUENCE_CONSTRAINT_MILLI(new IntervalConstraint(NumConstants::MILLI(), 1 - NumConstants::MILLI(), false, false));
shared_ptr<IntervalConstraint> FrequencySetInterface::FREQUENCE_CONSTRAINT_CENTI(new IntervalConstraint(NumConstants::CENTI(), 1 - NumConstants::CENTI(), false, false));
shared_ptr<IntervalConstraint> FrequencySetInterface::FREQUENCE_CONSTRAINT_SMALL(new IntervalConstraint(NumConstants::SMALL(), 1 - NumConstants::SMALL(), false, false));

// ///////////////////////////////////////
// AbstractFrequencySet


void AbstractFrequencySet::setFrequenciesFromAlphabetStatesFrequencies(const map<int, double>& frequencies)
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

  if (x != 0)
    for (size_t i = 0; i < s; ++i)
    {
      freq[i] /= x;
    }
  else
    for (size_t i = 0; i < s; ++i)
    {
      freq[i] = 1.0 / (double)s;
    }

  setFrequencies(freq);
}

const std::map<int, double> AbstractFrequencySet::getAlphabetStatesFrequencies() const
{
  map<int, double> fmap;
  for (size_t i = 0; i < stateMap_->getNumberOfModelStates(); ++i)
  {
    fmap[stateMap_->getAlphabetStateAsInt(i)] += freq_[i];
  }
  return fmap;
}

// ////////////////////////////
// FullFrequencySet

FullFrequencySet::FullFrequencySet(
    shared_ptr<const StateMapInterface> stateMap,
    bool allowNullFreqs,
    unsigned short method,
    const string& name) :
  AbstractFrequencySet(
    stateMap,
    "Full.",
    name),
  sFreq_(stateMap->getNumberOfModelStates(), method, allowNullFreqs, "Full.")
{
  vector<double> vd;
  double r = 1. / static_cast<double>(stateMap->getNumberOfModelStates());

  for (size_t i = 0; i < stateMap->getNumberOfModelStates(); i++)
  {
    vd.push_back(r);
  }

  sFreq_.setFrequencies(vd);
  addParameters_(sFreq_.getParameters());
  updateFreq_();
}

FullFrequencySet::FullFrequencySet(
    shared_ptr<const StateMapInterface> stateMap,
    const vector<double>& initFreqs,
    bool allowNullFreqs,
    unsigned short method,
    const string& name) :
  AbstractFrequencySet(
    stateMap,
    "Full.",
    name),
  sFreq_(stateMap->getNumberOfModelStates(), method, allowNullFreqs, "Full.")
{
  sFreq_.setFrequencies(initFreqs);
  addParameters_(sFreq_.getParameters());
  updateFreq_();
}

void FullFrequencySet::setFrequencies(const vector<double>& frequencies)
{
  sFreq_.setFrequencies(frequencies);
  setParametersValues(sFreq_.getParameters());

  updateFreq_();
}

void FullFrequencySet::setNamespace(const std::string& nameSpace)
{
  sFreq_.setNamespace(nameSpace);
  AbstractFrequencySet::setNamespace(nameSpace);
}

void FullFrequencySet::fireParameterChanged(const ParameterList& parameters)
{
  sFreq_.matchParametersValues(parameters);
  updateFreq_();
}

void FullFrequencySet::updateFreq_()
{
  for (size_t i = 0; i < getAlphabet()->getSize(); i++)
  {
    getFreq_(i) = sFreq_.prob(i);
  }
}

// ///////////////////////////////////////////
// / FixedFrequencySet

FixedFrequencySet::FixedFrequencySet(
    shared_ptr<const StateMapInterface> stateMap,
    const vector<double>& initFreqs,
    const string& name) :
  AbstractFrequencySet(
    stateMap,
    "Fixed.",
    name)
{
  if (stateMap->getNumberOfModelStates() != initFreqs.size())
    throw Exception("FixedFrequencySet::constructor. size of init vector does not match the number of states in the model.");
  setFrequencies(initFreqs);
}

FixedFrequencySet::FixedFrequencySet(
    shared_ptr<const StateMapInterface> stateMap,
    const string& name) :
  AbstractFrequencySet(
    stateMap,
    "Fixed.",
    name)
{
  size_t n = stateMap->getNumberOfModelStates();
  setFrequencies_(std::vector<double>(n, 1. / (double)n));
}

void FixedFrequencySet::setFrequencies(const vector<double>& frequencies)
{
  if (frequencies.size() != getNumberOfFrequencies())
    throw DimensionException("FixedFrequencySet::setFrequencies", frequencies.size(), getNumberOfFrequencies());
  double sum = 0.0;
  for (size_t i = 0; i < frequencies.size(); i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.00001)
    throw Exception("FixedFrequencySet::setFrequencies. Frequencies sum must equal 1 (sum = " + TextTools::toString(sum) + ").");
  setFrequencies_(frequencies);
}

MarkovModulatedFrequencySet::MarkovModulatedFrequencySet(
    unique_ptr<FrequencySetInterface> freqSet,
    const std::vector<double>& rateFreqs) :
  AbstractFrequencySet(
    make_shared<MarkovModulatedStateMap>(freqSet->stateMap(), static_cast<unsigned int>(rateFreqs.size())),
    "MarkovModulated.",
    "MarkovModulated." + freqSet->getName()),
  freqSet_(std::move(freqSet)),
  rateFreqs_(rateFreqs)
{
  freqSet_->setNamespace(std::string("MarkovModulated.") + freqSet_->getNamespace());
  addParameters_(freqSet_->getParameters());
  setFrequencies_(VectorTools::kroneckerMult(rateFreqs, freqSet_->getFrequencies()));
}

//////////////////////////////////
/// From Model


FromModelFrequencySet::FromModelFrequencySet(const FromModelFrequencySet& fmfs) :
  AbstractFrequencySet(fmfs),
  model_(fmfs.model_->clone())
{}

FromModelFrequencySet& FromModelFrequencySet::operator=(const FromModelFrequencySet& fmfs)
{
  AbstractFrequencySet::operator=(fmfs);
  model_.reset(fmfs.model_->clone());
  return *this;
}

FromModelFrequencySet::~FromModelFrequencySet() {}

FromModelFrequencySet::FromModelFrequencySet(
    shared_ptr<TransitionModelInterface> model) :
  AbstractFrequencySet(model->getStateMap(), "FromModel." + (model ? model->getNamespace() : ""), "FromModel"),
  model_(model)
{
  model_->setNamespace(getNamespace());
  addParameters_(model_->getParameters());
  setFrequencies_(model_->getFrequencies());
}


void FromModelFrequencySet::setNamespace(const std::string& name)
{
  AbstractParameterAliasable::setNamespace(name);
  model_->setNamespace(name);
}


void FromModelFrequencySet::setFrequencies(const std::vector<double>& frequencies)
{
  std::map<int, double> freq;
  for (size_t i = 0; i < getNumberOfFrequencies(); ++i)
  {
    freq[stateMap().getAlphabetStateAsInt(i)] += frequencies[i];
  }
  model_->setFreq(freq);
  matchParametersValues(model_->getParameters());
}

void FromModelFrequencySet::fireParameterChanged(const ParameterList& pl)
{
  model_->matchParametersValues(pl);
  setFrequencies_(model_->getFrequencies());
}


//////////////////////////////////
/// User

UserFrequencySet::UserFrequencySet(
    shared_ptr<const StateMapInterface> stateMap,
    const std::string& path,
    size_t nCol) :
  AbstractFrequencySet(
    stateMap,
    "Empirical.",
    "Empirical"),
  path_(path),
  nCol_(nCol)
{
  readFromFile_();
}

UserFrequencySet::UserFrequencySet(const UserFrequencySet& fmfs) :
  AbstractFrequencySet(fmfs),
  path_(fmfs.path_),
  nCol_(fmfs.nCol_)
{}

UserFrequencySet& UserFrequencySet::operator=(const UserFrequencySet& fmfs)
{
  AbstractFrequencySet::operator=(fmfs);
  path_ = fmfs.path_;
  nCol_ = fmfs.nCol_;
  return *this;
}

void UserFrequencySet::readFromFile_()
{
  if (!FileTools::fileExists(path_.c_str()))
    throw Exception("UserFrequencySet::readFromFile. Frequencies file not found : " +  path_);

  ifstream in(path_.c_str(), ios::in);

  // Read profile:

  for (unsigned int i = 0; i < getAlphabet()->getSize(); i++)
  {
    if (!in)
      throw Exception("UserFrequencySet::readFromFile. Missing frequencies in file : " +  path_);

    string line = FileTools::getNextLine(in);
    StringTokenizer st(line);
    double s(0);
    for (unsigned int j = 0; j < nCol_; j++)
    {
      if (!st.hasMoreToken())
        throw Exception("UserFrequencySet::readFromFile. Missing frequencies for column " + TextTools::toString(nCol_) + " in line " + TextTools::toString(i));
      s = TextTools::toDouble(st.nextToken());
    }
    getFreq_(i) = s;
  }

  double sf = VectorTools::sum(getFrequencies_());
  if (fabs(sf - 1) > 0.000001)
  {
    ApplicationTools::displayMessage("WARNING!!! Frequencies sum to " + TextTools::toString(sf) + ", frequencies have been scaled.");
    sf *= 1. / sf;
  }

  // Closing stream:
  in.close();
}

void UserFrequencySet::setFrequencies(const vector<double>& frequencies)
{
  if (frequencies.size() != getNumberOfFrequencies())
    throw DimensionException("UserFrequencySet::setFrequencies", frequencies.size(), getNumberOfFrequencies());
  double sum = 0.0;
  for (size_t i = 0; i < frequencies.size(); i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > 0.00001)
    throw Exception("FixedFrequencySet::setFrequencies. Frequencies sum must equal 1 (sum = " + TextTools::toString(sum) + ").");
  setFrequencies_(frequencies);
}
