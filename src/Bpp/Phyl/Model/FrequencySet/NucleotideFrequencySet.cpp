// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/NumConstants.h>

#include "NucleotideFrequencySet.h"

using namespace bpp;

#include <cmath>
using namespace std;

// ///////////////////////////////////////
// FullNucleotideFrequencySet


FullNucleotideFrequencySet::FullNucleotideFrequencySet(
    shared_ptr<const NucleicAlphabet> alphabet,
    bool allowNullFreqs,
    const string& name) :
  AbstractFrequencySet(
    make_shared<CanonicalStateMap>(alphabet, false),
    "Full.",
    name)
{
  addParameter_(new Parameter(
                  "Full.theta", 0.5,
                  allowNullFreqs ?
                  Parameter::PROP_CONSTRAINT_IN :
                  FrequencySetInterface::FREQUENCE_CONSTRAINT_CENTI));
  addParameter_(new Parameter(
                  "Full.theta1", 0.5,
                  allowNullFreqs ?
                  Parameter::PROP_CONSTRAINT_IN :
                  FrequencySetInterface::FREQUENCE_CONSTRAINT_CENTI));
  addParameter_(new Parameter("Full.theta2", 0.5,
                              allowNullFreqs ?
                              Parameter::PROP_CONSTRAINT_IN :
                              FrequencySetInterface::FREQUENCE_CONSTRAINT_CENTI));
  getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
}

FullNucleotideFrequencySet::FullNucleotideFrequencySet(
    shared_ptr<const NucleicAlphabet> alphabet,
    double theta, 
    double theta1,
    double theta2,
    bool allowNullFreqs, 
    const string& name) :
  AbstractFrequencySet(
    make_shared<CanonicalStateMap>(alphabet, false),
    "Full.",
    name)
{
  addParameter_(new Parameter(
                  "Full.theta",
                  theta,
                  allowNullFreqs ?
                  Parameter::PROP_CONSTRAINT_IN :
                  FrequencySetInterface::FREQUENCE_CONSTRAINT_CENTI));
  addParameter_(new Parameter(
                  "Full.theta1",
                  theta1,
                  allowNullFreqs ?
                  Parameter::PROP_CONSTRAINT_IN :
                  FrequencySetInterface::FREQUENCE_CONSTRAINT_CENTI));
  addParameter_(new Parameter(
                  "Full.theta2",
                  theta2,
                  allowNullFreqs ?
                  Parameter::PROP_CONSTRAINT_IN :
                  Parameter::PROP_CONSTRAINT_EX));
  getFreq_(0) = theta1 * (1. - theta);
  getFreq_(1) = (1 - theta2) * theta;
  getFreq_(2) = theta2 * theta;
  getFreq_(3) = (1 - theta1) * (1. - theta);
}

void FullNucleotideFrequencySet::setFrequencies(const vector<double>& frequencies)
{
  if (frequencies.size() != 4)
    throw DimensionException(" FullNucleotideFrequencySet::setFrequencies", frequencies.size(), 4);
  double sum = 0.0;
  for (auto& f : frequencies) sum += f;
  if (fabs(1. - sum) > NumConstants::SMALL())
    throw Exception("FullNucleotideFrequencySet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  double theta = frequencies[1] + frequencies[2];
  getParameter_(0).setValue(theta);
  getParameter_(1).setValue(frequencies[0] / (1 - theta));
  getParameter_(2).setValue(frequencies[2] / theta);

  setFrequencies_(frequencies);
}

void FullNucleotideFrequencySet::fireParameterChanged(const ParameterList& parameters)
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
// GCFrequencySet

void GCFrequencySet::setFrequencies(const vector<double>& frequencies)
{
  if (frequencies.size() != 4)
    throw DimensionException("GCFrequencySet::setFrequencies", frequencies.size(), 4);
  double sum = 0.0;
  for (unsigned int i = 0; i < 4; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > NumConstants::SMALL())
    throw Exception("GCFrequencySet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  double theta = frequencies[1] + frequencies[2];
  // We set everything in one shot here:
  getParameter_(0).setValue(theta);
  getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
  getFreq_(1) = getFreq_(2) = theta / 2.;
}

void GCFrequencySet::fireParameterChanged(const ParameterList& parameters)
{
  double theta = getParameter_(0).getValue();
  getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
  getFreq_(1) = getFreq_(2) = theta / 2.;
}
