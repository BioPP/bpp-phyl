//
// File: NucleotideFrequenciesSet.cpp
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

#include "NucleotideFrequenciesSet.h"

#include <Bpp/Numeric/NumConstants.h>

using namespace bpp;

#include <cmath>
using namespace std;

// ///////////////////////////////////////
// FullNucleotideFrequenciesSet


FullNucleotideFrequenciesSet::FullNucleotideFrequenciesSet(
  const NucleicAlphabet* alphabet, bool allowNullFreqs,
  const string& name) :
  AbstractFrequenciesSet(new CanonicalStateMap(alphabet, false), "Full.", name)
{
  addParameter_(new Parameter(
    "Full.theta", 0.5,
    allowNullFreqs ?
    &Parameter::PROP_CONSTRAINT_IN :
    &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter(
    "Full.theta1", 0.5,
    allowNullFreqs ?
    &Parameter::PROP_CONSTRAINT_IN :
    &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter("Full.theta2", 0.5,
                    allowNullFreqs ?
                    &Parameter::PROP_CONSTRAINT_IN :
                    &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
}

FullNucleotideFrequenciesSet::FullNucleotideFrequenciesSet(
  const NucleicAlphabet* alphabet, double theta, double theta1, double theta2,
  bool allowNullFreqs, const string& name) :
  AbstractFrequenciesSet(new CanonicalStateMap(alphabet, false), "Full.", name)
{
  addParameter_(new Parameter(
    "Full.theta",
    theta,
    allowNullFreqs ?
    &Parameter::PROP_CONSTRAINT_IN :
    &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter(
    "Full.theta1",
    theta1,
    allowNullFreqs ?
    &Parameter::PROP_CONSTRAINT_IN :
    &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter(
    "Full.theta2",
    theta2,
    allowNullFreqs ?
    &Parameter::PROP_CONSTRAINT_IN :
    &Parameter::PROP_CONSTRAINT_EX));
  getFreq_(0) = theta1 * (1. - theta);
  getFreq_(1) = (1 - theta2) * theta;
  getFreq_(2) = theta2 * theta;
  getFreq_(3) = (1 - theta1) * (1. - theta);
}

void FullNucleotideFrequenciesSet::setFrequencies(const vector<double>& frequencies) 
{
  if (frequencies.size() != 4) throw DimensionException(" FullNucleotideFrequenciesSet::setFrequencies", frequencies.size(), 4);
  double sum = 0.0;
  for (unsigned int i = 0; i < 4; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > NumConstants::SMALL())
    throw Exception("FullNucleotideFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  double theta = frequencies[1] + frequencies[2];
  getParameter_(0).setValue(theta);
  getParameter_(1).setValue(frequencies[0] / (1 - theta));
  getParameter_(2).setValue(frequencies[2] / theta);

  setFrequencies_(frequencies);
}

void FullNucleotideFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
  AbstractFrequenciesSet::fireParameterChanged(parameters);
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

void GCFrequenciesSet::setFrequencies(const vector<double>& frequencies) 
{
  if (frequencies.size() != 4) throw DimensionException("GCFrequenciesSet::setFrequencies", frequencies.size(), 4);
  double sum = 0.0;
  for (unsigned int i = 0; i < 4; i++)
  {
    sum += frequencies[i];
  }
  if (fabs(1. - sum) > NumConstants::SMALL())
    throw Exception("GCFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  double theta = frequencies[1] + frequencies[2];
  // We set everything in one shot here:
  getParameter_(0).setValue(theta);
  getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
  getFreq_(1) = getFreq_(2) = theta / 2.;
}

void GCFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
  AbstractFrequenciesSet::fireParameterChanged(parameters);
  double theta = getParameter_(0).getValue();
  getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
  getFreq_(1) = getFreq_(2) = theta / 2.;
}


