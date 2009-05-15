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

using namespace bpp;

#include <cmath>

using namespace std;

FullFrequenciesSet::FullFrequenciesSet(const Alphabet * alphabet, const string& name, const string& prefix):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, prefix)
{
  for(unsigned int i = 0; i < alphabet->getSize(); i++)
  {
    Parameter p(prefix + name + alphabet->intToChar((int)i), 1. / alphabet->getSize(), &Parameter::PROP_CONSTRAINT_IN);
    addParameter_(p);
    getFreq_(i) = 1. / alphabet->getSize();
  }
}

FullFrequenciesSet::FullFrequenciesSet(const Alphabet * alphabet, const vector<double>& initFreqs, const string& name, const string& prefix) throw (Exception):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, prefix)
{
  if(initFreqs.size() != alphabet->getSize())
    throw Exception("FullFrequenciesSet(constructor). There must be " + TextTools::toString(alphabet->getSize()) + " frequencies.");
  double sum = 0.0;
  for(unsigned int i = 0; i < initFreqs.size(); i++)
  {
    sum += initFreqs[i];
  }
  if(fabs(1-sum) > 0.000001)
  {
    throw Exception("Root frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  }
  for(unsigned int i = 0; i < alphabet->getSize(); i++)
  {
    Parameter p(prefix + name + alphabet->intToChar((int)i), initFreqs[i], &Parameter::PROP_CONSTRAINT_IN);
    addParameter_(p);
    getFreq_(i) = initFreqs[i];
  }
}

FullNAFrequenciesSet::FullNAFrequenciesSet(const NucleicAlphabet * alphabet, const string& prefix):
  AbstractFrequenciesSet(4, alphabet, prefix)
{
  Parameter thetaP(prefix + "theta" , 0.5, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(thetaP);
  Parameter theta1P(prefix + "theta1", 0.5, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(theta1P);
  Parameter theta2P(prefix + "theta2", 0.5, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(theta2P);
  getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
}

FullNAFrequenciesSet::FullNAFrequenciesSet(const NucleicAlphabet * alphabet, double theta, double theta1, double theta2, const string& prefix):
  AbstractFrequenciesSet(4, alphabet, prefix)
{
  Parameter thetaP(prefix + "theta" , theta, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(thetaP);
  Parameter theta1P(prefix + "theta1", theta1, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(theta1P);
  Parameter theta2P(prefix + "theta2", theta2, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(theta2P);
  getFreq_(0) = theta1 * (1. - theta);
  getFreq_(1) = (1 - theta2) * theta;
  getFreq_(2) = theta2 * theta;
  getFreq_(3) = (1 - theta1) * (1. - theta);
}

void FullNAFrequenciesSet::fireParameterChanged(const ParameterList & pl)
{
  double theta  = getParameter_(0).getValue();
  double theta1 = getParameter_(1).getValue();
  double theta2 = getParameter_(2).getValue();
  getFreq_(0) = theta1 * (1. - theta);
  getFreq_(1) = (1 - theta2) * theta;
  getFreq_(2) = theta2 * theta;
  getFreq_(3) = (1 - theta1) * (1. - theta);}

FullProteinFrequenciesSet::FullProteinFrequenciesSet(const ProteicAlphabet * alphabet, const vector<double> & initFreqs, const string & prefix) throw (Exception):
  AbstractFrequenciesSet(20, alphabet, prefix)
{
  if(initFreqs.size() != 20)
    throw Exception("FullProteinFrequenciesSet(constructor). There must be 20 frequencies.");
  for(unsigned int i = 1; i < 20; i++)
  {
    Parameter p(prefix + "theta" + TextTools::toString(i) , 0.05 / (1-0.05*((double)i - 1.)), &Parameter::PROP_CONSTRAINT_EX);
    addParameter_(p);
  }
  setFrequencies(initFreqs);
}

FullProteinFrequenciesSet::FullProteinFrequenciesSet(const ProteicAlphabet * alphabet, const string& prefix):
  AbstractFrequenciesSet(20, alphabet, prefix)
{
  for(unsigned int i = 1; i < 20; i++)
  {
    Parameter p(prefix + "theta" + TextTools::toString(i) , 0.05 / (1-0.05*((double)i - 1.)), &Parameter::PROP_CONSTRAINT_EX);
    addParameter_(p);
    getFreq_ (i - 1) = 0.05;
  }
  getFreq_(19) = 0.05;
}

void FullProteinFrequenciesSet::setFrequencies(const vector<double> & frequencies)
  throw (DimensionException, Exception)
{
  if(frequencies.size() != 20) throw DimensionException("FullProteinFrequenciesSet::setFrequencies", frequencies.size(), 20);
  double sum = 0.0;
  for(unsigned int i = 0; i < 20; i++)
  {
    sum += frequencies[i];
  }
  if(fabs(1.-sum) > 0.000001)
  {
    throw Exception("FullProteinFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  }
  setFrequencies_(frequencies);
  double cumFreq = 1.;
  for(unsigned int i = 0; i < 19; i++)
  {
    double theta = getFreq_(i) / cumFreq;
    cumFreq -= getFreq_(i);
    getParameter_(i).setValue(theta);
  }
}

void FullProteinFrequenciesSet::fireParameterChanged(const ParameterList & pl)
{
  double cumTheta = 1.;
  for(unsigned int i = 0; i < 19; i++)
  {
    double theta = getParameter_(i).getValue();
    getFreq_(i) = cumTheta * theta;
    cumTheta *= (1. - theta);
  }
  getFreq_(19) = cumTheta;
}

FixedFrequenciesSet::FixedFrequenciesSet(const Alphabet * alphabet, const vector<double>& initFreqs, const string& prefix):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, prefix)
{
  setFrequencies(initFreqs);
}

