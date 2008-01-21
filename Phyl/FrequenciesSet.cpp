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

FullFrequenciesSet::FullFrequenciesSet(const Alphabet * alphabet, const string & prefix):
  AbstractFrequenciesSet(alphabet)
{
  _freq.resize(alphabet->getSize());
  for(unsigned int i = 0; i < alphabet->getSize(); i++)
  {
    _parameters.addParameter(Parameter(prefix + alphabet->intToChar((int)i), 1. / alphabet->getSize(), &Parameter::PROP_CONSTRAINT_IN));
    _freq[i] = 1. / alphabet->getSize();
  }
}

FullFrequenciesSet::FullFrequenciesSet(const Alphabet * alphabet, const vector<double> & initFreqs, const string & prefix) throw (Exception):
  AbstractFrequenciesSet(alphabet)
{
  if(initFreqs.size() != alphabet->getSize())
    throw Exception("FullFrequenciesSet(constructor). There must be " + TextTools::toString(alphabet->getSize()) + " frequencies.");
  double sum = 0.0;
  for(unsigned int i = 0; i < initFreqs.size(); i++)
  {
    sum += initFreqs[i];
  }
  if(fabs(1-sum) > 0.00000000000001)
  {
    throw Exception("Root frequencies must equal 1.");
  }
  _freq.resize(alphabet->getSize());
  for(unsigned int i = 0; i < alphabet->getSize(); i++)
  {
    _parameters.addParameter(Parameter(prefix + alphabet->intToChar((int)i), initFreqs[i], &Parameter::PROP_CONSTRAINT_IN));
    _freq[i] = initFreqs[i];
  }
}

FullNAFrequenciesSet::FullNAFrequenciesSet(const NucleicAlphabet * alphabet, const string & prefix):
  AbstractFrequenciesSet(alphabet)
{
  _freq.resize(4);
  _parameters.addParameter(Parameter(prefix + "theta" , 0.5, &Parameter::PROP_CONSTRAINT_EX));
  _parameters.addParameter(Parameter(prefix + "theta1", 0.5, &Parameter::PROP_CONSTRAINT_EX));
  _parameters.addParameter(Parameter(prefix + "theta2", 0.5, &Parameter::PROP_CONSTRAINT_EX));
  _freq[0] = _freq[1] = _freq[2] = _freq[3] = 0.25;
}

FullNAFrequenciesSet::FullNAFrequenciesSet(const NucleicAlphabet * alphabet, double theta, double theta1, double theta2, const string & prefix):
  AbstractFrequenciesSet(alphabet)
{
  _freq.resize(4);
  _parameters.addParameter(Parameter(prefix + "theta" , theta , &Parameter::PROP_CONSTRAINT_EX));
  _parameters.addParameter(Parameter(prefix + "theta1", theta1, &Parameter::PROP_CONSTRAINT_EX));
  _parameters.addParameter(Parameter(prefix + "theta2", theta2, &Parameter::PROP_CONSTRAINT_EX));
  _freq[0] = theta1 * (1. - theta);
  _freq[1] = (1 - theta2) * theta;
  _freq[2] = theta2 * theta;
  _freq[3] = (1 - theta1) * (1. - theta);
}

void FullNAFrequenciesSet::fireParameterChanged(const ParameterList & pl)
{
  double theta  = _parameters[0]->getValue();
  double theta1 = _parameters[1]->getValue();
  double theta2 = _parameters[2]->getValue();
  _freq[0] = theta1 * (1. - theta);
  _freq[1] = (1 - theta2) * theta;
  _freq[2] = theta2 * theta;
  _freq[3] = (1 - theta1) * (1. - theta);
}

FullProteinFrequenciesSet::FullProteinFrequenciesSet(const ProteicAlphabet * alphabet, const vector<double> & initFreqs, const string & prefix) throw (Exception):
  AbstractFrequenciesSet(alphabet)
{
  if(initFreqs.size() != 20)
    throw Exception("FullProteinFrequenciesSet(constructor). There must be 20 frequencies.");
  double sum = 0.0;
  for(unsigned int i = 0; i < 20; i++)
  {
    sum += initFreqs[i];
  }
  if(fabs(1-sum) > 0.00000000000001)
  {
    throw Exception("Root frequencies must equal 1.");
  }
  _freq = initFreqs;
  double cumFreqs = 1.;
  for(unsigned int i = 1; i < 20; i++)
  {
    double theta = _freq[i] / cumFreqs;
    cumFreqs *= (1. - _freq[i]);
    _parameters.addParameter(Parameter(prefix + "theta" + TextTools::toString(i) , theta, &Parameter::PROP_CONSTRAINT_EX));
  }
}

FullProteinFrequenciesSet::FullProteinFrequenciesSet(const ProteicAlphabet * alphabet, const string & prefix):
  AbstractFrequenciesSet(alphabet)
{
  _freq.resize(20);
  for(unsigned int i = 1; i < 20; i++)
  {
    _parameters.addParameter(Parameter(prefix + "theta" + TextTools::toString(i) , 0.05 / pow(0.95, (double)i - 1.), &Parameter::PROP_CONSTRAINT_EX));
    _freq[i-1] = 0.05;
  }
  _freq[19] = 0.05;
}

void FullProteinFrequenciesSet::fireParameterChanged(const ParameterList & pl)
{
  double cumTheta = 1.;
  for(unsigned int i = 1; i < 20; i++)
  {
    double theta = _parameters[i]->getValue();
    _freq[i-1] = cumTheta * theta;
    cumTheta *= (1. - theta);
  }
  _freq[19] = cumTheta;
}

FixedFrequenciesSet::FixedFrequenciesSet(const Alphabet * alphabet, const vector<double>& initFreqs, const string & prefix):
  AbstractFrequenciesSet(alphabet)
{
  double sum = 0.0;
  for(unsigned int i = 0; i < initFreqs.size(); i++)
  {
    sum += initFreqs[i];
  }
  if(fabs(1-sum) > 0.00000000000001)
  {
    throw Exception("Root frequencies must equal 1.");
  }
  _freq = initFreqs;
}

