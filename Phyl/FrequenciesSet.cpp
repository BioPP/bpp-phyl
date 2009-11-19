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

//From NumCalc:
#include <NumCalc/NumConstants.h>
#include <NumCalc/VectorTools.h>

/////////////////////////////////////////
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

//////////////////////////////
// FullFrequenciesSet

FullFrequenciesSet::FullFrequenciesSet(const Alphabet* alphabet):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Full.")
{
  unsigned int size = alphabet->getSize();

  for (unsigned int i = 0; i < alphabet->getSize() - 1; i++)
  {
    Parameter p("Full.theta_" + TextTools::toString(i + 1), 1. / (size - i), &Parameter::PROP_CONSTRAINT_IN);
    addParameter_(p);
    getFreq_(i) = 1. / size;
  }
  unsigned int i = alphabet->getSize() - 1;
  getFreq_(i) = 1. / size;
}

FullFrequenciesSet::FullFrequenciesSet(const Alphabet* alphabet, const vector<double>& initFreqs) throw (Exception):
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
    Parameter p("Full.theta_" + TextTools::toString(i + 1), initFreqs[i] / y, &Parameter::PROP_CONSTRAINT_IN);
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
  for (i = 0; i < alphabet->getSize()-1; i++)
  {
    getFreq_(i) = getParameter_("theta_" + TextTools::toString(i + 1)).getValue() * y;
    y *= 1 - getParameter_("theta_" + TextTools::toString(i + 1)).getValue();
  }
  
  i = alphabet->getSize() - 1;
  getFreq_(i) = y;
}



//////////////////////////////
// FullCodonFrequenciesSet

FullCodonFrequenciesSet::FullCodonFrequenciesSet(const CodonAlphabet* alphabet):
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
      j++;
      Parameter p("Full.theta_" + TextTools::toString(i+1), 1. / (size - j), &Parameter::PROP_CONSTRAINT_IN);
      addParameter_(p);
      getFreq_(i) = 1. / size;
    }
  }
  unsigned int i = alphabet->getSize() - 1;
  getFreq_(i) = (alphabet->isStop(i)) ? 0 : 1. / size;
}


FullCodonFrequenciesSet::FullCodonFrequenciesSet(const CodonAlphabet* alphabet, const vector<double>& initFreqs) throw (Exception):
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
      Parameter p("Full.theta_" + TextTools::toString(i+1), initFreqs[i] / y, &Parameter::PROP_CONSTRAINT_IN);
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
      getParameter_("theta_" + TextTools::toString(i+1)).setValue(frequencies[i] / y);
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
  for (i = 0; i < alphabet->getSize()-1; i++)
  {
    if (!(alphabet->isStop(i)))
      getFreq_(i) = getParameter_("theta_" + TextTools::toString(i+1)).getValue() * y;
    y *= 1 - getParameter_("theta_" + TextTools::toString(i+1)).getValue();
  }
  
  i = alphabet->getSize() - 1;
  getFreq_(i) = (alphabet->isStop(i)) ? 0 : y;
}



/////////////////////////////////////////
// FullNAFrequenciesSet


FullNAFrequenciesSet::FullNAFrequenciesSet(const NucleicAlphabet * alphabet):
  AbstractFrequenciesSet(4, alphabet, "Full.")
{
  Parameter thetaP("Full.theta" , 0.5, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(thetaP);
  Parameter theta1P("Full.theta_1", 0.5, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(theta1P);
  Parameter theta2P("Full.theta_2", 0.5, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(theta2P);
  getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
}

FullNAFrequenciesSet::FullNAFrequenciesSet(const NucleicAlphabet* alphabet, double theta, double theta_1, double theta_2):
  AbstractFrequenciesSet(4, alphabet, "Full.")
{
  Parameter thetaP("Full.theta" , theta, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(thetaP);
  Parameter theta1P("Full.theta_1", theta_1, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(theta1P);
  Parameter theta2P("Full.theta_2", theta_2, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(theta2P);
  getFreq_(0) = theta_1 * (1. - theta);
  getFreq_(1) = (1 - theta_2) * theta;
  getFreq_(2) = theta_2 * theta;
  getFreq_(3) = (1 - theta_1) * (1. - theta);
}

void FullNAFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if (frequencies.size() != 4) throw DimensionException(" FullNAFrequenciesSet::setFrequencies", frequencies.size(), 4);
  double sum = 0.0;
  for (unsigned int i = 0; i < 4; i++)
    sum += frequencies[i];
  if (fabs(1. - sum) > NumConstants::SMALL)
    throw Exception("FullNAFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  double theta = frequencies[1] + frequencies[2];
  getParameter_(0).setValue(theta);
  getParameter_(1).setValue(frequencies[0] / (1 - theta));
  getParameter_(2).setValue(frequencies[2] / theta);

  setFrequencies_(frequencies);
}

void FullNAFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
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
    sum += frequencies[i];
  if(fabs(1. - sum) > NumConstants::SMALL)
    throw Exception("GCFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  double theta = frequencies[1] + frequencies[2];
  //We set everything in one shot here:
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

FixedFrequenciesSet::FixedFrequenciesSet(const Alphabet* alphabet, const vector<double>& initFreqs):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.")
{
  setFrequencies(initFreqs);
}

FixedFrequenciesSet::FixedFrequenciesSet(const Alphabet* alphabet):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.")
{
  unsigned int size = alphabet->getSize();
  
  for(unsigned int i = 0; i < alphabet->getSize(); i++)
    getFreq_(i) = 1. / size;
}

void FixedFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if(frequencies.size() != getAlphabet()->getSize()) throw DimensionException("FixedFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  for(unsigned int i = 0; i < frequencies.size(); i++)
    sum += frequencies[i];
  if(fabs(1. - sum) > 0.000001)
    throw Exception("FixedFrequenciesSet::setFrequencies. Frequencies sum must equal 1 (sum = " + TextTools::toString(sum) + ").");
  setFrequencies_(frequencies);
}

/////////////////////////////////////////////
/// CodonFixedFrequenciesSet

CodonFixedFrequenciesSet::CodonFixedFrequenciesSet(const CodonAlphabet* alphabet, const vector<double>& initFreqs):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.")
{
  setFrequencies(initFreqs);
}

CodonFixedFrequenciesSet::CodonFixedFrequenciesSet(const CodonAlphabet* alphabet):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.")
{
  unsigned int size = alphabet->getSize() - alphabet->numberOfStopCodons();
  
  for(unsigned int i = 0; i < alphabet->getSize(); i++)
    getFreq_(i) = (alphabet->isStop(i)) ? 0 : 1. / size;
}

void CodonFixedFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if(frequencies.size() != getAlphabet()->getSize()) throw DimensionException("FixedFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  for(unsigned int i = 0; i < frequencies.size(); i++)
    sum += frequencies[i];
  if(fabs(1. - sum) > 0.000001)
    throw Exception("CodonFixedFrequenciesSet::setFrequencies. Frequencies sum must equal 1 (sum = " + TextTools::toString(sum) + ").");
  setFrequencies_(frequencies);
}

///////////////////////////////////////////////
// IndependentWordFrequenciesSet

IndependentWordFrequenciesSet::IndependentWordFrequenciesSet(const Vector<FrequenciesSet*>& freqvector) : AbstractFrequenciesSet(getSizeFromVector(freqvector),extract_alph(freqvector),"IndependentWord.")
{
  int i,j,k,t;
  int l=freqvector.size();

  unique_AbsFreq=0;
  i=0;
  j=1;
  while (!unique_AbsFreq && i<(l-1)){
    if (freqvector[i]==freqvector[j])
      unique_AbsFreq=1;
    else {
      j++;
      if (j==l){
        i++;
        j=i+1;
      }
    }
  }

  if (!unique_AbsFreq){
    for (i=0;i<l;i++){
      _VFreq.push_back(freqvector[i]);
      _VnestedPrefix.push_back(freqvector[i]->getNamespace());
      _VFreq[i]->setNamespace("IndependentWord."+TextTools::toString(i)+"_"+_VnestedPrefix[i]);
      addParameters_(_VFreq[i]->getParameters());
    }
  }
  else {
    string st="";
    for (i=0;i< l; i++){
      _VFreq.push_back(freqvector[0]);
      _VnestedPrefix.push_back(freqvector[0]->getNamespace());
      st+=TextTools::toString(i);
    }
    _VFreq[0]->setNamespace("IndependentWord."+st+"_"+_VnestedPrefix[0]);
    addParameters_(_VFreq[0]->getParameters());
  }

  vector<double> f[l];
  for (i=0;i<l;i++)
    f[i]=_VFreq[i]->getFrequencies();

  int s=getAlphabet()->getSize();

  for (i=0;i<s;i++){
    j=i;
    getFreq_(j)=1;
    for (k=l-1;k>=0;k--){
      t=_VFreq[k]->getAlphabet()->getSize();
      getFreq_(i)*=f[k][j%t];
      j/=t;
    }
  }

}

IndependentWordFrequenciesSet::IndependentWordFrequenciesSet(FrequenciesSet* pabsfreq, int num) : AbstractFrequenciesSet((int)pow((double) pabsfreq->getAlphabet()->getSize(),num),new WordAlphabet(pabsfreq->getAlphabet(), num),"IndependentWord."), unique_AbsFreq(1)
{
  int i,j,k,t;

  string st="";
  for (i=0;i< num; i++){
    _VFreq.push_back(pabsfreq);
    _VnestedPrefix.push_back(pabsfreq->getNamespace());
    st+=TextTools::toString(i);
  }
  _VFreq[0]->setNamespace("IndependentWord."+st+"_"+_VnestedPrefix[0]);
  addParameters_(_VFreq[0]->getParameters());
  
  vector<double> f[num];
  for (i=0;i<num;i++)
    f[i]=_VFreq[i]->getFrequencies();

  int s=getAlphabet()->getSize();
  
  for (i=0;i<s;i++){
    j=i;
    getFreq_(j)=1;
    for (k=num-1;k>=0;k--){
      t=_VFreq[k]->getAlphabet()->getSize();
      getFreq_(i)*=f[k][j%t];
      j/=t;
    }
  }
}

Alphabet* IndependentWordFrequenciesSet::extract_alph(const Vector<FrequenciesSet*>& freqvector)
{
  int i,l=freqvector.size();
  Vector<const Alphabet*> pval;

  for (i=0;i<l;i++)
    pval.push_back(freqvector[i]->getAlphabet());

  return (new WordAlphabet(pval));
}

unsigned int IndependentWordFrequenciesSet::getSizeFromVector(const Vector<FrequenciesSet*>& freqvector)
{
  unsigned int s = 1;
  unsigned int l = freqvector.size();

  for (unsigned int i = 0; i < l; i++)
    s *= freqvector[i]->getAlphabet()->getSize();

  return s;
}

IndependentWordFrequenciesSet::IndependentWordFrequenciesSet(const IndependentWordFrequenciesSet& iwfs) :
  AbstractFrequenciesSet(iwfs.getNumberOfFrequencies(), new WordAlphabet(*dynamic_cast<const WordAlphabet*>(iwfs.getAlphabet())), iwfs.getNamespace())
{
  unsigned int l = iwfs._VFreq.size();

  FrequenciesSet* pAFS = (unique_AbsFreq ? iwfs._VFreq[0]->clone() : 0);
  
  for (unsigned i = 0; i < l; i++)
  {
    _VFreq.push_back(unique_AbsFreq?pAFS:iwfs._VFreq[i]->clone());
    _VnestedPrefix.push_back(iwfs._VnestedPrefix[i]);
  }
}

void IndependentWordFrequenciesSet::fireParameterChanged(const ParameterList& pl)
{
  unsigned int l = _VFreq.size();

  bool f = 0;
  if (unique_AbsFreq)
    f = _VFreq[0]->matchParametersValues(pl);
  else
    for (unsigned int i = 0; i < l; i++)
      f |= _VFreq[i]->matchParametersValues(pl);

  if (f)
    updateFrequencies();
}

void IndependentWordFrequenciesSet::updateFrequencies()
{
  unsigned int l = _VFreq.size();
  unsigned int s = getAlphabet()->getSize();
  vector<double> f[l];

  unsigned int i, p, t, i2;
// [Julien, 30/10/09] : We should not need that, has frequencies are automatically updated when a parameter is changed...
//  if (unique_AbsFreq)
//    _VFreq[0]->updateFrequencies();
//  else
//    for (i = 0; i < l; i++)
//      _VFreq[i]->updateFrequencies();

  for (i = 0; i < l; i++)
    f[i] = _VFreq[i]->getFrequencies();
  
  for (i = 0; i < s; i++)
  {
    i2 = i;
    getFreq_(i2) = 1;
    for (p = l - 1; p >= 0; p--)
    {
      t = _VFreq[p]->getAlphabet()->getSize();
      getFreq_(i) *= f[p][i2%t];
      i2 /= t;
    }
  }
}

void IndependentWordFrequenciesSet::setFrequencies(const vector<double>& frequencies) throw (DimensionException, Exception)
{
  if(frequencies.size() != getAlphabet()->getSize()) throw DimensionException("IndependentWordFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  unsigned int size = frequencies.size();
  for (unsigned int i = 0; i < size; i++)
    sum += frequencies[i];
  if (fabs(1.-sum) > 0.000001)
    throw Exception("IndependentWordFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");

  unsigned int d,i,j,k,s,l=_VFreq.size();
  vector<double> freq;

  if (unique_AbsFreq)
  {
    s = _VFreq[0]->getAlphabet()->getSize();
    freq.resize(s);
    for (j = 0; j < s; j++)
      freq[j] = 0;
    d = size;
    for (i = 0; i < l; i++){
      d/=s;
      for (k = 0; k < size; k++)
        freq[(k/d)%s] += frequencies[k];
    }
    for (j = 0; j < s; j++)
      freq[j] /= l;
    _VFreq[0]->setFrequencies(freq);
  }
  else
  {
    d = size;
    for (i = 0; i < l; i++)
    {
      s = _VFreq[i]->getAlphabet()->getSize();
      freq.resize(s);
      d /= s;
      for (j = 0; j < s; j++)
        freq[j] = 0;
      for (k = 0; k < size; k++)
        freq[(k/d) % s] += frequencies[k];
      _VFreq[i]->setFrequencies(freq);
    }
  }

  // updating freq_

  vector<double> f[l];
  for (k = 0; k < l; k++)
    f[k] = _VFreq[k]->getFrequencies();
  
  for (i = 0; i < size; i++)
  {
    j = i;
    getFreq_(j)=1;
    for (k = l - 1; k >= 0; k--)
    {
      s = _VFreq[k]->getAlphabet()->getSize();
      getFreq_(i) *= f[k][j % s];
      j/=s;
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
  for (unsigned int i=0;i< _VFreq.size(); i++)
    _VFreq[i]->setNamespace(prefix +TextTools::toString(i)+"_"+_VnestedPrefix[i]);
}

string IndependentWordFrequenciesSet::getName() const
{
  string s= "IndependentWord : ";
  for (unsigned int i=0;i< _VFreq.size(); i++)
    s+= _VFreq[i]->getName();
  return s;
}
