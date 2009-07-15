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

/////////////////////////////////////////
// AbstractFrequenciesSet

AbstractFrequenciesSet::AbstractFrequenciesSet(unsigned int n, const Alphabet * alphabet, const string& prefix) :
  AbstractParametrizable(prefix), alphabet_(alphabet), freq_(n)
{
}

void AbstractFrequenciesSet::setFrequenciesFromMap(map<int,double>& mfreq)
{
  int s=getAlphabet()->getSize();
  vector<double> freq(s);
  double x=0;
  for (unsigned int i=0; i< s; i++){
    freq[i]=(mfreq.find(i)!=mfreq.end())?mfreq[i]:0;
    x+=freq[i];
  }
  for (unsigned int i=0; i< s; i++){
    freq[i]/=x;
  }
  setFrequencies(freq);
}

//////////////////////////////
// FullFrequenciesSet

FullFrequenciesSet::FullFrequenciesSet(const Alphabet * alphabet):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Full.")
{
  const CodonAlphabet* pca= AlphabetTools::isCodonAlphabet(alphabet)?dynamic_cast<const CodonAlphabet*>(alphabet):0;
  unsigned int i, size=alphabet->getSize();
  size-=(pca)?pca->numberOfStopCodons():0;
  
  for(i = 0; i < alphabet->getSize(); i++)
  {
    if (pca && pca->isStop(i)){
      getFreq_(i)=0;
    }
    else {
      Parameter p("Full." + alphabet->intToChar((int)i), 1. / size, &Parameter::PROP_CONSTRAINT_IN);
      addParameter_(p);
      getFreq_(i) = 1. / size;
    }
  }
}

FullFrequenciesSet::FullFrequenciesSet(const Alphabet * alphabet, const vector<double>& initFreqs) throw (Exception):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Full.")
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
    throw Exception("Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  }
  const CodonAlphabet* pca= AlphabetTools::isCodonAlphabet(getAlphabet())?dynamic_cast<const CodonAlphabet*>(getAlphabet()):0;

  for(unsigned int i = 0; i < alphabet->getSize(); i++)
  {
    if (!(pca && pca->isStop(i))){
      Parameter p("Full." + alphabet->intToChar((int)i), initFreqs[i], &Parameter::PROP_CONSTRAINT_IN);
      addParameter_(p);
    }
    getFreq_(i) = initFreqs[i];
  }

}

void FullFrequenciesSet::setFrequencies(const vector<double> & frequencies) throw (DimensionException, Exception)
{
  if(frequencies.size() != getAlphabet()->getSize()) throw DimensionException("FullFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  for(unsigned int i = 0; i < frequencies.size(); i++)
    sum += frequencies[i];
  if(fabs(1.-sum) > 0.000001)
    throw Exception("FullFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  
  const CodonAlphabet* pca= AlphabetTools::isCodonAlphabet(getAlphabet())?dynamic_cast<const CodonAlphabet*>(getAlphabet()):0;

  for(unsigned int i = 0; i < getAlphabet()->getSize(); i++)
    {
      if (!(pca && pca->isStop(i)))
        getParameter_(getAlphabet()->intToChar(i)).setValue(frequencies[i]);
      getFreq_(i) = frequencies[i];
    }
}

void FullFrequenciesSet::fireParameterChanged(const ParameterList & pl)
{
  const CodonAlphabet* pca= AlphabetTools::isCodonAlphabet(getAlphabet())?dynamic_cast<const CodonAlphabet*>(getAlphabet()):0;

  for(unsigned int i = 0; i < getAlphabet()->getSize(); i++)
    {
      if (!(pca && pca->isStop(i)))
        getFreq_(i) = getParameter_(getAlphabet()->intToChar(i)).getValue();
    }
}

/////////////////////////////////////////
// FullNAFrequenciesSet


FullNAFrequenciesSet::FullNAFrequenciesSet(const NucleicAlphabet * alphabet):
  AbstractFrequenciesSet(4, alphabet, "FullNA.")
{
  Parameter thetaP("FullNA.theta" , 0.5, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(thetaP);
  Parameter theta1P("FullNA.theta1", 0.5, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(theta1P);
  Parameter theta2P("FullNA.theta2", 0.5, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(theta2P);
  getFreq_(0) = getFreq_(1) = getFreq_(2) = getFreq_(3) = 0.25;
}

FullNAFrequenciesSet::FullNAFrequenciesSet(const NucleicAlphabet * alphabet, double theta, double theta1, double theta2):
  AbstractFrequenciesSet(4, alphabet, "FullNA.")
{
  Parameter thetaP("FullNA.theta" , theta, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(thetaP);
  Parameter theta1P("FullNA.theta1", theta1, &Parameter::PROP_CONSTRAINT_EX);
  addParameter_(theta1P);
  Parameter theta2P("FullNA.theta2", theta2, &Parameter::PROP_CONSTRAINT_EX);
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
  getFreq_(3) = (1 - theta1) * (1. - theta);
}

///////////////////////////////////////////
// GCFrequenciesSet

void GCFrequenciesSet::fireParameterChanged(const ParameterList & pl)
{
  double theta = getParameter_(0).getValue();
  getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
  getFreq_(1) = getFreq_(2) = theta / 2.;
}

void GCFrequenciesSet::setFrequencies(const vector<double> & frequencies) throw (DimensionException, Exception)
{
  if(frequencies.size() != 4) throw DimensionException("GCFrequenciesSet::setFrequencies", frequencies.size(), 4);
  double sum = 0.0;
  for(unsigned int i = 0; i < 4; i++)
    sum += frequencies[i];
  if(fabs(1.-sum) > 0.000001)
    throw Exception("GCFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");
  double theta=frequencies[1] + frequencies[2];
  getParameter_(0).setValue(theta);
  getFreq_(0) = getFreq_(3) = (1. - theta) / 2.;
  getFreq_(1) = getFreq_(2) = theta / 2.;

}


///////////////////////////////////////////
// FullProteinFrequenciesSet

FullProteinFrequenciesSet::FullProteinFrequenciesSet(const ProteicAlphabet * alphabet, const vector<double> & initFreqs) throw (Exception):
  AbstractFrequenciesSet(20, alphabet, "FullProtein.")
{
  if(initFreqs.size() != 20)
    throw Exception("FullProteinFrequenciesSet(constructor). There must be 20 frequencies.");
  for(unsigned int i = 1; i < 20; i++)
  {
    Parameter p("FullProtein.theta" + TextTools::toString(i) , 0.05 / (1-0.05*((double)i - 1.)), &Parameter::PROP_CONSTRAINT_EX);
    addParameter_(p);
  }
  setFrequencies(initFreqs);
}

FullProteinFrequenciesSet::FullProteinFrequenciesSet(const ProteicAlphabet * alphabet):
  AbstractFrequenciesSet(20, alphabet, "FullProtein.")
{
  for(unsigned int i = 1; i < 20; i++)
  {
    Parameter p("FullProtein.theta" + TextTools::toString(i) , 0.05 / (1-0.05*((double)i - 1.)), &Parameter::PROP_CONSTRAINT_EX);
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

/////////////////////////////////////////////
/// FixedFrequenciesSet

FixedFrequenciesSet::FixedFrequenciesSet(const Alphabet * alphabet, const vector<double>& initFreqs):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.")
{
  setFrequencies(initFreqs);
}

FixedFrequenciesSet::FixedFrequenciesSet(const Alphabet * alphabet):
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.")
{
  const CodonAlphabet* pca= AlphabetTools::isCodonAlphabet(alphabet)?dynamic_cast<const CodonAlphabet*>(alphabet):0;
  unsigned int i, size=alphabet->getSize();
  size-=(pca)?pca->numberOfStopCodons():0;
  
  for(i = 0; i < alphabet->getSize(); i++)
    {
      if (pca && pca->isStop(i)){
        getFreq_(i)=0;
      }
      else {
        getFreq_(i) = 1. / size;
      }
    }
}

///////////////////////////////////////////////
// IndependentWordFrequenciesSet

IndependentWordFrequenciesSet::IndependentWordFrequenciesSet(const Vector<AbstractFrequenciesSet*>& freqvector) : AbstractFrequenciesSet(getSize(freqvector),extract_alph(freqvector),"IndFreq.")
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
      _VAbsFreq.push_back(freqvector[i]);
      _VnestedPrefix.push_back(freqvector[i]->getNamespace());
      _VAbsFreq[i]->setNamespace("IndFreq."+TextTools::toString(i)+"_"+_VnestedPrefix[i]);
      addParameters_(_VAbsFreq[i]->getParameters());
    }
  }
  else {
    string st="";
    for (i=0;i< l; i++){
      _VAbsFreq.push_back(freqvector[0]);
      _VnestedPrefix.push_back(freqvector[0]->getNamespace());
      st+=TextTools::toString(i);
    }
    _VAbsFreq[0]->setNamespace("IndFreq."+st+"_"+_VnestedPrefix[0]);
    addParameters_(_VAbsFreq[0]->getParameters());
  }

  vector<double> f[l];
  for (i=0;i<l;i++)
    f[i]=_VAbsFreq[i]->getFrequencies();

  int s=getAlphabet()->getSize();

  for (i=0;i<s;i++){
    j=i;
    getFreq_(j)=1;
    for (k=l-1;k>=0;k--){
      t=_VAbsFreq[k]->getAlphabet()->getSize();
      getFreq_(i)*=f[k][j%t];
      j/=t;
    }
  }

}

IndependentWordFrequenciesSet::IndependentWordFrequenciesSet(AbstractFrequenciesSet* pabsfreq, int num) : AbstractFrequenciesSet((int)pow((double) pabsfreq->getAlphabet()->getSize(),num),new WordAlphabet(pabsfreq->getAlphabet(), num),"IndFreq."), unique_AbsFreq(1)
{
  int i,j,k,t;

  string st="";
  for (i=0;i< num; i++){
    _VAbsFreq.push_back(pabsfreq);
    _VnestedPrefix.push_back(pabsfreq->getNamespace());
    st+=TextTools::toString(i);
  }
  _VAbsFreq[0]->setNamespace("IndFreq."+st+"_"+_VnestedPrefix[0]);
  addParameters_(_VAbsFreq[0]->getParameters());
  
  vector<double> f[num];
  for (i=0;i<num;i++)
    f[i]=_VAbsFreq[i]->getFrequencies();

  int s=getAlphabet()->getSize();
  
  for (i=0;i<s;i++){
    j=i;
    getFreq_(j)=1;
    for (k=num-1;k>=0;k--){
      t=_VAbsFreq[k]->getAlphabet()->getSize();
      getFreq_(i)*=f[k][j%t];
      j/=t;
    }
  }
}

Alphabet* IndependentWordFrequenciesSet::extract_alph(const Vector<AbstractFrequenciesSet*>& freqvector)
{
  int i,l=freqvector.size();
  Vector<const Alphabet*> pval;

  for (i=0;i<l;i++)
    pval.push_back(freqvector[i]->getAlphabet());

  return (new WordAlphabet(pval));
}

unsigned int IndependentWordFrequenciesSet::getSize(const Vector<AbstractFrequenciesSet*>& freqvector)
{
  int s=1;
  int i,l=freqvector.size();

  for (i=0;i<l;i++)
    s*=freqvector[i]->getAlphabet()->getSize();

  return s;
}

IndependentWordFrequenciesSet::IndependentWordFrequenciesSet(const IndependentWordFrequenciesSet& iwfs) : AbstractFrequenciesSet(iwfs), unique_AbsFreq(iwfs.unique_AbsFreq)
{
  int i,l=iwfs._VAbsFreq.size();
  
  for (i=0;i<l;i++){
    _VAbsFreq.push_back(iwfs._VAbsFreq[i]);
    _VnestedPrefix.push_back(iwfs._VnestedPrefix[i]);
  }
}

void IndependentWordFrequenciesSet::fireParameterChanged(const ParameterList& pl)
{
  int l=_VAbsFreq.size();
  
  int i,p,t,i2;
  if (unique_AbsFreq)
    _VAbsFreq[0]->matchParametersValues(pl);
  else
    for (i=0;i<l;i++)
      _VAbsFreq[i]->matchParametersValues(pl);

  int s=getAlphabet()->getSize();
  vector<double> f[l];

  for (i=0;i<l;i++)
    f[i]=_VAbsFreq[i]->getFrequencies();
  
  for (i=0;i<s;i++){
    i2=i;
    getFreq_(i2)=1;
    for (p=l-1;p>=0;p--){
      t=_VAbsFreq[p]->getAlphabet()->getSize();
      getFreq_(i)*=f[p][i2%t];
      i2/=t;
    }
  }
  
}

void IndependentWordFrequenciesSet::setFrequencies(const vector<double> & frequencies) throw (DimensionException, Exception)
{
  if(frequencies.size() != getAlphabet()->getSize()) throw DimensionException("IndependentWordFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  double sum = 0.0;
  int size=frequencies.size();
  for(unsigned int i = 0; i < size; i++)
    sum += frequencies[i];
  if(fabs(1.-sum) > 0.000001)
    throw Exception("IndependentWordFrequenciesSet::setFrequencies. Frequencies must equal 1 (sum = " + TextTools::toString(sum) + ").");

  int d,i,j,k,s,l=_VAbsFreq.size();
  vector<double> freq;

  if (unique_AbsFreq){
    s=_VAbsFreq[0]->getAlphabet()->getSize();
    freq.resize(s);
    for (j=0;j<s;j++)
      freq[j]=0;
    d=size;
    for (i=0;i<l;i++){
      d/=s;
      for (k=0;k<size;k++)
        freq[(k/d)%s]+=frequencies[k];
    }
    for (j=0;j<s;j++)
      freq[j]/=l;
    _VAbsFreq[0]->setFrequencies(freq);
  }
  else {
    d=size;
    for (i=0;i<l;i++){
      s=_VAbsFreq[i]->getAlphabet()->getSize();
      freq.resize(s);
      d/=s;
      for (j=0;j<s;j++)
        freq[j]=0;
      for (k=0;k<size;k++)
        freq[(k/d)%s]+=frequencies[k];
      _VAbsFreq[i]->setFrequencies(freq);
    }
  }

  // updating freq_

  vector<double> f[l];
  for (k=0;k<l;k++)
    f[k]=_VAbsFreq[k]->getFrequencies();
  
  for (i=0;i<size;i++){
    j=i;
    getFreq_(j)=1;
    for (k=l-1;k>=0;k--){
      s=_VAbsFreq[k]->getAlphabet()->getSize();
      getFreq_(i)*=f[k][j%s];
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

  for (unsigned int i=0;i< _VAbsFreq.size(); i++)
    _VAbsFreq[i]->setNamespace(prefix +TextTools::toString(i)+"_"+_VnestedPrefix[i]);
}

string IndependentWordFrequenciesSet::getName() const
{
  string s= "IndependentWord : ";
  for (unsigned int i=0;i< _VAbsFreq.size(); i++)
    s+= _VAbsFreq[i]->getName();
  return s;
}
