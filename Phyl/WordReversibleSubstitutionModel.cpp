//
// File: WordReversibleSubstitutionModel.cpp
// Created by:  Laurent Gueguen
// Created on: Jan 2009
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)
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

#include "WordReversibleSubstitutionModel.h"
#include "FrequenciesSet.h"

// From SeqLib:
#include <Seq/SequenceContainerTools.h>
#include <Seq/WordAlphabet.h>

// From NumCalc:
#include <NumCalc/MatrixTools.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/EigenValue.h>

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/******************************************************************************/

WordReversibleSubstitutionModel::WordReversibleSubstitutionModel(const std::vector<SubstitutionModel*>& modelVector,
                                                                 const std::string& st) : AbstractWordReversibleSubstitutionModel(modelVector, (st=="")?"Word.":st)
{
  int i,nbmod=_VSubMod.size();
  
  // relative rates
  for (i=0; i< nbmod-1; i++){
    addParameter_(Parameter("Word.relrate"+TextTools::toString(i) , 1.0/(nbmod-i),&Parameter::PROP_CONSTRAINT_EX));
  }

  WordReversibleSubstitutionModel::updateMatrices();
}

WordReversibleSubstitutionModel::WordReversibleSubstitutionModel(const Alphabet* alph,
                                                                 const std::string& st) : AbstractWordReversibleSubstitutionModel(alph, (st=="")?"Word.":st)

{
}

WordReversibleSubstitutionModel::WordReversibleSubstitutionModel(SubstitutionModel* pmodel, unsigned int num, const std::string& st) : AbstractWordReversibleSubstitutionModel(pmodel, num,  (st=="")?"Word.":st)
{
  unsigned int i;
  
  // relative rates
  for (i = 0; i < num - 1; i++)
  {
    addParameter_(Parameter("Word.relrate" + TextTools::toString(i) , 1.0 / (num -i ), &Parameter::PROP_CONSTRAINT_EX));
  }

  WordReversibleSubstitutionModel::updateMatrices();
}

void WordReversibleSubstitutionModel::updateMatrices()
{
  int i,k, nbmod=_VSubMod.size();
  double x;
  for (k=nbmod-1;k>=0;k--){
    x=1.0;
    for (i=0;i<k;i++)
      x*=1-getParameterValue("relrate"+TextTools::toString(i));
    if (k!=nbmod-1)
      x*=getParameterValue("relrate"+TextTools::toString(k));
    _rate[k]=x;
  }

  AbstractWordReversibleSubstitutionModel::updateMatrices();
}

void WordReversibleSubstitutionModel::completeMatrices()
{
  int nbmod=_VSubMod.size();
  int i,p,j,m;
  int salph=getAlphabet()->getSize();
  
  // _freq for this generator
  
  for (i=0;i<salph;i++){
    freq_[i]=1;
    j=i;
    for (p=nbmod-1;p>=0;p--){
      m=_VSubMod[p]->getNumberOfStates();
      freq_[i]*=_VSubMod[p]->getFrequencies()[j%m];
      j/=m;
    }
  }
}

double WordReversibleSubstitutionModel::Pij_t(unsigned int i, unsigned int j, double d) const
{
  double x = 1;
  unsigned int nbmod = _VSubMod.size();
  unsigned int t;
  int p;

  unsigned int i2 = i;
  unsigned int j2 = j;

  for (p = nbmod - 1; p >= 0; p--)
  {
    t = _VSubMod[p]->getNumberOfStates();
    x *= _VSubMod[p]->Pij_t(i2%t, j2%t, d*_rate[p]);
    i2 /= t;
    j2 /= t;
  }
  
  return(x);
}

const RowMatrix<double>& WordReversibleSubstitutionModel::getPij_t(double d) const
{
  unsigned int nbStates = getNumberOfStates();
  unsigned int i, j;

  for (i = 0; i < nbStates; i++)
  {
    for (j = 0; j < nbStates; j++)
    {  
      _p(i,j) = Pij_t(i,j,d);
    }
  }
  return _p;
}

double WordReversibleSubstitutionModel::dPij_dt(unsigned int i, unsigned int j, double d) const
{
  double r, x;
   int nbmod = _VSubMod.size();
   int t;
  int p,q;
  
   int i2 = i;
   int j2 = j;

  r = 0;
  for (q = 0; q < nbmod; q++)
  {
    i2 = i;
    j2 = j;
    x = 1;
    for (p = nbmod - 1; p >= 0; p--)
    {
      t=_VSubMod[p]->getNumberOfStates();
      if (q!=p)
        x*=_VSubMod[p]->Pij_t(i2%t,j2%t,d*_rate[p]);
      else
        x*=_rate[p]*_VSubMod[p]->dPij_dt(i2%t,j2%t,d*_rate[p]);
      i2/=t;
      j2/=t;
    }
    r+=x;
  }
  return(r);
}

const RowMatrix<double>& WordReversibleSubstitutionModel::getdPij_dt(double d) const
{
  unsigned int nbetats = getNumberOfStates();
  unsigned int i, j;

  for (i = 0; i < nbetats; i++)
    for (j = 0; j < nbetats; j++){
      _p(i,j) = dPij_dt(i,j,d);
    }
  return _p;
}

double WordReversibleSubstitutionModel::d2Pij_dt2(unsigned int i, unsigned int j, double d) const
{
  double r, x;
   int nbmod = _VSubMod.size();
   int b, q, t;
  int p;
  
   int i2 = i;
   int j2 = j;

  r = 0;

  for (q = 1; q < nbmod; q++)
  {
    for (b = 0; b < q; b++)
    {
      x = 1;
      i2 = i;
      j2 = j;
      for (p = nbmod - 1; p >= 0; p--)
      {
        t = _VSubMod[p]->getNumberOfStates();
        if (p == q) 
          x *= _rate[p] * _VSubMod[p]->dPij_dt(i2%t, j2%t, d*_rate[p]);
        else if (p==b)
          x *= _rate[p] * _VSubMod[p]->dPij_dt(i2%t, j2%t, d*_rate[p]);
        else
          x *= _VSubMod[p]->Pij_t(i2%t, j2%t, d*_rate[p]);
        
        i2 /= t;
        j2 /= t;
      }
      r += x;
    }
  }

  r *= 2;
  
  for (q = 0; q < nbmod; q++)
  {
    x = 1;
    i2 = i;
    j2 = j;
    for (p = nbmod - 1; p >= 0; p--)
    {
      t = _VSubMod[p]->getNumberOfStates();
      if (q != p)
        x *= _VSubMod[p]->Pij_t(i2%t, j2%t, d*_rate[p]);
      else
        x *= _rate[p] * _rate[p] * _VSubMod[p]->d2Pij_dt2(i2%t, j2%t, d * _rate[p]);
      
      i2 /= t;
      j2 /= t;
    }
    r += x;
  }
  
  return(r);
}

const RowMatrix<double>& WordReversibleSubstitutionModel::getd2Pij_dt2(double d) const
{
  unsigned int nbetats = getNumberOfStates();
  unsigned int i,j;

  for (i = 0; i < nbetats; i++)
    for (j = 0; j < nbetats; j++)
      _p(i,j) = Pij_t(i,j,d);

  return _p;
}


string WordReversibleSubstitutionModel::getName() const
{
  unsigned int nbmod = _VSubMod.size();
  string s = "WordReversibleSubstitutionModel model: " + _VSubMod[0]->getName();
  for (unsigned int i = 1; i < nbmod - 1; i++)
    s += " " + _VSubMod[i]->getName();

  return s;
}
    
  
