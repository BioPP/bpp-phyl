//
// File: AbstractWordReversibleSubstitutionModel.cpp
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

#include "AbstractWordReversibleSubstitutionModel.h"

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

AbstractWordReversibleSubstitutionModel::AbstractWordReversibleSubstitutionModel(const Vector<SubstitutionModel*>& modelVector,
                                                                                 const std::string& st) : AbstractReversibleSubstitutionModel(AbstractWordReversibleSubstitutionModel::extract_alph(modelVector),st),   new_alphabet_(1)
{
  enableEigenDecomposition(0);
  int i,j;
  int n=modelVector.size();

  //test whether two models are identical

  bool flag=0;
  i=0;
  j=1;
  while (!flag && i<(n-1)){
    if (modelVector[i]==modelVector[j])
      flag=1;
    else {
      j++;
      if (j==n){
        i++;
        j=i+1;
      }
    }
  }

  if (!flag){
    for (i=0;i< n; i++){
      _VAbsRevMod.push_back(modelVector[i]);
      _VnestedPrefix.push_back(modelVector[i]->getNamespace());
      _VAbsRevMod[i]->setNamespace(st+TextTools::toString(i)+"_"+_VnestedPrefix[i]);
      addParameters_(_VAbsRevMod[i]->getParameters());
    }
  }
  else {
    string t="";
    for (i=0;i< n; i++){
      _VAbsRevMod.push_back(modelVector[0]);
      _VnestedPrefix.push_back(modelVector[0]->getNamespace());
      t+=TextTools::toString(i);
    }
    _VAbsRevMod[0]->setNamespace(st+t+"_"+_VnestedPrefix[0]);
    addParameters_(_VAbsRevMod[0]->getParameters());
  }

  _rate=new double[n];
  for (i=0;i< n; i++)
    _rate[i]=1.0/n;

  _p.resize(getNumberOfStates(),getNumberOfStates());
}

AbstractWordReversibleSubstitutionModel::AbstractWordReversibleSubstitutionModel(const Alphabet* alph,
                                                                                 const std::string& st) : AbstractReversibleSubstitutionModel(alph,st), new_alphabet_(0)

{
  enableEigenDecomposition(0);
  _rate=0;
  _p.resize(getNumberOfStates(),getNumberOfStates());
}

AbstractWordReversibleSubstitutionModel::AbstractWordReversibleSubstitutionModel(SubstitutionModel* pmodel,
                                                                                 unsigned int num,
                                                                                 const std::string& st) : AbstractReversibleSubstitutionModel(new WordAlphabet(pmodel->getAlphabet(), num),st),   new_alphabet_(1)
{
  enableEigenDecomposition(0);
  int i;
  _rate=new double[num];

  string t="";
  for (i=0;i< num; i++){
    _VAbsRevMod.push_back(pmodel);
    _VnestedPrefix.push_back(pmodel->getNamespace());
    _rate[i]=1.0/num;
    t+=TextTools::toString(i);
  }
  
  pmodel->setNamespace(st+t+"_"+_VnestedPrefix[0]);
  addParameters_(pmodel->getParameters());
  
  _p.resize(getNumberOfStates(),getNumberOfStates());
}

AbstractWordReversibleSubstitutionModel::AbstractWordReversibleSubstitutionModel(const AbstractWordReversibleSubstitutionModel& wrsm) : AbstractReversibleSubstitutionModel(wrsm.alphabet_,wrsm.getNamespace()),   new_alphabet_(0)
{
  enableEigenDecomposition(0);
  int i;
  int num=wrsm._VAbsRevMod.size();

  _rate=new double[num];
  for (i=0;i< num; i++){
    _VAbsRevMod.push_back(wrsm._VAbsRevMod[i]);
    _VnestedPrefix.push_back(wrsm._VnestedPrefix[i]);
    _rate[i]=wrsm._rate[i];
  }

  addParameters_(wrsm.getParameters());
  _p.resize(getNumberOfStates(),getNumberOfStates());
}

AbstractWordReversibleSubstitutionModel::~AbstractWordReversibleSubstitutionModel()
{
  _VAbsRevMod.clear();

  if (new_alphabet_)
    delete alphabet_;
  delete _rate;
}

unsigned int AbstractWordReversibleSubstitutionModel::getNumberOfStates() const
{
  return getAlphabet()->getSize();
}

Alphabet* AbstractWordReversibleSubstitutionModel::extract_alph(const Vector<SubstitutionModel*>& modelVector)
{
  int i;

  Vector<const Alphabet*> VAlph;
  
  for (i=0;i< modelVector.size(); i++)
    VAlph.push_back(modelVector[i]->getAlphabet());

  return (new WordAlphabet(VAlph));
}

void AbstractWordReversibleSubstitutionModel::setNamespace(const string& prefix)
{
  AbstractReversibleSubstitutionModel::setNamespace(prefix);

  if (_VAbsRevMod.size()<2 || _VAbsRevMod[0]==_VAbsRevMod[1]){
    string t="";
    for (unsigned int i=0;i< _VAbsRevMod.size(); i++)
      t+=TextTools::toString(i);
    _VAbsRevMod[0]->setNamespace(prefix + t +"_"+_VnestedPrefix[0]);
  }
  else
    for (unsigned int i=0;i< _VAbsRevMod.size(); i++)
      _VAbsRevMod[i]->setNamespace(prefix +TextTools::toString(i)+"_"+_VnestedPrefix[i]);
}

void AbstractWordReversibleSubstitutionModel::fireParameterChanged(const ParameterList &parameters)
{
  for (unsigned int i=0;i< _VAbsRevMod.size(); i++)  
    _VAbsRevMod[i]->matchParametersValues(parameters);
  updateMatrices();
}

/******************************************************************************/

void AbstractWordReversibleSubstitutionModel::updateMatrices()
{
  int i,j,k,p,m;
  string s;
  ParameterList ParL;
  ParameterList::iterator itPar;
  double x;
  int nbmod=_VAbsRevMod.size();

  int salph=getNumberOfStates();

  for (i=0;i<salph;i++){
    freq_[i]=1;
    j=i;
    for (p=nbmod-1;p>=0;p--){
      m=_VAbsRevMod[p]->getNumberOfStates();
      freq_[i]*=_VAbsRevMod[p]->getFrequencies()[j%m];
      j/=m;
    }
  }

  // Generator

  Vector<int> vsize;
  int n,l;
  
  for (k=0;k<nbmod;k++)
    vsize.push_back(_VAbsRevMod[k]->getNumberOfStates());

  for (i=0;i<salph;i++)
    for (j=0;j<salph;j++){
      generator_(i,j)=0;
    }

  RowMatrix<double> gk;
  
  m=1;

  for (k=nbmod-1;k>=0;k--){
    gk=_VAbsRevMod[k]->getGenerator();
    for (i=0;i<vsize[k];i++)  
      for (j=0;j<vsize[k];j++) 
        if (i!=j){
          n=0;
          while (n<salph){ //loop on prefix
            for (l=0;l<m;l++){ //loop on suffix
              generator_(n+i*m+l,n+j*m+l)=gk(i,j)*_rate[k];
            }
            n+=m*vsize[k];
          }
        }
    m*=vsize[k];
  }
  
  for (i=0;i<salph;i++){
    x=0;
    for (j=0;j<salph;j++)
      if (j!=i)
        x+=generator_(i,j);
    generator_(i,i)=-x;
  }
  
  completeMatrices();
  
 // Eigen values:

  if (enableEigenDecomposition()){
    EigenValue<double> ev(generator_);
    eigenValues_=ev.getRealEigenValues();
    rightEigenVectors_ = ev.getV();
    MatrixTools::inv(rightEigenVectors_,leftEigenVectors_);

    vector<double> vi=ev.getImagEigenValues();

    // looking for the 0 eigenvector
    
    x=0;
    j=0;
    while (j<salph){
      if (abs(eigenValues_[j])<0.000001 && abs(vi[j])<0.000001){
        // selecting the 0 eigenvalue that does not match with
        // STOP codons
        i=0;
        while (i<salph){
          if (abs(leftEigenVectors_(j,i))>0.001){
            // i is a word in the eigenspace
            // check that i is not a stop codon
            k=0;
            while (k<salph){
              if ((i!=k) && generator_(k,i)>0.001)
                break;
              k++;
            }
            if (k<salph) // i is not a stop codon
              break;
          }
          i++;
        }
        if (i<salph) // good eigenvector
          break;
      }
      j++;
    }

    if (j>=salph){
      cerr << "AbstractWordReversibleSubstitutionModel::updateMatrices : Problem in eigenspace of" << getName() << endl;
      exit(0);
    }

    // stationnary distribution
    
    for (i=0;i<salph;i++){
      if ((i!=j) && (abs(eigenValues_[i])<0.000001)) // if i is in a STOP codon eigenspace, _freq[i]=0
        freq_[i]=0;
      else
        freq_[i]=leftEigenVectors_(j,i);
      
      x+=freq_[i];
    }

    for (i=0;i<salph;i++)
      freq_[i]/=x;

    // normalization
    
    x=0;
    for (i=0;i<salph;i++)
      x+=freq_[i]*generator_(i,i);
    
    MatrixTools::scale(generator_,-1/x);
    
    for (i=0;i<salph;i++)
      eigenValues_[i]/=-x;

  }
  
}


void AbstractWordReversibleSubstitutionModel::setFreq(map<int, double>& freqs)
{
  map<int, double> freq;
  int nbmod=_VAbsRevMod.size();

  int i,j,s,k,d,size;

  d=size=getNumberOfStates();
  
  for (i=0;i<nbmod;i++){
    freq.clear();
    s=_VAbsRevMod[i]->getAlphabet()->getSize();
    d/=s;
    for (j=0;j<s;j++)
      freq[j]=0;
    for (k=0;k<size;k++)
      freq[(k/d)%s]+=freqs[k];
    _VAbsRevMod[i]->setFreq(freq);
  }

  updateMatrices();
}
