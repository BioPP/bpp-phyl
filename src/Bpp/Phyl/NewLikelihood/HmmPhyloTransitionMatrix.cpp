//
// File: HmmPhyloTransitionMatrix.cpp
// Created by: Laurent Guéguen
// Created on: samedi 21 septembre 2013, à 14h 43
//

/*
Copyright or © or Copr. Bio++Development Team, (November 16, 2004)

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

#include "HmmPhyloTransitionMatrix.h"

#include <Bpp/Text/TextTools.h>

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>

using namespace bpp;
using namespace std;

// HmmPhyloTransitionMatrix::HmmPhyloTransitionMatrix(size_t size, const string& prefix) :
//   AbstractParametrizable(prefix),
//   vSimplex_(),
//   intAlph_(new HmmIntegerAlphabet((unsigned int)size)),
//   pij_(size, size),
//   tmpmat_(size, size),
//   eqFreq_(size),
//   upToDate_(false)
// {
//   for (size_t i=0; i<size; i++)
//   {
//     vSimplex_.push_back(Simplex(size,1));
//     ParameterList pl=vSimplex_[i].getParameters();
//     for (size_t j=0; j<pl.size(); j++)
//     {
//       Parameter* p=pl[j].clone();
//       p->setName(TextTools::toString(i+1)+"."+p->getName());
//       addParameter_(p);
//     }
//   }
// }

HmmPhyloTransitionMatrix::HmmPhyloTransitionMatrix(const HmmProcessAlphabet* procalph, const string& prefix) :
  AbstractParametrizable(prefix),
  vSimplex_(),
  procAlph_(procalph),
  pij_((size_t)procalph->getNumberOfStates(), (size_t)procalph->getNumberOfStates()),
  tmpmat_((size_t)procalph->getNumberOfStates(), (size_t)procalph->getNumberOfStates()),
  eqFreq_((size_t)procalph->getNumberOfStates()),
  upToDate_(false)
{
  size_t size=(size_t)procalph->getNumberOfStates();

  for (size_t i=0; i<size; i++)
    {
      vSimplex_.push_back(Simplex(size,1));
      ParameterList pl=vSimplex_[i].getParameters();
      for (size_t j=0; j<pl.size(); j++)
        {
          Parameter* p=pl[j].clone();
          p->setName(TextTools::toString(i+1)+"."+p->getName());
          addParameter_(p);
        }
    }
}

HmmPhyloTransitionMatrix::HmmPhyloTransitionMatrix(const HmmPhyloTransitionMatrix& hptm) :
  AbstractParametrizable(hptm),
  vSimplex_(hptm.vSimplex_),
  procAlph_(hptm.procAlph_),
  pij_(hptm.pij_),
  tmpmat_(hptm.tmpmat_),
  eqFreq_(hptm.eqFreq_),
  upToDate_(hptm.upToDate_)
{
}

HmmPhyloTransitionMatrix& HmmPhyloTransitionMatrix::operator=(const HmmPhyloTransitionMatrix& hptm)
{
  AbstractParametrizable::operator=(hptm);
  vSimplex_=hptm.vSimplex_;
  procAlph_=hptm.procAlph_;
  pij_=hptm.pij_;
  tmpmat_=hptm.tmpmat_;
  eqFreq_=hptm.eqFreq_;
  upToDate_=hptm.upToDate_;
  
  return *this;
}

void HmmPhyloTransitionMatrix::setHmmStateAlphabet(const HmmStateAlphabet* stateAlphabet) throw (HmmUnvalidAlphabetException)
{
  if (stateAlphabet==NULL)
    throw HmmUnvalidAlphabetException("Null alphabet in HmmPhyloTransitionMatrix::setHmmStateAlphabet");
  if (dynamic_cast<const HmmProcessAlphabet*>(stateAlphabet)==NULL)
    throw HmmUnvalidAlphabetException("Non Process alphabet in HmmPhyloTransitionMatrix::setHmmStateAlphabet");

  procAlph_=dynamic_cast<const HmmProcessAlphabet*>(stateAlphabet);
}

void HmmPhyloTransitionMatrix::setTransitionProbabilities(const Matrix<double>& mat)
{
  if (mat.getNumberOfRows()!=vSimplex_.size())
    throw BadSizeException("HmmPhyloTransitionMatrix::setTransitionProbabilities: Wrong number of rows in given Matrix", mat.getNumberOfRows(), vSimplex_.size());
  
  ParameterList pl;
  
  for (size_t i=0; i<mat.getNumberOfRows();i++)
  {
    vSimplex_[i].setFrequencies(mat.row(i));
    ParameterList pls=vSimplex_[i].getParameters();
    for (size_t j=0; j<pls.size(); j++)
    {
      Parameter* p=pls[j].clone();
      p->setName(TextTools::toString(i+1)+"."+p->getName());
      pl.addParameter(p);
    }
  }
  
  matchParametersValues(pl);
}


const Matrix<double>& HmmPhyloTransitionMatrix::getPij() const
 {
   if (!upToDate_){
     for (size_t i=0; i<vSimplex_.size(); i++)
       for (size_t j=0; j<vSimplex_[i].dimension(); j++)
         pij_(i,j)=vSimplex_[i].prob(j);
     upToDate_=true;
   }
   
   return pij_;
 }

const std::vector<double>& HmmPhyloTransitionMatrix::getEquilibriumFrequencies() const
{
  size_t salph=getNumberOfStates();
  
  if (!upToDate_){
    pij_=getPij();

    MatrixTools::pow(pij_, 256, tmpmat_);

    for (size_t i = 0; i < salph; i++)
      eqFreq_[i] = tmpmat_(0,i);
    
    upToDate_=true;
  }

  return eqFreq_;
}

void HmmPhyloTransitionMatrix::fireParameterChanged(const ParameterList& parameters)
{
  size_t salph=getNumberOfStates();

  for (size_t i=0; i< salph; i++)
  {
    ParameterList pl=vSimplex_[i].getParameters();

    for (size_t j=0; j<pl.size(); j++)
      pl[j].setValue(getParameterValue(TextTools::toString(i+1)+"."+pl[j].getName()));

    vSimplex_[i].matchParametersValues(pl);
  }
  
  upToDate_=false;
}


