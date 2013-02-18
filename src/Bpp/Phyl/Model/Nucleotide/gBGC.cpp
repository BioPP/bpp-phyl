//
// File: gBGC.cpp
// Created by: Laurent Gueguen
// Created on: lundi 13 février 2012, à 09h 42
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)
   This software is a computer program whose purpose is to provide
   classes for phylogenetic data analysis.

   This software is governed by the CeCILL license under French law and
   abiding by the rules of distribution of free software. You can use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and rights to copy,
   modify and redistribute granted by the license, users are provided
   only with a limited warranty and the software's author, the holder of
   the economic rights, and the successive licensors have only limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading, using, modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean that it is complicated to manipulate, and that also
   therefore means that it is reserved for developers and experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards
   their requirements in conditions enabling the security of their
   systems and/or data to be ensured and, more generally, to use and
   operate it in the same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#include "gBGC.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/EigenValue.h>

/******************************************************************************/

gBGC::gBGC(const NucleicAlphabet* alph, NucleotideSubstitutionModel* const pm, double gamma) :
  AbstractParameterAliasable("gBGC."),
  AbstractSubstitutionModel(alph,"gBGC."),
  pmodel_(pm->clone()),
  nestedPrefix_(pm->getNamespace()),
  gamma_(gamma)
{
  pmodel_->setNamespace("gBGC."+nestedPrefix_);
  pmodel_->enableEigenDecomposition(0);
  addParameters_(pmodel_->getParameters());
  addParameter_(new Parameter("gBGC.gamma", gamma_,new IntervalConstraint(-999, 10, true, true), true));

  updateMatrices();
}

gBGC::gBGC(const gBGC& gbgc) :
  AbstractParameterAliasable(gbgc),
  AbstractSubstitutionModel(gbgc),
  pmodel_(gbgc.pmodel_->clone()),
  nestedPrefix_(gbgc.nestedPrefix_),
  gamma_(gbgc.gamma_)
{
}

gBGC& gBGC::operator=(const gBGC& gbgc)
{
  AbstractParameterAliasable::operator=(gbgc);
  AbstractSubstitutionModel::operator=(gbgc);
  pmodel_ = gbgc.pmodel_->clone();
  nestedPrefix_ = gbgc.nestedPrefix_;
  gamma_=gbgc.gamma_;
  return *this;
}

void gBGC::fireParameterChanged(const ParameterList& parameters)
{
  pmodel_->matchParametersValues(parameters);
  AbstractSubstitutionModel::matchParametersValues(parameters);
  updateMatrices();
}

void gBGC::updateMatrices()
{
  gamma_=getParameterValue("gamma");
  unsigned int i,j;
  // Generator:
  double eg=exp(gamma_);
  for ( i = 0; i < 4; i++)
    for ( j = 0; j < 4; j++)
      generator_(i,j)=pmodel_->Qij(i,j);

  generator_(0,1) *= eg;
  generator_(0,2) *= eg;
  generator_(3,1) *= eg;
  generator_(3,2) *= eg;

  generator_(0,0) -= (generator_(0,1)+generator_(0,2))*(1-1/eg);
  generator_(3,3) -= (generator_(3,1)+generator_(3,2))*(1-1/eg);

  // calcul spectral

  EigenValue<double> ev(generator_);
  eigenValues_ = ev.getRealEigenValues();
  
  rightEigenVectors_ = ev.getV();
  MatrixTools::inv(rightEigenVectors_,leftEigenVectors_);

  iEigenValues_ = ev.getImagEigenValues();

  // frequence stationnaire
  double x = 0;
  j = 0;
  while (j < 4){
    if (abs(eigenValues_[j]) < 0.000001 && abs(iEigenValues_[j]) < 0.000001) {
      eigenValues_[j]=0; //to avoid approximation problems in the future
      iEigenValues_[j]=0; //to avoid approximation problems in the future
      for (i = 0; i < 4; i++)
        {
          freq_[i] = leftEigenVectors_(j,i);
          x += freq_[i];
        }
      break;
    }
    j++;
  }

  for (i = 0; i < 4; i++)
    freq_[i] /= x;

  // mise a l'echelle

  x = 0;
  for (i = 0; i < 4; i++)
    x += freq_[i] * generator_(i,i);

  MatrixTools::scale(generator_,-1 / x);

  for (i = 0; i < 4; i++)
    eigenValues_[i] /= -x;
  
  isDiagonalizable_=true;
  for (i = 0; i < size_ && isDiagonalizable_; i++)
    if (abs(iEigenValues_[i])> NumConstants::SMALL()){
      isDiagonalizable_=false;
      break;
    }

  // and the exchangeability_
  for ( i = 0; i < size_; i++)
    for ( j = 0; j < size_; j++)
      exchangeability_(i,j) = generator_(i,j) / freq_[j];

}

void gBGC::setNamespace(const std::string& prefix)
{
  AbstractSubstitutionModel::setNamespace(prefix);
  // We also need to update the namespace of the nested model:
  pmodel_->setNamespace(prefix + nestedPrefix_);
}


std::string gBGC::getName() const
{
  return "gBGC(model=" + pmodel_->getName()+")";
}

