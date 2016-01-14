//
// File: gBGC.cpp
// Created by: Laurent Gueguen
// Created on: lundi 13 février 2012, à 09h 42
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)
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

gBGC::gBGC(const NucleicAlphabet* alph, NucleotideSubstitutionModel* const pm, double B) :
  AbstractParameterAliasable("gBGC."),
  AbstractNucleotideSubstitutionModel(alph, new CanonicalStateMap(alph, false), "gBGC."),
  model_(pm),
  nestedPrefix_(pm->getNamespace()),
  B_(B)
{
  model_->setNamespace("gBGC."+nestedPrefix_);
  model_->enableEigenDecomposition(0);
  addParameters_(model_->getParameters());
  addParameter_(new Parameter("gBGC.B", B_, new IntervalConstraint(-999, 10, true, true), true));

  updateMatrices();
}

gBGC::gBGC(const gBGC& gbgc) :
  AbstractParameterAliasable(gbgc),
  AbstractNucleotideSubstitutionModel(gbgc),
  model_(gbgc.model_->clone()),
  nestedPrefix_(gbgc.nestedPrefix_),
  B_(gbgc.B_)
{
}

gBGC& gBGC::operator=(const gBGC& gbgc)
{
  AbstractParameterAliasable::operator=(gbgc);
  AbstractSubstitutionModel::operator=(gbgc);
  model_ = auto_ptr<NucleotideSubstitutionModel>(gbgc.model_.get());
  nestedPrefix_ = gbgc.nestedPrefix_;
  B_=gbgc.B_;
  return *this;
}

void gBGC::fireParameterChanged(const ParameterList& parameters)
{
  AbstractSubstitutionModel::fireParameterChanged(parameters);
  model_->matchParametersValues(parameters);
  updateMatrices();
}

void gBGC::updateMatrices()
{
  B_=getParameterValue("B");
  unsigned int i,j;
  // Generator:

  for ( i = 0; i < 4; i++)
    for ( j = 0; j < 4; j++)
      generator_(i,j)=model_->Qij(i,j);

  if (B_!=0)
  {    
    double bp=B_/(1-exp(-B_));
    double bm=B_/(exp(B_)-1);
    
    generator_(0,0) -= (generator_(0,1)+generator_(0,2))*(bp-1);
    generator_(1,1) -= (generator_(1,0)+generator_(1,3))*(bm-1);
    generator_(2,2) -= (generator_(2,0)+generator_(2,3))*(bm-1);
    generator_(3,3) -= (generator_(3,1)+generator_(3,2))*(bp-1);

    generator_(0,1) *= bp;
    generator_(0,2) *= bp;
    generator_(3,1) *= bp;
    generator_(3,2) *= bp;
    generator_(1,0) *= bm;
    generator_(2,0) *= bm;
    generator_(1,3) *= bm;
    generator_(2,3) *= bm;
  }

  if (enableEigenDecomposition())
  {
    // calcul spectral
    
    EigenValue<double> ev(generator_);
    eigenValues_ = ev.getRealEigenValues();
    iEigenValues_ = ev.getImagEigenValues();
  
    rightEigenVectors_ = ev.getV();
    try
    {
      MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
      isNonSingular_ = true;
      isDiagonalizable_ = true;
      
      for (i = 0; i < 4 && isDiagonalizable_; i++)
      {
        if (abs(iEigenValues_[i]) > NumConstants::TINY())
          isDiagonalizable_ = false;
      }

      // frequence stationnaire

      if (isDiagonalizable_)
      {
        size_t nulleigen = 0;
        double val;
        isNonSingular_ = false;
        
        while (nulleigen < 4){
          if (abs(eigenValues_[nulleigen]) < 0.000001 && abs(iEigenValues_[nulleigen]) < 0.000001) {
            val = rightEigenVectors_(0, nulleigen);
            i=1;
            while (i < 4)
            {
              if (abs(rightEigenVectors_(i, nulleigen) - val) > NumConstants::SMALL())
                break;
              i++;
            }
            
            if (i < 4)
              nulleigen++;
            else
            {
              isNonSingular_ = true;
              break;
            }
          }
          else
            nulleigen++;
        }

        if (isNonSingular_)
        {
          eigenValues_[nulleigen] = 0; // to avoid approximation errors on long long branches
          iEigenValues_[nulleigen] = 0; // to avoid approximation errors on long long branches
          
          for (i = 0; i < 4; i++)
            freq_[i] = leftEigenVectors_(nulleigen, i);
          
          val = 0;
          for (i = 0; i < 4; i++)
            val += freq_[i];
          
          for (i = 0; i < 4; i++)
            freq_[i] /= val;
        }
        else
        {
          ApplicationTools::displayMessage("Unable to find eigenvector for eigenvalue 1 in gBGC. Taylor series used instead.");
          isDiagonalizable_ = false;
        }
      }
    }
    catch (ZeroDivisionException& e)
    {
      ApplicationTools::displayMessage("Singularity during diagonalization of gBGC in gBGC. Taylor series used instead.");
      isNonSingular_ = false;
      isDiagonalizable_ = false;
    }

    if (!isNonSingular_)
    {
      double min = generator_(0, 0);
      for (i = 1; i < 4; i++)
      {
        if (min > generator_(i, i))
          min = generator_(i, i);
      }

      MatrixTools::scale(generator_, -1 / min);

      if (vPowGen_.size() == 0)
        vPowGen_.resize(30);

      MatrixTools::getId(4, tmpMat_);    // to compute the equilibrium frequency  (Q+Id)^256
      MatrixTools::add(tmpMat_, generator_);
      MatrixTools::pow(tmpMat_, 4, vPowGen_[0]);

      for (i = 0; i < 4; i++)
        freq_[i] = vPowGen_[0](0, i);

      MatrixTools::getId(4, vPowGen_[0]);
    }

    // mise a l'echelle

    double x = 0;
    for (i = 0; i < 4; i++)
      x += freq_[i] * generator_(i,i);

    MatrixTools::scale(generator_,-1 / x);
    for (i = 0; i < 4; i++)
    {
      eigenValues_[i] /= -x;
      iEigenValues_[i] /= -x;
    }

    if (!isNonSingular_)
      MatrixTools::Taylor(generator_, 30, vPowGen_);

  }
}

void gBGC::setNamespace(const std::string& prefix)
{
  AbstractSubstitutionModel::setNamespace(prefix);
  // We also need to update the namespace of the nested model:
  model_->setNamespace(prefix + nestedPrefix_);
}


std::string gBGC::getName() const
{
  return model_->getName()+"+gBGC";
}

