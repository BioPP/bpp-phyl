// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "gBGC.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/EigenValue.h>

/******************************************************************************/

gBGC::gBGC(
    shared_ptr<const NucleicAlphabet> alph,
    unique_ptr<NucleotideSubstitutionModelInterface> pm,
    double B) :
  AbstractParameterAliasable("gBGC."),
  AbstractNucleotideSubstitutionModel(alph, pm->getStateMap(), "gBGC."),
  model_(std::move(pm)),
  nestedPrefix_(model_->getNamespace()),
  B_(B)
{
  model_->setNamespace("gBGC." + nestedPrefix_);
  model_->enableEigenDecomposition(0);
  model_->computeFrequencies(false);

  addParameters_(model_->getParameters());
  addParameter_(new Parameter("gBGC.B", B_, std::make_shared<IntervalConstraint>(-999, 10, true, true)));

  computeFrequencies(true);
  updateMatrices_();
}

gBGC::gBGC(const gBGC& gbgc) :
  AbstractParameterAliasable(gbgc),
  AbstractNucleotideSubstitutionModel(gbgc),
  model_(gbgc.model_->clone()),
  nestedPrefix_(gbgc.nestedPrefix_),
  B_(gbgc.B_)
{}

gBGC& gBGC::operator=(const gBGC& gbgc)
{
  AbstractParameterAliasable::operator=(gbgc);
  AbstractSubstitutionModel::operator=(gbgc);
  model_.reset(gbgc.model_->clone());
  nestedPrefix_ = gbgc.nestedPrefix_;
  B_ = gbgc.B_;
  return *this;
}

void gBGC::fireParameterChanged(const ParameterList& parameters)
{
  AbstractSubstitutionModel::fireParameterChanged(parameters);
  model_->matchParametersValues(parameters);
  updateMatrices_();
}

void gBGC::updateMatrices_()
{
  B_ = getParameterValue("B");
  unsigned int i, j;
  // Generator:

  for (i = 0; i < 4; ++i)
  {
    for (j = 0; j < 4; ++j)
    {
      generator_(i, j) = model_->Qij(i, j);
    }
  }

  if (B_ != 0)
  {
    double bp = B_ / (1 - exp(-B_));
    double bm = B_ / (exp(B_) - 1);

    generator_(0, 0) -= (generator_(0, 1) + generator_(0, 2)) * (bp - 1);
    generator_(1, 1) -= (generator_(1, 0) + generator_(1, 3)) * (bm - 1);
    generator_(2, 2) -= (generator_(2, 0) + generator_(2, 3)) * (bm - 1);
    generator_(3, 3) -= (generator_(3, 1) + generator_(3, 2)) * (bp - 1);

    generator_(0, 1) *= bp;
    generator_(0, 2) *= bp;
    generator_(3, 1) *= bp;
    generator_(3, 2) *= bp;
    generator_(1, 0) *= bm;
    generator_(2, 0) *= bm;
    generator_(1, 3) *= bm;
    generator_(2, 3) *= bm;
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

        while (nulleigen < 4)
        {
          if (abs(eigenValues_[nulleigen]) < 0.000001 && abs(iEigenValues_[nulleigen]) < 0.000001)
          {
            val = rightEigenVectors_(0, nulleigen);
            i = 1;
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
          {
            freq_[i] = leftEigenVectors_(nulleigen, i);
          }

          val = 0;
          for (i = 0; i < 4; i++)
          {
            val += freq_[i];
          }

          for (i = 0; i < 4; i++)
          {
            freq_[i] /= val;
          }
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

      setScale(-1 / min);

      if (vPowGen_.size() == 0)
        vPowGen_.resize(30);

      MatrixTools::getId(4, tmpMat_);    // to compute the equilibrium frequency  (Q+Id)^256
      MatrixTools::add(tmpMat_, generator_);
      MatrixTools::pow(tmpMat_, 4, vPowGen_[0]);

      for (i = 0; i < 4; ++i)
      {
        freq_[i] = vPowGen_[0](0, i);
      }

      MatrixTools::getId(4, vPowGen_[0]);
    }

    // mise a l'echelle

    normalize();


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
