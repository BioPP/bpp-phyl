// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "YpR.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/EigenValue.h>

#include <Bpp/Seq/Alphabet/AlphabetTools.h>

#include <Bpp/Text/TextTools.h>

/******************************************************************************/

YpR::YpR(
    std::shared_ptr<const RNY> alph,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pm,
    const string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractSubstitutionModel(alph, make_shared<CanonicalStateMap>(alph, false), prefix),
  pmodel_(pm->clone()),
  nestedPrefix_(pm->getNamespace())
{
  pmodel_->setNamespace(prefix + nestedPrefix_);
  pmodel_->enableEigenDecomposition(0);
  pmodel_->computeFrequencies(false);

  addParameters_(pmodel_->getParameters());
  computeFrequencies(true);
}

YpR::YpR(const YpR& ypr) :
  AbstractParameterAliasable(ypr),
  AbstractSubstitutionModel(ypr),
  pmodel_(ypr.pmodel_->clone()),
  nestedPrefix_(ypr.getNestedPrefix())
{
}

void YpR::updateMatrices_()
{
  updateMatrices_(0, 0, 0, 0, 0, 0, 0, 0);
}


void YpR::updateMatrices_(double CgT, double cGA,
    double TgC, double tGA,
    double CaT, double cAG,
    double TaC, double tAC)
{
  checkModel(*pmodel_);

  // Generator:
  unsigned int i, j, i1, i2, i3, j1, j2, j3;

  std::vector<double> a(4);  // a[A], a[G], a[C], a[T]
  std::vector<double> b(4);  // b[A], b[G], b[C], b[T]

  // From M

  a[0] = pmodel_->Qij(2, 0);
  a[1] = pmodel_->Qij(0, 2);
  a[2] = pmodel_->Qij(3, 1);
  a[3] = pmodel_->Qij(1, 3);
  b[0] = (pmodel_->Qij(1, 0) + pmodel_->Qij(3, 0))/2;   // To limit numerical issues
  b[1] = (pmodel_->Qij(1, 2) + pmodel_->Qij(3, 2))/2;
  b[2] = (pmodel_->Qij(0, 1) + pmodel_->Qij(2, 1))/2;
  b[3] = (pmodel_->Qij(0, 3) + pmodel_->Qij(2, 3))/2;

  

  // M_1 on R C T
  RowMatrix<double> M1(3, 3);

  M1(0, 0) = 0;
  M1(0, 1) = b[2];
  M1(0, 2) = b[3];
  M1(1, 0) = pmodel_->Qij(1, 0) + pmodel_->Qij(1, 2);
  M1(1, 1) = 0;
  M1(1, 2) = pmodel_->Qij(1, 3);
  M1(2, 0) = pmodel_->Qij(3, 0) + pmodel_->Qij(3, 2);
  M1(2, 1) = pmodel_->Qij(3, 1);
  M1(2, 2) = 0;

  // M_2 on A G C T
  RowMatrix<double> M2(4, 4);

  M2(0, 0) = 0;
  M2(0, 1) = pmodel_->Qij(0, 2);
  M2(0, 2) = pmodel_->Qij(0, 1);
  M2(0, 3) = pmodel_->Qij(0, 3);
  M2(1, 0) = pmodel_->Qij(2, 0);
  M2(1, 1) = 0;
  M2(1, 2) = pmodel_->Qij(2, 1);
  M2(1, 3) = pmodel_->Qij(2, 3);
  M2(2, 0) = pmodel_->Qij(1, 0);
  M2(2, 1) = pmodel_->Qij(1, 2);
  M2(2, 2) = 0;
  M2(2, 3) = pmodel_->Qij(1, 3);
  M2(3, 0) = pmodel_->Qij(3, 0);
  M2(3, 1) = pmodel_->Qij(3, 2);
  M2(3, 2) = pmodel_->Qij(3, 1);
  M2(3, 3) = 0;   

  // M_3 on A G Y
  RowMatrix<double> M3(3, 3);

  M3(0, 0) = 0;
  M3(0, 1) = pmodel_->Qij(0, 2);
  M3(0, 2) = pmodel_->Qij(0, 1) + pmodel_->Qij(0, 3);
  M3(1, 0) = pmodel_->Qij(2, 0);
  M3(1, 1) = 0;
  M3(1, 2) = pmodel_->Qij(2, 1) + pmodel_->Qij(2, 3);
  M3(2, 0) = b[0];
  M3(2, 1) = b[1];
  M3(2, 2) = 0;


  for (i1 = 0; i1 < 3; i1++)
  {
    for (i2 = 0; i2 < 4; i2++)
    {
      for (i3 = 0; i3 < 3; i3++)
      {
        i = 12 * i1 + 3 * i2 + i3;
        for (j1 = 0; j1 < 3; j1++)
        {
          for (j2 = 0; j2 < 4; j2++)
          {
            for (j3 = 0; j3 < 3; j3++)
            {
              j = 12 * j1 + 3 * j2 + j3;
              if ((i1 == j1) && (i2 == j2))
                generator_(i, j) = M3(i3, j3);
              else if ((i1 == j1) && (i3 == j3))
                generator_(i, j) = M2(i2, j2);
              else if ((i2 == j2) && (i3 == j3))
                generator_(i, j) = M1(i1, j1);
              else
                generator_(i, j) = 0;
            }
          }
        }
      }
    }
  }

  // Introduction des dependances

  for (i3 = 0; i3 < 3; i3++)
  {
    generator_(15 + i3, 12 + i3) += cGA * a[0]; // CG -> CA
    generator_(12 * i3 + 7, 12 * i3 + 6) += cGA * a[0];

    generator_(15 + i3, 27 + i3) += CgT * a[3]; // CG -> TG
    generator_(12 * i3 + 7, 12 * i3 + 10) += CgT * a[3];

    generator_(27 + i3, 24 + i3) += tGA * a[0]; // TG -> TA
    generator_(12 * i3 + 10, 12 * i3 + 9) += tGA * a[0];

    generator_(27 + i3, 15 + i3) += TgC * a[2]; // TG -> CG
    generator_(12 * i3 + 10, 12 * i3 + 7) += TgC * a[2];

    generator_(12 + i3, 24 + i3) += CaT * a[3]; // CA -> TA
    generator_(12 * i3 + 6, 12 * i3 + 9) += CaT * a[3];

    generator_(12 + i3, 15 + i3) += cAG * a[1]; // CA -> CG
    generator_(12 * i3 + 6, 12 * i3 + 7) += cAG * a[1];

    generator_(24 + i3, 27 + i3) += tAC * a[1]; // TA -> TG
    generator_(12 * i3 + 9, 12 * i3 + 10) += tAC * a[1];

    generator_(24 + i3, 12 + i3) += TaC * a[2]; // TA -> CA
    generator_(12 * i3 + 9, 12 * i3 + 6) += TaC * a[2];
  }

  double x;

  for (i = 0; i < 36; ++i)
  {
    x = 0;
    for (j = 0; j < 36; ++j)
    {
      if (j != i)
        x += generator_(i, j);
    }
    generator_(i, i) = -x;
  }

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
    for (i = 0; i < size_; i++)
    {
      if (abs(iEigenValues_[i]) > NumConstants::TINY())
      {
        isDiagonalizable_ = false;
      }
    }

    // frequence stationnaire

    x = 0;
    j = 0;
    while (j < 36)
    {
      if (abs(eigenValues_[j]) < NumConstants::SMALL() &&
          abs(iEigenValues_[j]) < NumConstants::SMALL())
      {
        eigenValues_[j] = 0; // to avoid approximation problems in the future
        for (i = 0; i < 36; i++)
        {
          freq_[i] = leftEigenVectors_(j, i);
          x += freq_[i];
        }
        break;
      }
      j++;
    }
    for (i = 0; i < 36; i++)
    {
      freq_[i] /= x;
    }
  }
  catch (ZeroDivisionException& e)
  {
    ApplicationTools::displayMessage("Singularity during  diagonalization. Taylor series used instead.");
    isNonSingular_ = false;
    isDiagonalizable_ = false;

    if (vPowGen_.size() == 0)
      vPowGen_.resize(30);

    double min = generator_(0, 0);
    for (i = 1; i < 36; i++)
    {
      if (min > generator_(i, i))
        min = generator_(i, i);
    }

    setScale(-1 / min);

    MatrixTools::getId(36, tmpMat_);    // to compute the equilibrium frequency  (Q+Id)^256

    MatrixTools::add(tmpMat_, generator_);
    MatrixTools::pow(tmpMat_, 256, vPowGen_[0]);

    for (i = 0; i < 36; i++)
    {
      freq_[i] = vPowGen_[0](0, i);
    }

    MatrixTools::getId(36, vPowGen_[0]);
  }

  // mise a l'echelle

  x = 0;
  for (i1 = 0; i1 < 3; i1++)
  {
    for (i2 = 0; i2 < 4; i2++)
    {
      for (i3 = 0; i3 < 3; i3++)
      {
        i = 12 * i1 + 3 * i2 + i3;
        for (j2 = 0; j2 < 4; j2++)
        {
          if (j2 != i2)
          {
            j = 12 * i1 + 3 * j2 + i3;
            x += freq_[i] * generator_(i, j);
          }
        }
      }
    }
  }

  setScale(1 / x);

  if (isScalable())
  {
    for (i = 0; i < 36; i++)
    {
      eigenValues_[i] /= x;
      iEigenValues_[i] /= x;
    }
  }

  if (!isNonSingular_)
    MatrixTools::Taylor(generator_, 30, vPowGen_);

  // and the exchangeability_
  for (i = 0; i < size_; i++)
  {
    for (j = 0; j < size_; j++)
    {
      exchangeability_(i, j) = generator_(i, j) / freq_[j];
    }
  }
}

void YpR::checkModel(const SubstitutionModelInterface& pm) const
{
  auto alph = pm.getAlphabet();
  if (!AlphabetTools::isNucleicAlphabet(*alph))
    throw Exception("Need a DNA model");

  // Check that the model is good for YpR, ie transversion rates do
  // not depend on the origin state

  if ((pm.Qij(0, 1) != pm.Qij(2, 1)) || (pm.Qij(0, 3) != pm.Qij(2, 3))
      || (pm.Qij(1, 0) != pm.Qij(3, 0)) || (pm.Qij(1, 2) != pm.Qij(3, 2)))
      throw Exception("Not R/Y Model " + pm.getName());
}

void YpR::setNamespace(const std::string& prefix)
{
  AbstractSubstitutionModel::setNamespace(prefix);
  
  // We also need to update the namespace of the nested model:
  pmodel_->setNamespace(prefix + nestedPrefix_);
}

// ///////////////////////////////////////////////
// ///////////////////////////////////////////////


/******************************************************************************/

YpR_Sym::YpR_Sym(
    std::shared_ptr<const RNY> alph,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pm,
    double CgT, double TgC,
    double CaT, double TaC) :
  AbstractParameterAliasable("YpR_Sym."),
  YpR(alph, std::move(pm), "YpR_Sym.")
{
  addParameter_(new Parameter("YpR_Sym.rCgT", CgT, Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Sym.rTgC", TgC, Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Sym.rCaT", CaT, Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Sym.rTaC", TaC, Parameter::R_PLUS));

  updateMatrices_();
}

void YpR_Sym::updateMatrices_()
{
  double rCgT = getParameterValue("rCgT");
  double rTgC = getParameterValue("rTgC");
  double rCaT = getParameterValue("rCaT");
  double rTaC = getParameterValue("rTaC");

  YpR::updateMatrices_(rCgT, rCgT, rTgC, rTgC, rCaT, rCaT, rTaC, rTaC);
}

YpR_Sym::YpR_Sym(const YpR_Sym& ypr) :
  AbstractParameterAliasable(ypr),
  YpR(ypr)
{}

/******************************************************************************/

std::string YpR_Sym::getName() const
{
  return "YpR_Sym";
}


// ///////////////////////////////////////////////
// ///////////////////////////////////////////////

/******************************************************************************/

YpR_Gen::YpR_Gen(
    std::shared_ptr<const RNY> alph,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pm,
    double CgT, double cGA,
    double TgC, double tGA,
    double CaT, double cAG,
    double TaC, double tAG) :
  AbstractParameterAliasable("YpR_Gen."),
  YpR(alph, std::move(pm), "YpR_Gen.")
{
  addParameter_(new Parameter("YpR_Gen.rCgT", CgT, Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rcGA", cGA, Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rTgC", TgC, Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rtGA", tGA, Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rCaT", CaT, Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rcAG", cAG, Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rTaC", TaC, Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rtAG", tAG, Parameter::R_PLUS));

  updateMatrices_();
}

void YpR_Gen::updateMatrices_()
{
  double rCgT = getParameterValue("rCgT");
  double rcGA = getParameterValue("rcGA");
  double rTgC = getParameterValue("rTgC");
  double rtGA = getParameterValue("rtGA");
  double rCaT = getParameterValue("rCaT");
  double rcAG = getParameterValue("rcAG");
  double rTaC = getParameterValue("rTaC");
  double rtAG = getParameterValue("rtAG");

  YpR::updateMatrices_(rCgT, rcGA, rTgC, rtGA, rCaT, rcAG, rTaC, rtAG);
}

YpR_Gen::YpR_Gen(const YpR_Gen& ypr) :
  AbstractParameterAliasable(ypr),
  YpR(ypr)
{
  updateMatrices_();
}

/******************************************************************************/

std::string YpR_Gen::getName() const
{
  return "YpR_Gen";
}
