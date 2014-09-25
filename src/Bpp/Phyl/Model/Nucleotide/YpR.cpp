//
// File: YpR.cpp
// Created by: Laurent Gueguen
// Created on: Thu August 2 2007
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

#include "YpR.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/EigenValue.h>

#include <Bpp/Text/TextTools.h>

/******************************************************************************/

YpR::YpR(const RNY* alph, SubstitutionModel* const pm, const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractSubstitutionModel(alph, new CanonicalStateMap(alph, false), prefix),
  pmodel_(pm->clone()),
  _nestedPrefix(pm->getNamespace())
{
  pmodel_->setNamespace(prefix + _nestedPrefix);
  pmodel_->enableEigenDecomposition(0);
  addParameters_(pmodel_->getParameters());
}

YpR::YpR(const YpR& ypr, const std::string& prefix) :
  AbstractParameterAliasable(ypr),
  AbstractSubstitutionModel(ypr),
  pmodel_(ypr.pmodel_->clone()),
  _nestedPrefix(ypr.getNestedPrefix())

{
  pmodel_->setNamespace(prefix + _nestedPrefix);
}

YpR::YpR(const YpR& ypr) :
  AbstractParameterAliasable(ypr),
  AbstractSubstitutionModel(ypr),
  pmodel_(ypr.pmodel_->clone()),
  _nestedPrefix(ypr.getNestedPrefix())
{}

void YpR::updateMatrices()
{
  updateMatrices(0, 0, 0, 0, 0, 0, 0, 0);
}


void YpR::updateMatrices(double CgT, double cGA,
                         double TgC, double tGA,
                         double CaT, double cAG,
                         double TaC, double tAC)
{
  //  check_model(pmodel_);

  // Generator:
  const Alphabet* alph = pmodel_->getAlphabet();
  std::vector<size_t> l(4);

  l[0] = alph->getStateIndex("A");
  l[1] = alph->getStateIndex("G");
  l[2] = alph->getStateIndex("C");
  l[3] = alph->getStateIndex("T");

  unsigned int i, j, i1, i2, i3, j1, j2, j3;

  std::vector<double> a(4);  // a[A], a[G], a[C], a[T]
  std::vector<double> b(4);  // b[A], b[G], b[C], b[T]

  for (i = 0; i < 2; i++)
  {
    a[i] = pmodel_->Qij(l[1 - i], l[i]);
    b[i] = pmodel_->Qij(l[3 - i], l[i]);
    a[i + 2] = pmodel_->Qij(l[3 - i], l[i + 2]);
    b[i + 2] = pmodel_->Qij(l[1 - i], l[i + 2]);
  }

  // M_1
  RowMatrix<double> M1(3, 3);

  M1(0, 0) = 0;
  M1(0, 1) = b[2];
  M1(0, 2) = b[3];
  M1(1, 0) = b[0] + b[1];
  M1(1, 1) = 0;
  M1(1, 2) = a[3];
  M1(2, 0) = b[0] + b[1];
  M1(2, 1) = a[2];
  M1(2, 2) = 0;

  // M_2
  RowMatrix<double> M2(4, 4);

  M2(0, 0) = 0;
  M2(0, 1) = a[1];
  M2(0, 2) = b[2];
  M2(0, 3) = b[3];
  M2(1, 0) = a[0];
  M2(1, 1) = 0;
  M2(1, 2) = b[2];
  M2(1, 3) = b[3];
  M2(2, 0) = b[0];
  M2(2, 1) = b[1];
  M2(2, 2) = 0;
  M2(2, 3) = a[3];
  M2(3, 0) = b[0];
  M2(3, 1) = b[1];
  M2(3, 2) = a[2];
  M2(3, 3) = 0;

  // M_3
  RowMatrix<double> M3(3, 3);

  M3(0, 0) = 0;
  M3(0, 1) = a[1];
  M3(0, 2) = b[2] + b[3];
  M3(1, 0) = a[0];
  M3(1, 1) = 0;
  M3(1, 2) = b[2] + b[3];
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

  for (i = 0; i < 36; i++)
  {
    x = 0;
    for (j = 0; j < 36; j++)
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

    MatrixTools::scale(generator_, -1 / min);

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

  MatrixTools::scale(generator_, 1 / x);

  if (!isNonSingular_)
    MatrixTools::Taylor(generator_, 30, vPowGen_);

  for (i = 0; i < 36; i++)
  {
    eigenValues_[i] /= x;
    iEigenValues_[i] /= x;
  }

  // and the exchangeability_
  for (i = 0; i < size_; i++)
  {
    for (j = 0; j < size_; j++)
    {
      exchangeability_(i, j) = generator_(i, j) / freq_[j];
    }
  }
}

void YpR::check_model(SubstitutionModel* const pm) const
throw (Exception)
{
  if (!pm)
    throw Exception("No Model ");

  const Alphabet* alph = pm->getAlphabet();
  if (alph->getAlphabetType() != "DNA alphabet")
    throw Exception("Need a DNA model");

  std::vector<size_t> l(4);

  l[0] = alph->getStateIndex("A");
  l[1] = alph->getStateIndex("G");
  l[2] = alph->getStateIndex("C");
  l[3] = alph->getStateIndex("T");

  // Check that the model is good for YpR

  for (size_t i = 0; i < 2; ++i)
  {
    if (pm->Qij(l[2], l[i]) != pm->Qij(l[3], l[i]))
      throw Exception("Not R/Y Model " + pm->getName());
  }
  for (size_t i = 2; i < 4; ++i)
  {
    if (pm->Qij(l[0], l[i]) != pm->Qij(l[1], l[i]))
      throw Exception("Not R/Y Model " + pm->getName());
  }
}

void YpR::setNamespace(const std::string& prefix)
{
  AbstractSubstitutionModel::setNamespace(prefix);
  // We also need to update the namespace of the nested model:
  pmodel_->setNamespace(prefix + _nestedPrefix);
}

// ///////////////////////////////////////////////
// ///////////////////////////////////////////////


/******************************************************************************/

YpR_Sym::YpR_Sym(const RNY* alph,
                 SubstitutionModel* const pm,
                 double CgT, double TgC,
                 double CaT, double TaC) : AbstractParameterAliasable("YpR_Sym."),
  YpR(alph, pm, "YpR_Sym.")
{
  addParameter_(new Parameter("YpR_Sym.rCgT", CgT, &Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Sym.rTgC", TgC, &Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Sym.rCaT", CaT, &Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Sym.rTaC", TaC, &Parameter::R_PLUS));

  updateMatrices();
}

void YpR_Sym::updateMatrices()
{
  double rCgT = getParameterValue("rCgT");
  double rTgC = getParameterValue("rTgC");
  double rCaT = getParameterValue("rCaT");
  double rTaC = getParameterValue("rTaC");

  YpR::updateMatrices(rCgT, rCgT, rTgC, rTgC, rCaT, rCaT, rTaC, rTaC);
}

YpR_Sym::YpR_Sym(const YpR_Sym& ypr) : AbstractParameterAliasable(ypr),
  YpR(ypr, "YpR_Sym.")
{}

/******************************************************************************/

std::string YpR_Sym::getName() const
{
  return "YpR_Sym";
}


// ///////////////////////////////////////////////
// ///////////////////////////////////////////////

/******************************************************************************/

YpR_Gen::YpR_Gen(const RNY* alph,
                 SubstitutionModel* const pm,
                 double CgT, double cGA,
                 double TgC, double tGA,
                 double CaT, double cAG,
                 double TaC, double tAG) : AbstractParameterAliasable("YpR_Gen."),
  YpR(alph, pm, "YpR_Gen.")
{
  addParameter_(new Parameter("YpR_Gen.rCgT", CgT, &Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rcGA", cGA, &Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rTgC", TgC, &Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rtGA", tGA, &Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rCaT", CaT, &Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rcAG", cAG, &Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rTaC", TaC, &Parameter::R_PLUS));
  addParameter_(new Parameter("YpR_Gen.rtAG", tAG, &Parameter::R_PLUS));

  updateMatrices();
}

void YpR_Gen::updateMatrices()
{
  double rCgT = getParameterValue("rCgT");
  double rcGA = getParameterValue("rcGA");
  double rTgC = getParameterValue("rTgC");
  double rtGA = getParameterValue("rtGA");
  double rCaT = getParameterValue("rCaT");
  double rcAG = getParameterValue("rcAG");
  double rTaC = getParameterValue("rTaC");
  double rtAG = getParameterValue("rtAG");

  YpR::updateMatrices(rCgT, rcGA, rTgC, rtGA, rCaT, rcAG, rTaC, rtAG);
}

YpR_Gen::YpR_Gen(const YpR_Gen& ypr) : AbstractParameterAliasable(ypr),
  YpR(ypr, "YpR_Gen.")
{
  updateMatrices();
}

/******************************************************************************/

std::string YpR_Gen::getName() const
{
  return "YpR_Gen";
}
