// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/EigenValue.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>

#include "WordSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/WordAlphabet.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/******************************************************************************/

WordSubstitutionModel::WordSubstitutionModel(
    ModelList& modelList,
    const std::string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Word." : prefix),
  AbstractWordSubstitutionModel(
      modelList,
      (prefix == "") ? "Word." : prefix)
{
  size_t i, nbmod = VSubMod_.size();

  // relative rates
  for (i = 0; i < nbmod - 1; i++)
  {
    addParameter_(new Parameter("Word.relrate" + TextTools::toString(i + 1), 1.0 / static_cast<int>(nbmod - i), Parameter::PROP_CONSTRAINT_EX));
  }

  enableEigenDecomposition(false); // the product of the position
                                   // specific transition probabilities

  computeFrequencies(false); // it is done in AbstractWordSubstitutionModel
  WordSubstitutionModel::updateMatrices_();
}

WordSubstitutionModel::WordSubstitutionModel(
    shared_ptr<const Alphabet> alph,
    shared_ptr<const StateMapInterface> stateMap,
    const string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Word." : prefix),
  AbstractWordSubstitutionModel(alph, stateMap, (prefix == "") ? "Word." : prefix)
{
  enableEigenDecomposition(false); // the product of the position
                                   // specific transition probabilities
  computeFrequencies(false); // it is done in AbstractWordSubstitutionModel
}

WordSubstitutionModel::WordSubstitutionModel(
    unique_ptr<SubstitutionModelInterface> pmodel,
    unsigned int num,
    const std::string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Word." : prefix),
  AbstractWordSubstitutionModel(std::move(pmodel),
      num,
      (prefix == "") ? "Word." : prefix)
{
  size_t i;

  // relative rates
  for (i = 0; i < num - 1; i++)
  {
    addParameter_(new Parameter("Word.relrate" + TextTools::toString(i + 1), 1.0 / static_cast<int>(num - i ), Parameter::PROP_CONSTRAINT_EX));
  }

  enableEigenDecomposition(false); // the product of the position
                                   // specific transition probabilities
  computeFrequencies(false); // it is done in AbstractWordSubstitutionModel
  WordSubstitutionModel::updateMatrices_();
}

void WordSubstitutionModel::updateMatrices_()
{
  size_t i, nbmod = VSubMod_.size();
  double x, y;
  x = 1.0;

  for (i = 0; i < nbmod - 1; i++)
  {
    y = getParameterValue("relrate" + TextTools::toString(i + 1));
    Vrate_[i] = x * y;
    x *= 1 - y;
  }
  Vrate_[nbmod - 1] = x;

  AbstractWordSubstitutionModel::updateMatrices_();
}

void WordSubstitutionModel::completeMatrices_()
{
  size_t nbmod = VSubMod_.size();
  size_t i, p, j, m;
  size_t salph = getAlphabet()->getSize();

  // freq_ for this generator

  for (i = 0; i < salph; i++)
  {
    freq_[i] = 1;
    j = i;
    for (p = nbmod; p > 0; p--)
    {
      m = VSubMod_[p - 1]->getNumberOfStates();
      freq_[i] *= VSubMod_[p - 1]->getFrequencies()[j % m];
      j /= m;
    }
  }
}

const RowMatrix<double>& WordSubstitutionModel::getPij_t(double d) const
{
  vector<const Matrix<double>*> vM;
  size_t nbmod = VSubMod_.size();
  size_t i, j;

  for (i = 0; i < nbmod; i++)
  {
    vM.push_back(&VSubMod_[i]->getPij_t(d * Vrate_[i] * rate_));
  }

  size_t t;
  double x;
  size_t i2, j2;
  size_t nbStates = getNumberOfStates();
  size_t p;

  for (i = 0; i < nbStates; i++)
  {
    for (j = 0; j < nbStates; j++)
    {
      x = 1.;
      i2 = i;
      j2 = j;
      for (p = nbmod; p > 0; p--)
      {
        t = VSubMod_[p - 1]->getNumberOfStates();
        x *= (*vM[p - 1])(i2 % t, j2 % t);
        i2 /= t;
        j2 /= t;
      }
      pijt_(i, j) = x;
    }
  }

  return pijt_;
}

const RowMatrix<double>& WordSubstitutionModel::getdPij_dt(double d) const
{
  vector<const Matrix<double>*> vM, vdM;
  size_t nbmod = VSubMod_.size();
  size_t i, j;

  for (i = 0; i < nbmod; i++)
  {
    vM.push_back(&VSubMod_[i]->getPij_t(d * Vrate_[i] * rate_));
    vdM.push_back(&VSubMod_[i]->getdPij_dt(d * Vrate_[i] * rate_));
  }

  size_t t;
  double x, r;
  size_t i2, j2;
  size_t nbStates = getNumberOfStates();
  size_t p, q;

  for (i = 0; i < nbStates; i++)
  {
    for (j = 0; j < nbStates; j++)
    {
      r = 0;
      for (q = 0; q < nbmod; q++)
      {
        i2 = i;
        j2 = j;
        x = 1;
        for (p = nbmod; p > 0; p--)
        {
          t = VSubMod_[p - 1]->getNumberOfStates();
          if (q != p - 1)
            x *= (*vM[p - 1])(i2 % t, j2 % t);
          else
            x *= rate_ * Vrate_[p - 1] * (*vdM[p - 1])(i2 % t, j2 % t);
          i2 /= t;
          j2 /= t;
        }
        r += x;
      }
      dpijt_(i, j) = r;
    }
  }
  return dpijt_;
}

const RowMatrix<double>& WordSubstitutionModel::getd2Pij_dt2(double d) const

{
  vector<const Matrix<double>*> vM, vdM, vd2M;
  size_t nbmod = VSubMod_.size();
  size_t i, j;

  for (i = 0; i < nbmod; i++)
  {
    vM.push_back(&VSubMod_[i]->getPij_t(d * Vrate_[i] * rate_));
    vdM.push_back(&VSubMod_[i]->getdPij_dt(d * Vrate_[i] * rate_));
    vd2M.push_back(&VSubMod_[i]->getd2Pij_dt2(d * Vrate_[i] * rate_));
  }


  double r, x;
  size_t p, b, q, t;

  size_t i2, j2;
  size_t nbStates = getNumberOfStates();


  for (i = 0; i < nbStates; i++)
  {
    for (j = 0; j < nbStates; j++)
    {
      r = 0;
      for (q = 1; q < nbmod; q++)
      {
        for (b = 0; b < q; b++)
        {
          x = 1;
          i2 = i;
          j2 = j;
          for (p = nbmod; p > 0; p--)
          {
            t = VSubMod_[p - 1]->getNumberOfStates();
            if ((p - 1 == q) || (p - 1 == b))
              x *= rate_ * Vrate_[p - 1] * (*vdM[p - 1])(i2 % t, j2 % t);
            else
              x *= (*vM[p - 1])(i2 % t, j2 % t);

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
        for (p = nbmod; p > 0; p--)
        {
          t = VSubMod_[p - 1]->getNumberOfStates();
          if (q != p - 1)
            x *= (*vM[p - 1])(i2 % t, j2 % t);
          else
            x *= rate_ * rate_ * Vrate_[p - 1] * Vrate_[p - 1] * (*vd2M[p - 1])(i2 % t, j2 % t);

          i2 /= t;
          j2 /= t;
        }
        r += x;
      }
      d2pijt_(i, j) = r;
    }
  }
  return d2pijt_;
}

string WordSubstitutionModel::getName() const
{
  return "Word";
}
