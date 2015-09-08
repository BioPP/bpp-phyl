//
// File: WordSubstitutionModel.cpp
// Created by:  Laurent Gueguen
// Created on: Jan 2009
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)
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

#include "WordSubstitutionModel.h"
#include "FrequenciesSet/WordFrequenciesSet.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/EigenValue.h>

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
    addParameter_(new Parameter("Word.relrate" + TextTools::toString(i + 1), 1.0 / static_cast<int>(nbmod - i), &Parameter::PROP_CONSTRAINT_EX));
  }

  WordSubstitutionModel::updateMatrices();
}

WordSubstitutionModel::WordSubstitutionModel(
  const Alphabet* alph,
  StateMap* stateMap,
  const std::string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Word." : prefix),
  AbstractWordSubstitutionModel(alph, stateMap, (prefix == "") ? "Word." : prefix)
{}

WordSubstitutionModel::WordSubstitutionModel(
  SubstitutionModel* pmodel,
  unsigned int num,
  const std::string& prefix) :
  AbstractParameterAliasable((prefix == "") ? "Word." : prefix),
  AbstractWordSubstitutionModel(pmodel,
                                num,
                                (prefix == "") ? "Word." : prefix)
{
  size_t i;

  // relative rates
  for (i = 0; i < num - 1; i++)
  {
    addParameter_(new Parameter("Word.relrate" + TextTools::toString(i + 1), 1.0 / static_cast<int>(num - i ), &Parameter::PROP_CONSTRAINT_EX));
  }

  WordSubstitutionModel::updateMatrices();
}

void WordSubstitutionModel::updateMatrices()
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

  AbstractWordSubstitutionModel::updateMatrices();
}

void WordSubstitutionModel::completeMatrices()
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
      m = VSubMod_[p-1]->getNumberOfStates();
      freq_[i] *= VSubMod_[p-1]->getFrequencies()[j % m];
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


