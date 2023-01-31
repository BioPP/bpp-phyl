//
// File: AbstractWordSubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: 2009-01-08 00:00:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include <Bpp/Numeric/Matrix/EigenValue.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>

#include "AbstractWordSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/WordAlphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

#include <Bpp/App/ApplicationTools.h>

using namespace bpp;

// From the STL:
#include <cmath>
#include <complex>

using namespace std;

/******************************************************************************/

AbstractWordSubstitutionModel::AbstractWordSubstitutionModel(
  ModelList& modelList,
  const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractSubstitutionModel(
    modelList.getWordAlphabet(),
    std::shared_ptr<const StateMapInterface>(new CanonicalStateMap(modelList.getWordAlphabet(), false)),
    prefix),
  newAlphabet_(true),
  VSubMod_      (),
  VnestedPrefix_(),
  Vrate_        (modelList.size())
{
  size_t i, j;
  size_t n = modelList.size();

  // test whether two models are identical

  bool flag = 0;
  i = 0;
  j = 1;
  while (!flag && i < (n - 1))
  {
    if (modelList.getModel(i) == modelList.getModel(j))
      flag = 1;
    else
    {
      j++;
      if (j == n)
      {
        i++;
        j = i + 1;
      }
    }
  }

  if (!flag)
  {
    for (i = 0; i < n; i++)
    {
      VSubMod_.push_back(modelList.getModel(i));
      VnestedPrefix_.push_back(modelList.getModel(i)->getNamespace());
      VSubMod_[i]->setNamespace(prefix + TextTools::toString(i + 1) + "_" + VnestedPrefix_[i]);
      addParameters_(VSubMod_[i]->getParameters());
    }
  }
  else
  {
    string t = "";
    for (i = 0; i < n; i++)
    {
      VSubMod_.push_back(modelList.getModel(0));
      VnestedPrefix_.push_back(modelList.getModel(0)->getNamespace());
      t += TextTools::toString(i + 1);
    }
    VSubMod_[0]->setNamespace(prefix + t + "_" + VnestedPrefix_[0]);
    addParameters_(VSubMod_[0]->getParameters());
  }

  for (i = 0; i < n; i++)
  {
    Vrate_[i] = 1.0 / static_cast<double>(n);
  }
}

AbstractWordSubstitutionModel::AbstractWordSubstitutionModel(
  shared_ptr<const Alphabet> alph,
  shared_ptr<const StateMapInterface> stateMap,
  const string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractSubstitutionModel(alph, stateMap, prefix),
  newAlphabet_(false),
  VSubMod_      (),
  VnestedPrefix_(),
  Vrate_        (0)
{}

AbstractWordSubstitutionModel::AbstractWordSubstitutionModel(
  unique_ptr<SubstitutionModelInterface>& pmodel,
  unsigned int num,
  const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractSubstitutionModel(make_unique<WordAlphabet>(pmodel->getAlphabet(), num), 0, prefix),
  newAlphabet_(true),
  VSubMod_      (),
  VnestedPrefix_(),
  Vrate_        (num, 1.0 / num)
{
  stateMap_ = std::shared_ptr<const StateMapInterface>(new CanonicalStateMap(getAlphabet(), false));

  size_t i;

  string t = "";
  shared_ptr<SubstitutionModelInterface> pmodel2 = move(pmodel);
  for (i = 0; i < num; i++)
  {
    VSubMod_.push_back(pmodel2);
    VnestedPrefix_.push_back(pmodel2->getNamespace());
    t += TextTools::toString(i + 1);
  }

  pmodel->setNamespace(prefix + t + "_" + VnestedPrefix_[0]);
  addParameters_(pmodel->getParameters());
}

AbstractWordSubstitutionModel::AbstractWordSubstitutionModel(
  const AbstractWordSubstitutionModel& wrsm) :
  AbstractParameterAliasable(wrsm),
  AbstractSubstitutionModel(wrsm),
  newAlphabet_(wrsm.newAlphabet_),
  VSubMod_      (),
  VnestedPrefix_(wrsm.VnestedPrefix_),
  Vrate_        (wrsm.Vrate_)
{
  size_t i;
  size_t num = wrsm.VSubMod_.size();

  if (wrsm.newAlphabet_)
    alphabet_ = make_shared<WordAlphabet>(dynamic_cast<const WordAlphabet&>(wrsm.alphabet()));

  shared_ptr<SubstitutionModelInterface> pSM = nullptr;
  if ((num > 1) & (wrsm.VSubMod_[0] == wrsm.VSubMod_[1]))
    pSM.reset(wrsm.VSubMod_[0]->clone());

  for (i = 0; i < num; ++i)
  {
    VSubMod_.push_back(pSM ? pSM : shared_ptr<SubstitutionModelInterface>(wrsm.VSubMod_[i]->clone()));
  }
}

AbstractWordSubstitutionModel& AbstractWordSubstitutionModel::operator=(
  const AbstractWordSubstitutionModel& model)
{
  AbstractParameterAliasable::operator=(model);
  AbstractSubstitutionModel::operator=(model);
  newAlphabet_   = model.newAlphabet_;
  VnestedPrefix_ = model.VnestedPrefix_;
  Vrate_         = model.Vrate_;

  size_t i;
  size_t num = model.VSubMod_.size();

  if (model.newAlphabet_)
    alphabet_ = make_shared<WordAlphabet>(dynamic_cast<const WordAlphabet&>(model.alphabet()));

  shared_ptr<SubstitutionModelInterface> pSM = nullptr;
  if ((num > 1) & (model.VSubMod_[0] == model.VSubMod_[1]))
    pSM.reset(model.VSubMod_[0]->clone());

  for (i = 0; i < num; ++i)
  {
    VSubMod_[i] =  (pSM ? pSM : shared_ptr<SubstitutionModelInterface>(model.VSubMod_[i]->clone()));
  }

  return *this;
}

void AbstractWordSubstitutionModel::setNamespace(const std::string& prefix)
{
  AbstractSubstitutionModel::setNamespace(prefix);

  if (VSubMod_.size() < 2 || VSubMod_[0] == VSubMod_[1])
  {
    string t = "";
    for (size_t i = 0; i < VSubMod_.size(); i++)
    {
      t += TextTools::toString(i + 1);
    }
    VSubMod_[0]->setNamespace(prefix + t + "_" + VnestedPrefix_[0]);
  }
  else
  {
    for (size_t i = 0; i < VSubMod_.size(); i++)
    {
      VSubMod_[i]->setNamespace(prefix + TextTools::toString(i + 1) + "_" + VnestedPrefix_[i]);
    }
  }
}

/******************************************************************************/

void AbstractWordSubstitutionModel::updateMatrices_()
{
  // First we update position specific models. This need to be done
  // here and not in fireParameterChanged, as some parameter aliases
  // might have been defined and need to be resolved first.
  if (VSubMod_.size() < 2 || VSubMod_[0] == VSubMod_[1])
    VSubMod_[0]->matchParametersValues(getParameters());
  else
    for (size_t i = 0; i < VSubMod_.size(); i++)
    {
      VSubMod_[i]->matchParametersValues(getParameters());
    }

  size_t nbmod = VSubMod_.size();
  vector<bool> vnull; // vector of the indices of lines with only zeros

  // Generator

  size_t i, j, n, l, k, m;

  vector<size_t> vsize;

  for (k = 0; k < nbmod; k++)
  {
    vsize.push_back(VSubMod_[k]->getNumberOfStates());
  }

  RowMatrix<double> gk, exch;

  // First fill of the generator from simple position generators

  this->fillBasicGenerator_();

  // modification of generator_

  this->completeMatrices_();

  // sets diagonal terms

  setDiagonal();

  // at that point generator_ (and possibly freq_) are done for models
  // without enableEigenDecomposition

  // Eigen values:

  if (enableEigenDecomposition())
  {
    AbstractSubstitutionModel::updateMatrices_();
  }
  else  // compute freq_ if no eigenDecomposition
  {
    if (computeFrequencies())
    {
      size_t salph = getNumberOfStates();
      for (auto& fr : freq_)
      {
        fr = 1;
      }

      m = 1;
      for (k = nbmod; k > 0; k--)
      {
        auto& pSM = VSubMod_[k - 1];
        for (j = 0; j < vsize[k - 1]; ++j)
        {
          n = 0;
          while (n < salph)
          { // loop on prefix
            for (l = 0; l < m; ++l)
            { // loop on suffix
              freq_[n + j * m + l] *=  pSM->freq(j);
            }
            n += m * vsize[k - 1];
          }
        }
        m *= vsize[k - 1];
      }
      // normalization
      normalize();
    }
  }

  // compute the exchangeability_
  for (i = 0; i < size_; i++)
  {
    for (j = 0; j < size_; j++)
    {
      exchangeability_(i, j) = generator_(i, j) / freq_[j];
    }
  }
}


void AbstractWordSubstitutionModel::fillBasicGenerator_()
{
  size_t nbmod = VSubMod_.size();
  size_t salph = getNumberOfStates();

// Generator

  RowMatrix<double> gk;

  vector<size_t> vsize;

  for (size_t k = 0; k < nbmod; k++)
  {
    vsize.push_back(VSubMod_[k]->getNumberOfStates());
  }

  size_t m = 1;

  for (size_t k = nbmod; k > 0; k--)
  {
    gk = VSubMod_[k - 1]->getGenerator();
    for (size_t i = 0; i < vsize[k - 1]; i++)
    {
      const vector<double>& row_gi = gk.getRow(i);

      for (size_t j = 0; j < vsize[k - 1]; j++)
      {
        if (i != j)
        {
          size_t n = 0;
          while (n < salph)
          { // loop on prefix
            for (size_t l = 0; l < m; l++)
            { // loop on suffix
              generator_(n + i * m + l, n + j * m + l) = row_gi[j] * Vrate_[k - 1];
            }
            n += m * vsize[k - 1];
          }
        }
      }
    }
    m *= vsize[k - 1];
  }
}


void AbstractWordSubstitutionModel::setFreq(std::map<int, double>& freqs)
{
  map<int, double> tmpFreq;
  size_t nbmod = VSubMod_.size();

  size_t i, j, s, k, d, size;
  d = size = getNumberOfStates();

  if (VSubMod_.size() < 2 || VSubMod_[0] == VSubMod_[1])
  {
    s = VSubMod_[0]->getAlphabet()->getSize();
    for (j = 0; j < s; j++)
    {
      tmpFreq[static_cast<int>(j)] = 0;
    }

    for (i = 0; i < nbmod; i++)
    {
      d /= s;
      for (k = 0; k < size; k++)
      {
        tmpFreq[static_cast<int>((k / d) % s)] += freqs[static_cast<int>(k)];
      }
    }

    for (k = 0; k < s; k++)
    {
      tmpFreq[static_cast<int>(k)] /= static_cast<double>(nbmod);
    }

    VSubMod_[0]->setFreq(tmpFreq);
    matchParametersValues(VSubMod_[0]->getParameters());
  }
  else
    for (i = 0; i < nbmod; i++)
    {
      tmpFreq.clear();
      s = VSubMod_[i]->getAlphabet()->getSize();
      d /= s;
      for (j = 0; j < s; j++)
      {
        tmpFreq[static_cast<int>(j)] = 0;
      }
      for (k = 0; k < size; k++)
      {
        tmpFreq[static_cast<int>((k / d) % s)] += freqs[static_cast<int>(k)];
      }
      VSubMod_[i]->setFreq(tmpFreq);
      matchParametersValues(VSubMod_[i]->getParameters());
    }
}
