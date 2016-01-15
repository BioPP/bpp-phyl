//
// File: AbstractWordSubstitutionModel.cpp
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

#include "AbstractWordSubstitutionModel.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Matrix/EigenValue.h>
#include <Bpp/Numeric/VectorTools.h>

// From SeqLib:
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
      new CanonicalStateMap(modelList.getWordAlphabet(), false),
      prefix),
  new_alphabet_ (true),
  VSubMod_      (),
  VnestedPrefix_(),
  Vrate_        (modelList.size())
{
  enableEigenDecomposition(false);
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
  const Alphabet* alph,
  StateMap* stateMap,
  const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractSubstitutionModel(alph, stateMap, prefix),
  new_alphabet_ (false),
  VSubMod_      (),
  VnestedPrefix_(),
  Vrate_         (0)
{
  enableEigenDecomposition(false);
}

AbstractWordSubstitutionModel::AbstractWordSubstitutionModel(
  SubstitutionModel* pmodel,
  unsigned int num,
  const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractSubstitutionModel(new WordAlphabet(pmodel->getAlphabet(), num), pmodel->getStateMap().clone(), prefix),
  new_alphabet_ (true),
  VSubMod_      (),
  VnestedPrefix_(),
  Vrate_         (num)
{
  enableEigenDecomposition(false);
  size_t i;

  string t = "";
  for (i = 0; i < num; i++)
  {
    VSubMod_.push_back(pmodel);
    VnestedPrefix_.push_back(pmodel->getNamespace());
    Vrate_[i] = 1.0 / num;
    t += TextTools::toString(i + 1);
  }

  pmodel->setNamespace(prefix + t + "_" + VnestedPrefix_[0]);
  addParameters_(pmodel->getParameters());
}

AbstractWordSubstitutionModel::AbstractWordSubstitutionModel(
  const AbstractWordSubstitutionModel& wrsm) :
  AbstractParameterAliasable(wrsm),
  AbstractSubstitutionModel(wrsm),
  new_alphabet_ (wrsm.new_alphabet_),
  VSubMod_      (),
  VnestedPrefix_(wrsm.VnestedPrefix_),
  Vrate_         (wrsm.Vrate_)
{
  size_t i;
  size_t num = wrsm.VSubMod_.size();

  if (wrsm.new_alphabet_)
    alphabet_ = new WordAlphabet(*(dynamic_cast<const WordAlphabet*>(wrsm.getAlphabet())));

  SubstitutionModel* pSM = 0;
  if ((num > 1) & (wrsm.VSubMod_[0] == wrsm.VSubMod_[1]))
    pSM = wrsm.VSubMod_[0]->clone();

  for (i = 0; i < num; i++)
  {
    VSubMod_.push_back(pSM ? pSM : wrsm.VSubMod_[i]->clone());
  }
}

AbstractWordSubstitutionModel& AbstractWordSubstitutionModel::operator=(
  const AbstractWordSubstitutionModel& model)
{
  AbstractParameterAliasable::operator=(model);
  AbstractSubstitutionModel::operator=(model);
  new_alphabet_  = model.new_alphabet_;
  VnestedPrefix_ = model.VnestedPrefix_;
  Vrate_         = model.Vrate_;

  size_t i;
  size_t num = model.VSubMod_.size();

  if (model.new_alphabet_)
    alphabet_ = new WordAlphabet(*(dynamic_cast<const WordAlphabet*>(model.getAlphabet())));

  SubstitutionModel* pSM = 0;
  if ((num > 1) & (model.VSubMod_[0] == model.VSubMod_[1]))
    pSM = model.VSubMod_[0]->clone();

  for (i = 0; i < num; i++)
  {
    VSubMod_[i] =  (pSM ? pSM : model.VSubMod_[i]->clone());
  }

  return *this;
}

AbstractWordSubstitutionModel::~AbstractWordSubstitutionModel()
{
  if ((VSubMod_.size() > 1) && (VSubMod_[0] == VSubMod_[1]))
  {
    if (VSubMod_[0])
      delete VSubMod_[0];
  }
  else
    for (size_t i = 0; i < VSubMod_.size(); i++)
    {
      if (VSubMod_[i])
        delete VSubMod_[i];
    }
  if (new_alphabet_)
    delete alphabet_;
}

size_t AbstractWordSubstitutionModel::getNumberOfStates() const
{
  return getAlphabet()->getSize();
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

void AbstractWordSubstitutionModel::updateMatrices()
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
  size_t salph = getNumberOfStates();
  size_t nbStop = 0;
  vector<bool> vnull; // vector of the indices of lines with only zeros

  // Generator

  size_t i, j, n, l, k, m;

  vector<size_t> vsize;

  for (k = 0; k < nbmod; k++)
  {
    vsize.push_back(VSubMod_[k]->getNumberOfStates());
  }

  RowMatrix<double> gk, exch;

  m = 1;
  
  for (k = nbmod; k > 0; k--)
  {
    gk = VSubMod_[k - 1]->getGenerator();
    for (i = 0; i < vsize[k - 1]; i++)
    {
      for (j = 0; j < vsize[k - 1]; j++)
      {
        if (i != j)
        {
          n = 0;
          while (n < salph)
          { // loop on prefix
            for (l = 0; l < m; l++)
            { // loop on suffix
              generator_(n + i * m + l, n + j * m + l) = gk(i, j) * Vrate_[k - 1];
            }
            n += m * vsize[k - 1];
          }
        }
      }
    }
    m *= vsize[k - 1];
  }

  // modification of generator_

  this->completeMatrices();

  double x;

  for (i = 0; i < salph; i++)
  {
    x = 0;
    for (j = 0; j < salph; j++)
    {
      if (j != i)
        x += generator_(i, j);
    }
    generator_(i, i) = -x;
  }

  // at that point generator_ and freq_ are done for models without
  // enableEigenDecomposition

  // Eigen values:
  
  if (enableEigenDecomposition())
  {
    for (i = 0; i < salph; i++)
    {
      bool flag = true;
      for (j = 0; j < salph; j++)
      {
        if ((i != j) && abs(generator_(i, j)) > NumConstants::TINY())
        {
          flag = false;
          break;
        }
      }
      if (flag)
        nbStop++;
      vnull.push_back(flag);
    }

    if (nbStop != 0)
    {
      size_t gi = 0, gj = 0;

      gk.resize(salph - nbStop, salph - nbStop);
      for (i = 0; i < salph; i++)
      {
        if (!vnull[i])
        {
          gj = 0;
          for (j = 0; j < salph; j++)
          {
            if (!vnull[j])
            {
              gk(i - gi, j - gj) = generator_(i, j);
            }
            else
              gj++;
          }
        }
        else
          gi++;
      }

      EigenValue<double> ev(gk);
      eigenValues_ = ev.getRealEigenValues();
      iEigenValues_ = ev.getImagEigenValues();

      for (i = 0; i < nbStop; i++)
      {
        eigenValues_.push_back(0);
        iEigenValues_.push_back(0);
      }

      RowMatrix<double> rev = ev.getV();
      rightEigenVectors_.resize(salph, salph);
      gi = 0;
      for (i = 0; i < salph; i++)
      {
        if (vnull[i])
        {
          gi++;
          for (j = 0; j < salph; j++)
          {
            rightEigenVectors_(i, j) = 0;
          }

          rightEigenVectors_(i, salph - nbStop + gi - 1) = 1;
        }
        else
        {
          for (j = 0; j < salph - nbStop; j++)
          {
            rightEigenVectors_(i, j) = rev(i - gi, j);
          }

          for (j = salph - nbStop; j < salph; j++)
          {
            rightEigenVectors_(i, j) = 0;
          }
        }
      }
    }
    else
    {
      EigenValue<double> ev(generator_);
      eigenValues_ = ev.getRealEigenValues();
      iEigenValues_ = ev.getImagEigenValues();
      rightEigenVectors_ = ev.getV();
      nbStop = 0;
    }

    try
    {
      MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);

      // is it diagonalizable ?

      isDiagonalizable_ = true;
      for (i = 0; i < size_ && isDiagonalizable_; i++)
      {
        if (abs(iEigenValues_[i]) > NumConstants::SMALL())
          isDiagonalizable_ = false;
      }

      // is it singular?

      // looking for the 0 eigenvector for which the non-stop right
      // eigen vector elements are equal.
      //

      size_t nulleigen = 0;
      double val;

      isNonSingular_ = false;
      while (nulleigen < salph - nbStop)
      {
        if ((abs(eigenValues_[nulleigen]) < NumConstants::SMALL()) && (abs(iEigenValues_[nulleigen]) < NumConstants::SMALL()))
        {
          i = 0;
          while (vnull[i])
            i++;
          
          val = rightEigenVectors_(i, nulleigen);
          i++;
          while (i < salph)
          {
            if (!vnull[i])
            {
              if (abs(rightEigenVectors_(i, nulleigen) - val) > NumConstants::SMALL())
                break;
            }
            i++;
          }
          
          if (i < salph)
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
        
        for (i = 0; i < salph; i++)
          freq_[i] = leftEigenVectors_(nulleigen, i);
        
        x = 0;
        for (i = 0; i < salph; i++)
            x += freq_[i];
        
        for (i = 0; i < salph; i++)
          freq_[i] /= x;
      }
      
      else
      {
        ApplicationTools::displayMessage("Unable to find eigenvector for eigenvalue 1. Taylor series used instead.");
        isDiagonalizable_ = false;
      }
    }
    
    
    // if rightEigenVectors_ is singular
    catch (ZeroDivisionException& e)
    {
      ApplicationTools::displayMessage("Singularity during  diagonalization. Taylor series used instead.");
      isNonSingular_ = false;
      isDiagonalizable_ = false;
    }

    if (!isNonSingular_)
    {
      double min = generator_(0, 0);
      for (i = 1; i < salph; i++)
      {
        if (min > generator_(i, i))
          min = generator_(i, i);
      }

      MatrixTools::scale(generator_, -1 / min);

      if (vPowGen_.size() == 0)
        vPowGen_.resize(30);

      MatrixTools::getId(salph, tmpMat_);    // to compute the equilibrium frequency  (Q+Id)^256
      MatrixTools::add(tmpMat_, generator_);
      MatrixTools::pow(tmpMat_, 256, vPowGen_[0]);

      for (i = 0; i < salph; i++)
      {
        freq_[i] = vPowGen_[0](0, i);
      }

      MatrixTools::getId(salph, vPowGen_[0]);
    }

    // normalization

    x = 0;
    for (i = 0; i < salph; i++)
      x += freq_[i] * generator_(i, i);

    MatrixTools::scale(generator_, -1. / x);
    for (i = 0; i < salph; i++)
    {
      eigenValues_[i] /= -x;
      iEigenValues_[i] /= -x;
    }

    if (!isNonSingular_)
      MatrixTools::Taylor(generator_, 30, vPowGen_);
  }
  else  // compute freq_ is no eigenDecomposition
  {
    for (j = 0; j < size_; j++)
      freq_[j] = 1;
  
    m = 1;
    for (k = nbmod; k > 0; k--)
    {
      SubstitutionModel* pSM = VSubMod_[k - 1];
      for (j = 0; j < vsize[k - 1]; j++)
      {
        n = 0;
        while (n < salph)
        { // loop on prefix
          for (l = 0; l < m; l++)
          { // loop on suffix
            freq_[n + j * m + l] *=  pSM->freq(j);
          }
          n += m * vsize[k - 1];
        }
      }
      m *= vsize[k - 1];
    }
  }
  

  // compute the exchangeability_

  for (i = 0; i < size_; i++)
    for (j = 0; j < size_; j++)
      exchangeability_(i, j) = generator_(i, j) / freq_[j];
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
