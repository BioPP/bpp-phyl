// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "JCprot.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

#include <cmath>
#include <map>

using namespace std;

/******************************************************************************/

JCprot::JCprot(std::shared_ptr<const ProteicAlphabet> alpha) :
  AbstractParameterAliasable("JC69."),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "JC69."),
  exp_(), p_(size_, size_), freqSet_(nullptr), withFreq_(false)
{
  freqSet_.reset(new FixedProteinFrequencySet(alpha));
  updateMatrices_();
}

JCprot::JCprot(
    std::shared_ptr<const ProteicAlphabet> alpha,
    std::unique_ptr<ProteinFrequencySetInterface> freqSet,
    bool initFreqs) :
  AbstractParameterAliasable("JC69+F."),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "JC69+F."),
  exp_(), p_(size_, size_), freqSet_(std::move(freqSet)), withFreq_(true)
{
  freqSet_->setNamespace("JC69+F." + freqSet_->getNamespace());
  if (initFreqs)
    freqSet_->setFrequencies(freq_);
  else
    freq_ = freqSet_->getFrequencies();
  addParameters_(freqSet_->getParameters());

  updateMatrices_();
}

/******************************************************************************/

void JCprot::updateMatrices_()
{
  for (unsigned int i = 0; i < 20; ++i)
  {
    for (unsigned int j = 0; j < 20; ++j)
    {
      exchangeability_(i, j) = (i == j) ? -20. : 20. / 19.;
    }
  }

  if (!withFreq_)
  {
    // Frequencies:
    for (unsigned int i = 0; i < 20; ++i)
    {
      freq_[i] = 1. / 20.;
    }

    // Generator:
    for (unsigned int i = 0; i < 20; ++i)
    {
      for (unsigned int j = 0; j < 20; ++j)
      {
        generator_(i, j) = (i == j) ? -1. : 1. / 19.;
      }
    }

    // Eigen values:
    eigenValues_[0] = 0;
    for (unsigned int i = 1; i < 20; ++i)
    {
      eigenValues_[i] = -20. / 19.;
    }

    // Eigen vectors:
    for (unsigned int i = 0; i < 20; ++i)
    {
      leftEigenVectors_(0, i) = 1. / 20.;
    }
    for (unsigned int i = 1; i < 20; ++i)
    {
      for (unsigned int j = 0; j < 20; ++j)
      {
        leftEigenVectors_(i, j) = -1. / 20.;
      }
    }
    for (unsigned int i = 0; i < 19; ++i)
    {
      leftEigenVectors_(19 - i, i) = 19. / 20.;
    }

    for (unsigned int i = 0; i < 20; ++i)
    {
      rightEigenVectors_(i, 0) = 1.;
    }
    for (unsigned int i = 1; i < 20; ++i)
    {
      rightEigenVectors_(19, i) = -1.;
    }
    for (unsigned int i = 0; i < 19; ++i)
    {
      for (unsigned int j = 1; j < 20; ++j)
      {
        rightEigenVectors_(i, j) = 0.;
      }
    }
    for (unsigned int i = 1; i < 20; ++i)
    {
      rightEigenVectors_(19 - i, i) = 1.;
    }
  }
  else
    AbstractReversibleSubstitutionModel::updateMatrices_();
}

/******************************************************************************/

double JCprot::Pij_t(size_t i, size_t j, double d) const
{
  if (!withFreq_)
  {
    if (i == j)
      return 1. / 20. + 19. / 20. * exp(-rate_ * 20. / 19. * d);
    else
      return 1. / 20. -  1. / 20. * exp(-rate_ * 20. / 19. * d);
  }
  else
    return AbstractSubstitutionModel::Pij_t(i, j, d);
}

/******************************************************************************/

double JCprot::dPij_dt(size_t i, size_t j, double d) const
{
  if (!withFreq_)
  {
    if (i == j)
      return -rate_*        exp(-rate_ * 20. / 19. * d);
    else
      return rate_ * 1. / 19. * exp(-rate_ * 20. / 19. * d);
  }
  else
    return AbstractSubstitutionModel::dPij_dt(i, j, d);
}

/******************************************************************************/

double JCprot::d2Pij_dt2(size_t i, size_t j, double d) const
{
  if (!withFreq_)
  {
    if (i == j)
      return rate_ *  rate_ * 20. / 19.  * exp(-rate_ * 20. / 19. * d);
    else
      return -rate_ *  rate_ * 20. / 361. * exp(-rate_ * 20. / 19. * d);
  }
  else
    return AbstractSubstitutionModel::d2Pij_dt2(i, j, d);
}

/******************************************************************************/

const Matrix<double>& JCprot::getPij_t(double d) const
{
  if (!withFreq_)
  {
    exp_ = exp(-rate_ * 20. / 19. * d);
    for (unsigned int i = 0; i < size_; i++)
    {
      for (unsigned int j = 0; j < size_; j++)
      {
        p_(i, j) = (i == j) ? 1. / 20. + 19. / 20. * exp_ : 1. / 20. - 1. / 20. * exp_;
      }
    }
    return p_;
  }
  else
    return AbstractSubstitutionModel::getPij_t(d);
}

const Matrix<double>& JCprot::getdPij_dt(double d) const
{
  if (!withFreq_)
  {
    exp_ = exp(-rate_ * 20. / 19. * d);
    for (unsigned int i = 0; i < size_; i++)
    {
      for (unsigned int j = 0; j < size_; j++)
      {
        p_(i, j) =  rate_ * ((i == j) ? -exp_ : 1. / 19. * exp_);
      }
    }
    return p_;
  }
  else
    return AbstractSubstitutionModel::getdPij_dt(d);
}

const Matrix<double>& JCprot::getd2Pij_dt2(double d) const
{
  if (!withFreq_)
  {
    exp_ = exp( rate_ * -20. / 19. * d);
    for (unsigned int i = 0; i < size_; i++)
    {
      for (unsigned int j = 0; j < size_; j++)
      {
        p_(i, j) =  rate_ *  rate_ * ((i == j) ? 20. / 19. * exp_ : -20. / 361. * exp_);
      }
    }
    return p_;
  }
  else
    return AbstractSubstitutionModel::getd2Pij_dt2(d);
}

/******************************************************************************/

void JCprot::setFreqFromData(const SequenceDataInterface& data, double pseudoCount)
{
  map<int, double> counts;
  SequenceContainerTools::getFrequencies(data, counts, pseudoCount);
  for (auto i : counts)
  {
    freq_[(size_t)i.first] = i.second;
  }

  freqSet_->setFrequencies(freq_);
  // Update parameters and re-compute generator and eigen values:
  matchParametersValues(freqSet_->getParameters());
}

/******************************************************************************/
