// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "EquiprobableSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

#include <cmath>
#include <map>

using namespace std;

/******************************************************************************/

EquiprobableSubstitutionModel::EquiprobableSubstitutionModel(
    std::shared_ptr<const Alphabet> alpha) :
  AbstractParameterAliasable("Equi."),
  AbstractReversibleSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "Equi."),
  exp_(), p_(size_, size_), freqSet_(nullptr)
{
  freqSet_ = make_unique<FixedFrequencySet>(getStateMap(), freq_);
  updateMatrices_();
}

EquiprobableSubstitutionModel::EquiprobableSubstitutionModel(
    std::shared_ptr<const Alphabet> alpha,
    std::unique_ptr<FrequencySetInterface> freqSet,
    bool initFreqs) :
  AbstractParameterAliasable("Equi+F."),
  AbstractReversibleSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), "Equi+F."),
  exp_(), p_(size_, size_), freqSet_(std::move(freqSet))
{
  freqSet_->setNamespace("Equi+F." + freqSet_->getNamespace());
  if (initFreqs)
    freqSet_->setFrequencies(freq_);
  else
    freq_ = freqSet_->getFrequencies();
  addParameters_(freqSet_->getParameters());
  updateMatrices_();
}

/******************************************************************************/

void EquiprobableSubstitutionModel::updateMatrices_()
{
  // Frequencies:
  for (unsigned int i = 0; i < size_; ++i)
  {
    freq_[i] = 1. / (float)size_;
  }

  // Generator:
  for (unsigned int i = 0; i < size_; ++i)
  {
    for (unsigned int j = 0; j < size_; ++j)
    {
      generator_(i, j) = (i == j) ? -1. : 1. / (float)(size_ - 1);
      exchangeability_(i, j) = generator_(i, j) * (float)size_;
    }
  }

  // Eigen values:
  eigenValues_[0] = 0;
  for (unsigned int i = 1; i < size_; ++i)
  {
    eigenValues_[i] = -(float)size_ / (float)(size_ - 1);
  }

  // Eigen vectors:
  for (unsigned int i = 0; i < size_; ++i)
  {
    leftEigenVectors_(0, i) = 1. / (float)size_;
  }
  for (unsigned int i = 1; i < size_; ++i)
  {
    for (unsigned int j = 0; j < size_; ++j)
    {
      leftEigenVectors_(i, j) = -1. / (float)size_;
    }
  }
  for (unsigned int i = 0; i < (size_ - 1); ++i)
  {
    leftEigenVectors_((size_ - 1) - i, i) = (float)(size_ - 1) / (float)size_;
  }

  for (unsigned int i = 0; i < size_; ++i)
  {
    rightEigenVectors_(i, 0) = 1.;
  }
  for (unsigned int i = 1; i < size_; ++i)
  {
    rightEigenVectors_((size_ - 1), i) = -1.;
  }
  for (unsigned int i = 0; i < (size_ - 1); ++i)
  {
    for (unsigned int j = 1; j < size_; ++j)
    {
      rightEigenVectors_(i, j) = 0.;
    }
  }
  for (unsigned int i = 1; i < size_; ++i)
  {
    rightEigenVectors_((size_ - 1) - i, i) = 1.;
  }
}

/******************************************************************************/

double EquiprobableSubstitutionModel::Pij_t(size_t i, size_t j, double d) const
{
  if (i == j)
    return 1. / (float)size_ + (float)(size_ - 1) / (float)size_ * exp(-rate_ * (float)size_ / (float)(size_ - 1) * d);
  else
    return 1. / (float)size_ -  1. / (float)size_ * exp(-rate_ * (float)size_ / (float)(size_ - 1) * d);
}

/******************************************************************************/

double EquiprobableSubstitutionModel::dPij_dt(size_t i, size_t j, double d) const
{
  if (i == j)
    return -rate_*        exp(-rate_ * (float)size_ / (float)(size_ - 1) * d);
  else
    return rate_ * 1. / (float)(size_ - 1) * exp(-rate_ * (float)size_ / (float)(size_ - 1) * d);
}

/******************************************************************************/

double EquiprobableSubstitutionModel::d2Pij_dt2(size_t i, size_t j, double d) const
{
  if (i == j)
    return rate_ *  rate_ * (float)size_ / (float)(size_ - 1)  * exp(-rate_ * (float)size_ / (float)(size_ - 1) * d);
  else
    return -rate_ *  rate_ * (float)size_ / pow((float)(size_ - 1), 2) * exp(-rate_ * (float)size_ / (float)(size_ - 1) * d);
}

/******************************************************************************/

const Matrix<double>& EquiprobableSubstitutionModel::getPij_t(double d) const
{
  exp_ = exp(-rate_ * (float)size_ / (float)(size_ - 1) * d);
  for (unsigned int i = 0; i < size_; i++)
  {
    for (unsigned int j = 0; j < size_; j++)
    {
      p_(i, j) = (i == j) ? 1. / (float)size_ + (float)(size_ - 1) / (float)size_ * exp_ : 1. / (float)size_ - 1. / (float)size_ * exp_;
    }
  }
  return p_;
}

const Matrix<double>& EquiprobableSubstitutionModel::getdPij_dt(double d) const
{
  exp_ = exp(-rate_ * (float)size_ / (float)(size_ - 1) * d);
  for (unsigned int i = 0; i < size_; i++)
  {
    for (unsigned int j = 0; j < size_; j++)
    {
      p_(i, j) =  rate_ * ((i == j) ? -exp_ : 1. / (float)(size_ - 1) * exp_);
    }
  }
  return p_;
}

const Matrix<double>& EquiprobableSubstitutionModel::getd2Pij_dt2(double d) const
{
  exp_ = exp( rate_ * -(float)size_ / (float)(size_ - 1) * d);
  for (unsigned int i = 0; i < size_; i++)
  {
    for (unsigned int j = 0; j < size_; j++)
    {
      p_(i, j) =  rate_ *  rate_ * ((i == j) ? (float)size_ / (float)(size_ - 1) * exp_ : -(float)size_ / pow((float)(size_ - 1), 2) * exp_);
    }
  }
  return p_;
}

/******************************************************************************/

void EquiprobableSubstitutionModel::setFreq(std::map<int, double>& freqs)
{
  for (auto i : freqs)
  {
    freq_[(size_t)i.first] = i.second;
  }

  freqSet_->setFrequencies(freq_);
  // Update parameters and re-compute generator and eigen values:
  matchParametersValues(freqSet_->getParameters());
}

/******************************************************************************/
