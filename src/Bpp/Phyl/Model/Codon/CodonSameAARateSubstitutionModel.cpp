// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "CodonSameAARateSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

CodonSameAARateSubstitutionModel::CodonSameAARateSubstitutionModel(
  unique_ptr<ProteinSubstitutionModelInterface> pAAmodel,
  unique_ptr<CodonSubstitutionModelInterface> pCodonModel,
  unique_ptr<CodonFrequencySetInterface> pFreq,
  shared_ptr<const GeneticCode> pgencode) :
  AbstractParameterAliasable("SameAARate."),
  AbstractSubstitutionModel(pCodonModel->getAlphabet(), pCodonModel->getStateMap(), "SameAARate."),
  pAAmodel_(std::move(pAAmodel)),
  pCodonModel_(std::move(pCodonModel)),
  pFreq_(std::move(pFreq)),
  pgencode_(pgencode),
  X_(20, 20),
  phi_(20)
{
  pCodonModel_->enableEigenDecomposition(true);
  pAAmodel_->enableEigenDecomposition(true);

  pAAmodel_->setNamespace("SameAARate." + pAAmodel_->getNamespace());
  pCodonModel_->setNamespace("SameAARate." + pCodonModel_->getNamespace());

  //TODO (jdutheil on 30/12/22): if we want this, we need to use shared_ptr for FrequencySets
  if (pFreq_ && pFreq_.get() != &(dynamic_cast<const CoreCodonSubstitutionModelInterface*>(pCodonModel_.get())->codonFrequencySet()))
    pFreq_->setNamespace("SameAARate." + pFreq_->getNamespace());

  addParameters_(pAAmodel_->getParameters());
  addParameters_(pCodonModel_->getParameters());

  //TODO (jdutheil on 30/12/22): if we want this, we need to use shared_ptr for FrequencySets
  if (pFreq_ && pFreq_.get() != &(dynamic_cast<const CoreCodonSubstitutionModelInterface*>(pCodonModel_.get())->codonFrequencySet()))
    addParameters_(pFreq_->getParameters());

  compute_();
  updateMatrices_();
}

void CodonSameAARateSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  pAAmodel_->matchParametersValues(parameters);
  pCodonModel_->matchParametersValues(parameters);

  if (pFreq_)
    pFreq_->matchParametersValues(parameters);

  compute_();
  updateMatrices_();
}

void CodonSameAARateSubstitutionModel::compute_()
{
  const auto& freq = pFreq_ ? pFreq_->getFrequencies() : pCodonModel_->getFrequencies();

  const auto& gen = pCodonModel_->generator();

  std::fill(phi_.begin(), phi_.end(), 0);

  for (size_t i = 0; i < stateMap_->getNumberOfModelStates(); i++)
  {
    if (!pgencode_->isStop((int)i))
      phi_[pAAmodel_->getModelStates(pgencode_->translate((int)i))[0]] += freq[i];
  }


  for (size_t ai = 0; ai < 20; ai++)
  {
    std::fill(X_.getRow(ai).begin(), X_.getRow(ai).end(), 0);
  }

  for (size_t i = 0; i < stateMap_->getNumberOfModelStates(); i++)
  {
    if (pgencode_->isStop((int)i))
      continue;
    auto ai = pAAmodel_->getModelStates(pgencode_->translate((int)i))[0];

    auto& X_i = X_.getRow(ai);

    for (size_t j = 0; j < stateMap_->getNumberOfModelStates(); j++)
    {
      if (pgencode_->isStop((int)j))
        continue;
      auto aj = pAAmodel_->getModelStates(pgencode_->translate((int)j))[0];
      X_i[aj] += gen(i, j) * freq[i];
    }
  }

  for (size_t ai = 0; ai < 20; ai++)
  {
    auto& X_i = X_.getRow(ai);
    for (size_t aj = 0; aj < 20; aj++)
    {
      if (X_i[aj] != 0)
        X_i[aj] = phi_[ai] * pAAmodel_->generator()(ai, aj) / X_i[aj];
    }
  }


  for (size_t i = 0; i < stateMap_->getNumberOfModelStates(); i++)
  {
    if (pgencode_->isStop((int)i))
      continue;
    auto ai = pAAmodel_->getModelStates(pgencode_->translate((int)i))[0];

    const auto& X_i = X_.row(ai);

    for (size_t j = 0; j < stateMap_->getNumberOfModelStates(); j++)
    {
      if (pgencode_->isStop((int)j))
        continue;
      auto aj = pAAmodel_->getModelStates(pgencode_->translate((int)j))[0];
      generator_(i, j) = gen(i, j) * X_i[aj];
    }
  }

  // Set diagonal

  setDiagonal();
}
