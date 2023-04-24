//
// File: CodonSameAARateSubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: vendredi 30 octobre 2020, ÃÂ  17h 43
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
  pAAmodel_(move(pAAmodel)),
  pCodonModel_(move(pCodonModel)),
  pFreq_(move(pFreq)),
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

  const auto& gen = pCodonModel_->getGenerator();

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
        X_i[aj] = phi_[ai] * pAAmodel_->getGenerator()(ai, aj) / X_i[aj];
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
