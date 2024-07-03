// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractKroneckerCodonSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

AbstractKroneckerCodonSubstitutionModel::AbstractKroneckerCodonSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    const string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      std::shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      prefix),
  gCode_(gCode)
{
  enableEigenDecomposition(true);

  size_t i;
  for (i = 0; i < 3; i++)
  {
    VSubMod_.push_back(std::move(pmod));
    VnestedPrefix_.push_back(VSubMod_[i]->getNamespace());
  }

  VSubMod_[0]->setNamespace(prefix + "123_" + VnestedPrefix_[0]);
  VSubMod_[0]->enableEigenDecomposition(0);
  addParameters_(VSubMod_[0]->getParameters());

  initGenerators_();
}

AbstractKroneckerCodonSubstitutionModel::AbstractKroneckerCodonSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod,
    const std::vector<std::set< size_t>>& vPos,
    const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      std::shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      prefix),
  gCode_(gCode)
{
  enableEigenDecomposition(true);

  shared_ptr<NucleotideSubstitutionModelInterface> pmod2 = std::move(pmod);
  size_t i;
  for (i = 0; i < 3; i++)
  {
    VSubMod_.push_back(pmod2);
    VnestedPrefix_.push_back(pmod2->getNamespace());
  }
  pmod2->setNamespace(prefix + "123_" + VnestedPrefix_[0]);
  pmod2->enableEigenDecomposition(0);
  addParameters_(pmod2->getParameters());

  initGenerators_();
  setChangingPositions(vPos);
}

AbstractKroneckerCodonSubstitutionModel::AbstractKroneckerCodonSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    const string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      prefix),
  gCode_(gCode)
{
  enableEigenDecomposition(true);

  VSubMod_.push_back(std::move(pmod1));
  VnestedPrefix_.push_back(VSubMod_[0]->getNamespace());
  VSubMod_[0]->setNamespace(prefix + "1_" + VnestedPrefix_[0]);
  VSubMod_[0]->enableEigenDecomposition(0);
  addParameters_(VSubMod_[0]->getParameters());

  VSubMod_.push_back(std::move(pmod2));
  VnestedPrefix_.push_back(VSubMod_[1]->getNamespace());
  VSubMod_[1]->setNamespace(prefix + "2_" + VnestedPrefix_[1]);
  VSubMod_[1]->enableEigenDecomposition(0);
  addParameters_(VSubMod_[1]->getParameters());

  VSubMod_.push_back(std::move(pmod3));
  VnestedPrefix_.push_back(VSubMod_[2]->getNamespace());
  VSubMod_[2]->setNamespace(prefix + "3_" + VnestedPrefix_[2]);
  VSubMod_[2]->enableEigenDecomposition(0);
  addParameters_(VSubMod_[2]->getParameters());

  initGenerators_();
}

AbstractKroneckerCodonSubstitutionModel::AbstractKroneckerCodonSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
    std::unique_ptr<NucleotideSubstitutionModelInterface> pmod3,
    const std::vector<std::set< size_t>>& vPos,
    const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractKroneckerWordSubstitutionModel(
      gCode->getSourceAlphabet(),
      shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      prefix),
  gCode_(gCode)
{
  enableEigenDecomposition(true);

  VSubMod_.push_back(std::move(pmod1));
  VnestedPrefix_.push_back(VSubMod_[0]->getNamespace());
  VSubMod_[0]->setNamespace(prefix + "1_" + VnestedPrefix_[0]);
  VSubMod_[0]->enableEigenDecomposition(0);
  addParameters_(VSubMod_[0]->getParameters());

  VSubMod_.push_back(std::move(pmod2));
  VnestedPrefix_.push_back(VSubMod_[1]->getNamespace());
  VSubMod_[1]->setNamespace(prefix + "2_" + VnestedPrefix_[1]);
  VSubMod_[1]->enableEigenDecomposition(0);
  addParameters_(VSubMod_[1]->getParameters());

  VSubMod_.push_back(std::move(pmod3));
  VnestedPrefix_.push_back(VSubMod_[2]->getNamespace());
  VSubMod_[2]->setNamespace(prefix + "3_" + VnestedPrefix_[2]);
  VSubMod_[2]->enableEigenDecomposition(0);
  addParameters_(VSubMod_[2]->getParameters());

  initGenerators_();
  setChangingPositions(vPos);
}


void AbstractKroneckerCodonSubstitutionModel::completeMatrices_()
{
  size_t i, j;
  size_t salph = getNumberOfStates();

  for (i = 0; i < salph; i++)
  {
    for (j = 0; j < salph; j++)
    {
      if (gCode_->isStop(static_cast<int>(i)) || gCode_->isStop(static_cast<int>(j)))
      {
        generator_(i, j) = 0;
      }
      else
        generator_(i, j) *= getCodonsMulRate(i, j);
    }
  }
}
