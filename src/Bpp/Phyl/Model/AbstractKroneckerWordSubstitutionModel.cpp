// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>

#include "AbstractKroneckerWordSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractKroneckerWordSubstitutionModel::AbstractKroneckerWordSubstitutionModel(
    ModelList& modelList,
    const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractWordSubstitutionModel(modelList, prefix),
  sChangingPos_(),
  vGenerators_()
{
  enableEigenDecomposition(true);
  for (auto& submod : VSubMod_)
  {
    submod->enableEigenDecomposition(false);
  }

  initGenerators_();
}

AbstractKroneckerWordSubstitutionModel::AbstractKroneckerWordSubstitutionModel(
    ModelList& modelList,
    const std::vector<std::set<size_t>>& vPos,
    const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractWordSubstitutionModel(modelList, prefix),
  sChangingPos_(vPos),
  vGenerators_()
{
  if (!checkChangingPositions_())
    throw Exception("AbstractKroneckerWordSubstitutionModel::AbstractKroneckerWordSubstitutionModel: Bad set of changing positions ");

  enableEigenDecomposition(true);
  for (auto& submod : VSubMod_)
  {
    submod->enableEigenDecomposition(false);
  }

  initGenerators_();
}

AbstractKroneckerWordSubstitutionModel::AbstractKroneckerWordSubstitutionModel(
    unique_ptr<SubstitutionModelInterface> pmodel,
    unsigned int num,
    const string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractWordSubstitutionModel(std::move(pmodel), num, prefix),
  sChangingPos_(),
  vGenerators_()
{
  enableEigenDecomposition(true);
  for (auto& submod : VSubMod_)
  {
    submod->enableEigenDecomposition(false);
  }

  initGenerators_();
}

AbstractKroneckerWordSubstitutionModel::AbstractKroneckerWordSubstitutionModel(
    unique_ptr<SubstitutionModelInterface> pmodel,
    unsigned int num,
    const std::vector<std::set<size_t>>& vPos,
    const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractWordSubstitutionModel(std::move(pmodel), num, prefix),
  sChangingPos_(vPos),
  vGenerators_()
{
  if (!checkChangingPositions_())
    throw Exception("AbstractKroneckerWordSubstitutionModel::AbstractKroneckerWordSubstitutionModel: Bad set of changing positions ");

  enableEigenDecomposition(true);
  for (auto& submod : VSubMod_)
  {
    submod->enableEigenDecomposition(false);
  }

  initGenerators_();
}

AbstractKroneckerWordSubstitutionModel::AbstractKroneckerWordSubstitutionModel(
    std::shared_ptr<const Alphabet> alph,
    std::shared_ptr<const StateMapInterface> stateMap,
    const string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractWordSubstitutionModel(alph, stateMap, prefix),
  sChangingPos_(),
  vGenerators_()
{
  enableEigenDecomposition(true);
  for (auto& submod : VSubMod_)
  {
    submod->enableEigenDecomposition(false);
  }
}

AbstractKroneckerWordSubstitutionModel::AbstractKroneckerWordSubstitutionModel(
    const AbstractKroneckerWordSubstitutionModel& wrsm) :
  AbstractParameterAliasable(wrsm),
  AbstractWordSubstitutionModel(wrsm),
  sChangingPos_(wrsm.sChangingPos_),
  vGenerators_(wrsm.vGenerators_)
{
  initGenerators_();
}


AbstractKroneckerWordSubstitutionModel& AbstractKroneckerWordSubstitutionModel::operator=(
    const AbstractKroneckerWordSubstitutionModel& model)
{
  AbstractParameterAliasable::operator=(model);
  AbstractWordSubstitutionModel::operator=(model);
  sChangingPos_ = model.sChangingPos_;
  vGenerators_ = model.vGenerators_;

  return *this;
}


void AbstractKroneckerWordSubstitutionModel::setChangingPositions(const vector<std::set<size_t>>& vPos)
{
  sChangingPos_ = vPos;

  if (!checkChangingPositions_())
    throw Exception("AbstractKroneckerWordSubstitutionModel::setChangingPositions: Bad set of changing positions ");
}


void AbstractKroneckerWordSubstitutionModel::initGenerators_()
{
  vGenerators_.clear();

  size_t nbmod = VSubMod_.size();

  size_t size = 1;

  for (size_t k = 0; k < nbmod; k++)
  {
    size *= VSubMod_[k]->getNumberOfStates();
    vGenerators_.push_back(RowMatrix<double>(size, size));
  }
}

bool AbstractKroneckerWordSubstitutionModel::checkChangingPositions_()
{
  size_t nbmod = VSubMod_.size();

  for (const auto& i : sChangingPos_)
  {
    if (*(--(i.end())) > nbmod || *(i.begin()) == 0)
      return false;
  }

  return true;
}

void AbstractKroneckerWordSubstitutionModel::fillBasicGenerator_()
{
  size_t nbmod = VSubMod_.size();

// Generator

  if (sChangingPos_.size() == 0)
  {
    vGenerators_[0] = VSubMod_[0]->generator();

    for (size_t k = 1; k < nbmod - 1; k++)
    {
      MatrixTools::kroneckerMult(vGenerators_[k - 1], VSubMod_[k]->generator(), 1., 1., vGenerators_[k], false);
    }

    MatrixTools::kroneckerMult(vGenerators_[nbmod - 2], VSubMod_[nbmod - 1]->generator(), 1., 1., generator_, false);
  }
  else
  {
    bool begin(true);

    for (const auto& i : sChangingPos_)
    {
      size_t pos = 0;

      for (const auto& iPos : i)
      {
        size_t posok = iPos - 1; // position of the next generator to multiply

        if (pos == 0)
        {
          size_t ss = 1;
          while (pos < posok)
          {
            ss *= VSubMod_[pos]->getNumberOfStates();
            pos++;
          }
          if (ss != 1)
          {
            MatrixTools::getId(ss, vGenerators_[posok - 1]);
            MatrixTools::kroneckerMult(vGenerators_[posok - 1], VSubMod_[posok]->generator(), 1., 0., vGenerators_[posok], false);
          }
          else
            vGenerators_[posok] = VSubMod_[0]->generator();

          MatrixTools::fillDiag(vGenerators_[posok], 0.);
        }
        else
        {
          while (pos < posok)
          {
            MatrixTools::kroneckerMult(vGenerators_[pos - 1], VSubMod_[pos]->getNumberOfStates(), 1., vGenerators_[pos], false);
            pos++;
          }
          MatrixTools::kroneckerMult(vGenerators_[posok - 1], VSubMod_[posok]->generator(), 0., 0., vGenerators_[posok], false);
        }

        pos++;
      }

      while (pos < nbmod)
      {
        MatrixTools::kroneckerMult(vGenerators_[pos - 1], VSubMod_[pos]->getNumberOfStates(), 1., vGenerators_[pos], false);
        pos++;
      }

      if (begin)
      {
        begin = false;
        generator_ = vGenerators_[nbmod - 1];
      }
      else
        MatrixTools::add(generator_, vGenerators_[nbmod - 1]);
    }
  }
}
