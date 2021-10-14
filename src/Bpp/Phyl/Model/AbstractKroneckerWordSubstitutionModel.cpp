//
// File: AbstractKroneckerWordSubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: mardi 26 juillet 2016, ÃÂ  10h 56
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
  const std::vector<std::set<size_t> >& vPos,
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
  SubstitutionModel* pmodel,
  unsigned int num,
  const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractWordSubstitutionModel(pmodel, num, prefix),
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
  SubstitutionModel* pmodel,
  unsigned int num,
  const std::vector<std::set<size_t> >& vPos,
  const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractWordSubstitutionModel(pmodel, num, prefix),
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
  const Alphabet* alph,
  std::shared_ptr<const StateMap> stateMap,
  const std::string& prefix) :
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


void AbstractKroneckerWordSubstitutionModel::setChangingPositions(const std::vector<std::set<size_t> >& vPos)
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

void AbstractKroneckerWordSubstitutionModel::fillBasicGenerator()
{
  size_t nbmod = VSubMod_.size();

// Generator

  if (sChangingPos_.size() == 0)
  {
    vGenerators_[0] = VSubMod_[0]->getGenerator();

    for (size_t k = 1; k < nbmod - 1; k++)
    {
      MatrixTools::kroneckerMult(vGenerators_[k - 1], VSubMod_[k]->getGenerator(), 1., 1., vGenerators_[k], false);
    }

    MatrixTools::kroneckerMult(vGenerators_[nbmod - 2], VSubMod_[nbmod - 1]->getGenerator(), 1., 1., generator_, false);
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
            MatrixTools::kroneckerMult(vGenerators_[posok - 1], VSubMod_[posok]->getGenerator(), 1., 0., vGenerators_[posok], false);
          }
          else
            vGenerators_[posok] = VSubMod_[0]->getGenerator();

          MatrixTools::fillDiag(vGenerators_[posok], 0.);
        }
        else
        {
          while (pos < posok)
          {
            MatrixTools::kroneckerMult(vGenerators_[pos - 1], VSubMod_[pos]->getNumberOfStates(), 1., vGenerators_[pos], false);
            pos++;
          }
          MatrixTools::kroneckerMult(vGenerators_[posok - 1], VSubMod_[posok]->getGenerator(), 0., 0., vGenerators_[posok], false);
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
