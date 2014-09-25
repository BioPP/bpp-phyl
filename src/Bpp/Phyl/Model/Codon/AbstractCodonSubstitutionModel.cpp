//
// File: CodonSubstitutionModel.cpp
// Created by:  Laurent Gueguen
// Created on: Feb 2009
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

#include "AbstractCodonSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractCodonSubstitutionModel::AbstractCodonSubstitutionModel(
  const GeneticCode* gCode,
  NucleotideSubstitutionModel* pmod,
  const std::string& prefix,
  bool paramRates) :
  AbstractParameterAliasable(prefix),
  AbstractWordSubstitutionModel(gCode->getSourceAlphabet(), new CanonicalStateMap(gCode->getSourceAlphabet(), false), prefix),
  hasParametrizedRates_(paramRates),
  gCode_(gCode)
{
  enableEigenDecomposition(true);

  size_t i;
  for (i = 0; i < 3; i++)
  {
    VSubMod_.push_back(pmod);
    VnestedPrefix_.push_back(pmod->getNamespace());
  }

  pmod->setNamespace(prefix + "123_" + VnestedPrefix_[0]);
  pmod->enableEigenDecomposition(0);
  addParameters_(pmod->getParameters());

  Vrate_.resize(3);
  for (i = 0; i < 3; i++)
  {
    Vrate_[i] = 1.0 / 3;
  }


  if (hasParametrizedRates_)
  {
    // relative rates
    for (i = 0; i < 2; i++)
    {
      addParameter_(new Parameter(prefix + "relrate" + TextTools::toString(i + 1), 1.0 / static_cast<double>(3 - i), &Parameter::PROP_CONSTRAINT_EX));
    }
  }
}

AbstractCodonSubstitutionModel::AbstractCodonSubstitutionModel(
  const GeneticCode* gCode,
  NucleotideSubstitutionModel* pmod1,
  NucleotideSubstitutionModel* pmod2,
  NucleotideSubstitutionModel* pmod3,
  const std::string& prefix,
  bool paramRates) :
  AbstractParameterAliasable(prefix),
  AbstractWordSubstitutionModel(gCode->getSourceAlphabet(), new CanonicalStateMap(gCode->getSourceAlphabet(), false), prefix),
  hasParametrizedRates_(paramRates),
  gCode_(gCode)
{
  enableEigenDecomposition(1);

  if ((pmod1 == pmod2) || (pmod2 == pmod3) || (pmod1 == pmod3))
  {
    for (size_t i = 0; i < 3; ++i)
    {
      VSubMod_.push_back(pmod1);
      VnestedPrefix_.push_back(pmod1->getNamespace());
    }

    pmod1->setNamespace(prefix + "123_" + VnestedPrefix_[0]);
    pmod1->enableEigenDecomposition(0);
    addParameters_(pmod1->getParameters());
  }
  else
  {
    VSubMod_.push_back(pmod1);
    VnestedPrefix_.push_back(pmod1->getNamespace());
    VSubMod_[0]->setNamespace(prefix + "1_" + VnestedPrefix_[0]);
    VSubMod_[0]->enableEigenDecomposition(0);
    addParameters_(pmod1->getParameters());

    VSubMod_.push_back(pmod2);
    VnestedPrefix_.push_back(pmod2->getNamespace());
    VSubMod_[1]->setNamespace(prefix + "2_" + VnestedPrefix_[1]);
    VSubMod_[1]->enableEigenDecomposition(0);
    addParameters_(pmod2->getParameters());

    VSubMod_.push_back(pmod3);
    VnestedPrefix_.push_back(pmod3->getNamespace());
    VSubMod_[2]->setNamespace(prefix + "3_" + VnestedPrefix_[2]);
    VSubMod_[2]->enableEigenDecomposition(0);
    addParameters_(pmod3->getParameters());
  }

  Vrate_.resize(3);
  for (size_t i = 0; i < 3; ++i)
  {
    Vrate_[i] = 1.0 / 3;
  }

  if (hasParametrizedRates_)
  {
    // relative rates
    for (int i = 0; i < 2; ++i)
    {
      addParameter_(new Parameter(prefix + "relrate" + TextTools::toString(i + 1), 1.0 / (3.0 - i), &Parameter::PROP_CONSTRAINT_EX));
    }
  }
}

void AbstractCodonSubstitutionModel::updateMatrices()
{
  if (hasParametrizedRates_)
  {
    size_t i, nbmod = VSubMod_.size();
    double x;
    size_t k;
    for (k = nbmod - 1; k > 0; k--)
    {
      x = 1.0;
      for (i = 0; i < k; i++)
      {
        x *= 1 - getParameterValue("relrate" + TextTools::toString(i + 1));
      }
      if (k != nbmod - 1)
        x *= getParameterValue("relrate" + TextTools::toString(k + 1));
      Vrate_[k] = x;
    }
  }

  AbstractWordSubstitutionModel::updateMatrices();
}

void AbstractCodonSubstitutionModel::completeMatrices()
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


