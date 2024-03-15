// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "AbstractDFPSubstitutionModel.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

AbstractDFPSubstitutionModel::AbstractDFPSubstitutionModel(
    shared_ptr<const GeneticCode> gCode,
    const string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractSubstitutionModel(
      gCode->getSourceAlphabet(),
      shared_ptr<const StateMapInterface>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)),
      prefix),
  gCode_(gCode),
  tr_(1), trr_(1), tvv_(1), trv_(1), tsub_(1)
{
  enableEigenDecomposition(true);
  addParameter_(new Parameter(prefix + "tr", 1, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter(prefix + "trr", 1, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter(prefix + "tvv", 1, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter(prefix + "trv", 1, Parameter::R_PLUS_STAR));
  addParameter_(new Parameter(prefix + "tsub", 1, Parameter::R_PLUS_STAR));
}

/******************************************************************************/

void AbstractDFPSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  tr_ = getParameterValue("tr");
  trr_ = getParameterValue("trr");
  tvv_ = getParameterValue("tvv");
  trv_ = getParameterValue("trv");
  tsub_ = getParameterValue("tsub");

  updateMatrices_();
}

/******************************************************************************/

void AbstractDFPSubstitutionModel::updateMatrices_()
{
  size_t i, j;

  for (i = 0; i < 64; ++i)
  {
    for (j = 0; j < 64; ++j)
    {
      if (i == j || gCode_->isStop(static_cast<int>(i)) || gCode_->isStop(static_cast<int>(j)))
      {
        generator_(i, j) = 0;
      }
      else
      {
        generator_(i, j) = getCodonsMulRate(i, j);
      }
    }
  }

  setDiagonal();
  AbstractSubstitutionModel::updateMatrices_();
}

/******************************************************************************/

double AbstractDFPSubstitutionModel::getCodonsMulRate(size_t i, size_t j) const
{
  int nts(0);
  // sum nb of transitions & transversions between codons
  // ts = 4 * ntv + ntr

  for (size_t pos = 0; pos < 3; pos++)
  {
    int pi = (int) (pos == 0 ? i / 16 :
                    (pos == 1 ? (i / 4) % 4
                     : i % 4));
    int pj = (int) (pos == 0 ? j / 16 :
                    (pos == 1 ? (j / 4) % 4
                     : j % 4));

    nts += (pi == pj ? 0 : (abs(pi - pj) == 2 ? 1 : 4));
  }

  switch (nts)
  {
  case 0:
    break;
  case 1:
    return tr_;
    break;
  case 2:
    return trr_;
    break;
  case 4:
    return 1;
    break;
  case 5:
    return trv_;
    break;
  case 8:
    return tvv_;
    break;
  default:
    return tsub_;
    break;
  }

  return 1.;
}

/******************************************************************************/

