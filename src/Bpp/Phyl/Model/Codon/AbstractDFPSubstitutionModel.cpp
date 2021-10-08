//
// File: AbstractDFPSubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: jeudi 29 octobre 2020, Ã  16h 22
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
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


#include "AbstractDFPSubstitutionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractDFPSubstitutionModel::AbstractDFPSubstitutionModel(
  const GeneticCode* gCode,
  const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractSubstitutionModel(gCode->getSourceAlphabet(), std::shared_ptr<const StateMap>(new CanonicalStateMap(gCode->getSourceAlphabet(), false)), prefix),
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


void AbstractDFPSubstitutionModel::fireParameterChanged(const ParameterList& parameters)
{
  tr_ = getParameterValue("tr");
  trr_ = getParameterValue("trr");
  tvv_ = getParameterValue("tvv");
  trv_ = getParameterValue("trv");
  tsub_ = getParameterValue("tsub");

  updateMatrices();
}


void AbstractDFPSubstitutionModel::updateMatrices()
{
  size_t i, j;

  for (i = 0; i < 64; i++)
  {
    for (j = 0; j < 64; j++)
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
  AbstractSubstitutionModel::updateMatrices();
}


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
