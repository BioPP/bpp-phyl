//
// File: AbstractBiblioMixedTransitionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: lundi 18 juillet 2011, Ã  15h 27
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


#include "AbstractBiblioMixedTransitionModel.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

AbstractBiblioMixedTransitionModel::AbstractBiblioMixedTransitionModel(const std::string& prefix) :
  AbstractBiblioTransitionModel(prefix),
  pmixmodel_()
{}

AbstractBiblioMixedTransitionModel::AbstractBiblioMixedTransitionModel(const AbstractBiblioMixedTransitionModel& mod2) :
  AbstractBiblioTransitionModel(mod2),
  pmixmodel_(mod2.pmixmodel_->clone())
{}

AbstractBiblioMixedTransitionModel& AbstractBiblioMixedTransitionModel::operator=(const AbstractBiblioMixedTransitionModel& mod2)
{
  AbstractBiblioTransitionModel::operator=(mod2);
  pmixmodel_.reset(mod2.pmixmodel_->clone());
  return *this;
}

AbstractBiblioMixedTransitionModel::~AbstractBiblioMixedTransitionModel()
{}

Vuint AbstractBiblioMixedTransitionModel::getSubmodelNumbers(const std::string& desc) const
{
  std::string desc2;

  StringTokenizer st(desc, ",");
  while (st.hasMoreToken())
  {
    string param = st.nextToken();

    desc2 += getParameterNameWithoutNamespace(param);
    if (st.hasMoreToken())
      desc2 += ",";
  }

  return getMixedModel().getSubmodelNumbers(desc2);
}
