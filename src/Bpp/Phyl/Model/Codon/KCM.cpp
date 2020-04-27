//
// File: KCM.cpp
// Created by:  Laurent Gueguen
// Created on: vendredi 23 septembre 2016, à 11h 39
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

#include "KCM.h"
#include "../Nucleotide/GTR.h"

#include "../FrequencySet/CodonFrequencySet.h"
#include <Bpp/Numeric/NumConstants.h>

using namespace bpp;

using namespace std;

/******************************************************************************/

KCM::KCM(const GeneticCode* gc, bool oneModel) :
  AbstractBiblioTransitionModel("KCM"+string(oneModel?"7":"19")+"."),
  AbstractBiblioSubstitutionModel("KCM"+string(oneModel?"7":"19")+"."),
  pmodel_(),
  oneModel_(oneModel)
{
  const NucleicAlphabet* nalph=dynamic_cast<const CodonAlphabet*>(gc->getSourceAlphabet())->getNucleicAlphabet();
  
  if (oneModel)
    pmodel_.reset(new KroneckerCodonDistanceSubstitutionModel(gc, new GTR(nalph)));
  else
    pmodel_.reset(new KroneckerCodonDistanceSubstitutionModel(gc, new GTR(nalph), new GTR(nalph), new GTR(nalph)));

  string name="KCM"+string(oneModel?"7":"19")+".";

  pmodel_->setNamespace(name);

  addParameters_(pmodel_->getParameters());

  getParameter_("beta").setName(name+"omega"),
  
  lParPmodel_.addParameters(pmodel_->getParameters());

  vector<std::string> v = lParPmodel_.getParameterNames();

  for (size_t i = 0; i < v.size(); i++)
    mapParNamesFromPmodel_[v[i]] = getParameterNameWithoutNamespace(v[i]);

  mapParNamesFromPmodel_[name+"beta"] = "omega";

  updateMatrices();
}


KCM::KCM(const KCM& kcm) :
  AbstractBiblioTransitionModel(kcm),
  AbstractBiblioSubstitutionModel(kcm),
  pmodel_(new KroneckerCodonDistanceSubstitutionModel(*kcm.pmodel_)),
  oneModel_(kcm.oneModel_)
{}

KCM& KCM::operator=(const KCM& kcm)
{
  AbstractBiblioSubstitutionModel::operator=(kcm);

  oneModel_=kcm.oneModel_;
  pmodel_.reset(new KroneckerCodonDistanceSubstitutionModel(*kcm.pmodel_));
  return *this;
}

