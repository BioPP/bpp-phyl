//
// File: LGL08_CAT.h
// Created by: Mathieu Groussin
// Created on: Tuesday 11 December 2012
//

/*
  Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _LGL08_CAT_H_
#define _LGL08_CAT_H_

#include "../MixtureOfSubstitutionModels.h"
#include "ProteinSubstitutionModel.h"
#include "../AbstractSubstitutionModel.h"
#include "../AbstractBiblioMixedTransitionModel.h"

namespace bpp
{
/**
 * @brief The Le et al  (2008) CAT substitution model for proteins.
 * @author Mathieu Groussin
 *
 * This model is a mixture of N profiles empirically built with an EM algorithm
 * (see ref). The submodels are called C1, C2, ..., CN. For each model, exchangeabilities are equal (F81 model).
 *
 *
 * This model includes 2N-2 parameters :
 *
 * - relrate1, ..., relrate(N-1) are the relative rates of model C1, ..., C(N-1);
 * - relproba1, ..., relproba(N-1) are the relative proportions of model C1, ..., C(N-1);
 *
 * Important: See the relation between these parameters and the
 * rates and probabilities of the models in the description of
 * MixtureOfSubstitutionModels class.
 *
 * Reference:
 *
 * Le S.Q., Gascuel O. & Lartillot N. (2008) Bioinformatics 24:2317–2323.
 */

  class LGL08_CAT :
    public AbstractBiblioMixedTransitionModel
  {
  public:
    class EmbeddedModel :
      public AbstractReversibleProteinSubstitutionModel
    {
    private:
      double proportion_;
      std::string name_;

    public:
      EmbeddedModel(const ProteicAlphabet* alpha, std::string name, unsigned int nbCat = 10);
      virtual ~EmbeddedModel(){}
      EmbeddedModel* clone() const { return new EmbeddedModel(*this); }
      std::string getName() const { return name_;}
      double getProportion() const { return proportion_;}
    };

  public:
    /**
     * @brief Build a CAT model, with original equilibrium frequencies, probabilities and rates.
     *
     * @param alpha A proteic alphabet.
     * @param nbCat number of profiles
     *
     */
    LGL08_CAT(const ProteicAlphabet* alpha, unsigned int nbCat = 10);

    LGL08_CAT* clone() const { return new LGL08_CAT(*this); }

    LGL08_CAT(const LGL08_CAT& mod2) :
      AbstractBiblioMixedTransitionModel(mod2)
    {}

    LGL08_CAT& operator=(const LGL08_CAT& mod2)
    {
      AbstractBiblioMixedTransitionModel::operator=(mod2);
      return *this;
    }

    uint getNumberOfCategories() const
    {
      return (uint)pmixmodel_->getNumberOfModels();
    }
    
    std::string getName() const { return "LGL08_CAT";}

  };
} // end of namespace bpp.

#endif  // _LGL08_CAT_H_

