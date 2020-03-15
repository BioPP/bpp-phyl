//
// File: MixtureOfASubstitutionModel.h
// Created by: Laurent Gueguen
// Date: lundi 13 septembre 2010, à 21h 31
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

#ifndef _MIXTURE_OF_A_SUBSTITUTION_MODEL_H_
#define _MIXTURE_OF_A_SUBSTITUTION_MODEL_H_

#include <Bpp/Numeric/VectorTools.h>
#include "MixtureOfATransitionModel.h"

#include <vector>
#include <string>
#include <map>
#include <cstring> // C lib for string copy

namespace bpp
{

class MixtureOfASubstitutionModel :
  public MixtureOfATransitionModel
{
public:
  /**
   * @brief Constructor of a MixtureOfASubstitutionModel, where all
   * the models have rate 1 and equal probability.
   *
   * @param alpha pointer to the Alphabet
   * @param vpModel vector of pointers to ASubstitutionModel. All the
   *   ASubstitutionModel are owned by the instance.
   * @warning providing a vpModel with size 0 will generate a segmentation fault!
   */
  MixtureOfASubstitutionModel(const Alphabet* alpha,
                              SubstitutionModel* model,
                              std::map<std::string, DiscreteDistribution*> parametersDistributionsList,
                              int ffrom = -1,
                              int tto = -1) :
    AbstractParameterAliasable(model->getNamespace()),
    AbstractTransitionModel(alpha, model->shareStateMap(), model->getNamespace()),
    MixtureOfATransitionModel(alpha, model, parametersDistributionsList, ffrom, tto)
  {
  }

  MixtureOfASubstitutionModel(const MixtureOfASubstitutionModel& model) :
    AbstractParameterAliasable(model),
    AbstractTransitionModel(model),
    MixtureOfATransitionModel(model)
  {
  }
  

  MixtureOfASubstitutionModel& operator=(const MixtureOfASubstitutionModel& model)
  {
    MixtureOfATransitionModel::operator=(model);
    return *this;
  }
  
  MixtureOfASubstitutionModel* clone() const { return new MixtureOfASubstitutionModel(*this); }

  void updateMatrices()
  {
    MixtureOfATransitionModel::updateMatrices();
    // setting the rates, if to_ & from_ are different from -1

    if (to_ >= 0 && from_ >= 0)
    {
      Vdouble vd;

      for (size_t j = 0; j < modelsContainer_.size(); j++)
      {
        vd.push_back(1 / getSubNModel(j)->Qij(static_cast<size_t>(from_), static_cast<size_t>(to_)));
      }

      setVRates(vd);
    }
  }
  
  /**
   * @brief retrieve a pointer to the subsitution model with the given name.
   *
   * Return Null if not found.
   *
   */
  
  const SubstitutionModel* getSubModel(const std::string& name) const
  {
    return dynamic_cast<const SubstitutionModel*>(getModel(name));
  }

  const SubstitutionModel* getSubNModel(size_t i) const
  {
    return dynamic_cast<const SubstitutionModel*>(getNModel(i));
  }

};
} // end of namespace bpp.

#endif  // _MIXTUREOFASUBSTITUTIONMODEL_H_

