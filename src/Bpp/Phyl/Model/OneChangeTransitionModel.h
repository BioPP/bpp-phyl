//
// File: OneChangeTransitionModel.h
// Authors:
//   Laurent Gueguen
// Created: samedi 24 octobre 2015, Ã  18h 28
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

#ifndef BPP_PHYL_MODEL_ONECHANGETRANSITIONMODEL_H
#define BPP_PHYL_MODEL_ONECHANGETRANSITIONMODEL_H


#include "AbstractFromSubstitutionModelTransitionModel.h"

namespace bpp
{
/**
 * @brief From a model, compute transition probabilities given there
 * is at least a change in the branch.
 *
 * It has the same parameters as the SubModel.
 */

class OneChangeTransitionModel :
  public AbstractFromSubstitutionModelTransitionModel
{
public:
  OneChangeTransitionModel(std::unique_ptr<SubstitutionModelInterface> originalModel) :
    AbstractWrappedModel("Onechange."),
    AbstractWrappedTransitionModel("Onechange."),
    AbstractFromSubstitutionModelTransitionModel(std::move(originalModel), "OneChange.")
  {}

  OneChangeTransitionModel(const OneChangeTransitionModel& fmsm) :
    AbstractWrappedModel(fmsm),
    AbstractWrappedTransitionModel(fmsm),
    AbstractFromSubstitutionModelTransitionModel(fmsm)
  {}

  OneChangeTransitionModel& operator=(const OneChangeTransitionModel& fmsm)
  {
    AbstractWrappedModel::operator=(fmsm);
    AbstractWrappedTransitionModel::operator=(fmsm);
    AbstractFromSubstitutionModelTransitionModel::operator=(fmsm);
    return *this;
  }

  virtual ~OneChangeTransitionModel() {}

  OneChangeTransitionModel* clone() const override { return new OneChangeTransitionModel(*this); }

public:
  double Pij_t    (size_t i, size_t j, double t) const override;
  double dPij_dt  (size_t i, size_t j, double t) const override;
  double d2Pij_dt2(size_t i, size_t j, double t) const override;

  const Matrix<double>& getPij_t(double t) const override;

  const Matrix<double>& getdPij_dt(double t) const override;

  const Matrix<double>& getd2Pij_dt2(double t) const override;

  double freq(size_t i) const override { return transitionModel().freq(i); }

  const Vdouble& getFrequencies() const override { return transitionModel().getFrequencies(); }

  const FrequencySetInterface& frequencySet() const override
  {
    return transitionModel().frequencySet();
  }

  std::shared_ptr<const FrequencySetInterface> getFrequencySet() const override
  {
    return transitionModel().getFrequencySet();
  }

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount) override
  {
    transitionModel().setFreqFromData(data, pseudoCount);
  }

  virtual void setFreq(std::map<int, double>& m) override
  {
    transitionModel().setFreq(m);
  }

  double getRate() const override { return transitionModel().getRate(); }

  void setRate(double rate) override { return transitionModel().setRate(rate); }

  double getInitValue(size_t i, int state) const override { return model().getInitValue(i, state); }

  std::string getName() const override
  {
    return "OneChange";
  }

  /**
   * @}
   */
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_ONECHANGETRANSITIONMODEL_H
