//
// File: InMixedSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: vendredi 22 septembre 2017, Ã  09h 57
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

#ifndef BPP_PHYL_MODEL_INMIXEDSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_INMIXEDSUBSTITUTIONMODEL_H


#include "AbstractSubstitutionModel.h"
#include "AbstractWrappedModel.h"
#include "MixedTransitionModel.h"

namespace bpp
{
/**
 * @brief SubModel taken from a MixedTransitionModel, kept in the
 * context of the MixedTransitionModel (see
 * FromMixtureSubstitutionModel for an out of context subModel). So
 * "rate" and "scale" are set for the MixedTransitionModel.
 *
 * But method getRate returns the specific rate of the subModel.
 *
 * It owns the MixedTransitionModel.
 *
 * It has the same parameters as the MixedTransitionModel.
 */

class InMixedSubstitutionModel :
  public virtual AbstractWrappedSubstitutionModel
{
private:
  /**
   * @brief The MixedOfTransitionModels.
   */
  std::unique_ptr<MixedTransitionModelInterface> mixedModelPtr_;

  /**
   * @brief the number of the submodel
   */
  size_t subModelNumber_;

  /**
   * @brief The name of the mixture model (for io purpose).
   */
  std::string mixtName_;

public:
  InMixedSubstitutionModel(
      std::unique_ptr<MixedTransitionModelInterface> mixedModel,
      const std::string& subModelName,
      const std::string& mixtDesc);

  InMixedSubstitutionModel(
      std::unique_ptr<MixedTransitionModelInterface> mixedModel,
      size_t subModelNumber,
      const std::string& mixtDesc);

  InMixedSubstitutionModel(const InMixedSubstitutionModel& fmsm);

  InMixedSubstitutionModel& operator=(const InMixedSubstitutionModel& fmsm);

  InMixedSubstitutionModel* clone() const override { return new InMixedSubstitutionModel(*this); }

public:
  const MixedTransitionModelInterface& mixedModel() const
  {
    return *mixedModelPtr_;
  }

  const SubstitutionModelInterface& substitutionModel() const override
  {
    return dynamic_cast<const SubstitutionModelInterface&>(mixedModelPtr_->nModel(subModelNumber_));
  }

  size_t getSubModelNumber() const
  {
    return subModelNumber_;
  }

  bool computeFrequencies() const override
  {
    return mixedModel().computeFrequencies();
  }

  /**
   * @return Set if equilibrium frequencies should be computed from
   * the generator
   */
  void computeFrequencies(bool yn) override
  {
    mixedModel_().computeFrequencies(yn);
  }

protected:
  Vdouble& getFrequencies_() override
  {
    return mixedModel_().getFrequencies_();
  }

  /**
   * @}
   */
  MixedTransitionModelInterface& mixedModel_()
  {
    return *mixedModelPtr_;
  }

  SubstitutionModelInterface& substitutionModel_() override
  {
    return dynamic_cast<SubstitutionModelInterface&>(mixedModel_().nModel_(subModelNumber_));
  }

public:
  /**
   * @ brief Methods to supersede WrappedSubstitutionModel methods.
   *
   * @{
   */
  double freq(size_t i) const override { return transitionModel().freq(i); }

  double Pij_t    (size_t i, size_t j, double t) const override { return transitionModel().Pij_t(i, j, t); }
  double dPij_dt  (size_t i, size_t j, double t) const override { return transitionModel().dPij_dt (i, j, t); }
  double d2Pij_dt2(size_t i, size_t j, double t) const override { return transitionModel().d2Pij_dt2(i, j, t); }

  const Vdouble& getFrequencies() const override { return transitionModel().getFrequencies(); }

  const Matrix<double>& getPij_t(double t) const override { return transitionModel().getPij_t(t); }

  const Matrix<double>& getdPij_dt(double t) const override { return transitionModel().getdPij_dt(t); }

  const Matrix<double>& getd2Pij_dt2(double t) const override { return transitionModel().getd2Pij_dt2(t); }

  double getInitValue(size_t i, int state) const override
  {
    return model().getInitValue(i, state);
  }

  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override
  {
    mixedModel_().setFreqFromData(data, pseudoCount);
  }

  void setFreq(std::map<int, double>& frequencies) override
  {
    mixedModel_().setFreq(frequencies);
  }

  /**
   * @brief Methods to supersede SubstitutionModel methods.
   *
   * @{
   */
  double Qij(size_t i, size_t j) const override { return substitutionModel().Qij(i, j); }

  const Matrix<double>& getGenerator() const override { return substitutionModel().getGenerator(); }

  const Matrix<double>& getExchangeabilityMatrix() const override { return substitutionModel().getExchangeabilityMatrix(); }

  double Sij(size_t i, size_t j) const override { return substitutionModel().Sij(i, j); }

  void enableEigenDecomposition(bool yn) override { substitutionModel_().enableEigenDecomposition(yn); }

  bool enableEigenDecomposition() override { return substitutionModel_().enableEigenDecomposition(); }

  bool isDiagonalizable() const override { return substitutionModel().isDiagonalizable(); }

  bool isNonSingular() const override { return substitutionModel().isNonSingular(); }

  const Vdouble& getEigenValues() const override { return substitutionModel().getEigenValues(); }

  const Vdouble& getIEigenValues() const override { return substitutionModel().getIEigenValues(); }

  const Matrix<double>& getRowLeftEigenVectors() const override { return substitutionModel().getRowLeftEigenVectors(); }

  const Matrix<double>& getColumnRightEigenVectors() const override { return substitutionModel().getColumnRightEigenVectors(); }


  /**
   * @}
   */
  bool isScalable() const override
  {
    return substitutionModel().isScalable();
  }

  void setScalable(bool scalable) override 
  {
    substitutionModel_().setScalable(scalable);
  }


  double getScale() const override { return substitutionModel().getScale(); }

  void setScale(double scale) override { substitutionModel_().setScale(scale); }


  void normalize() override
  {
    substitutionModel_().normalize();
  }

  void setDiagonal() override
  {
    substitutionModel_().setDiagonal();
  }

  double getRate() const override
  {
    return transitionModel().getRate();
  }

  void setRate(double rate) override
  {
    return mixedModel_().setRate(rate);
  }

  void addRateParameter() override
  {
    mixedModel_().addRateParameter();
    addParameter_(new Parameter(getNamespace() + "rate", mixedModel().getRate(), Parameter::R_PLUS_STAR));
  }

  /**
   * @}
   */

  /**
   * @brief Methods to supersede AbstractSubstitutionnModel methods.
   *
   * @{
   */

  /**
   * @brief Tells the model that a parameter value has changed.
   *
   * This updates the matrices consequently.
   */
  void fireParameterChanged(const ParameterList& parameters) override
  {
    mixedModel_().matchParametersValues(parameters);
  }

  void setNamespace(const std::string& name) override
  {
    AbstractParameterAliasable::setNamespace(name);
    mixedModel_().setNamespace(name);
  }


  /**
   * @}
   */
  std::string getName() const override
  {
    return mixedModelPtr_->getName();
  }

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_INMIXEDSUBSTITUTIONMODEL_H
