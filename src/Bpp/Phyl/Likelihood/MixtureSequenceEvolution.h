// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_MIXTURESEQUENCEEVOLUTION_H
#define BPP_PHYL_LIKELIHOOD_MIXTURESEQUENCEEVOLUTION_H

#include <Bpp/Numeric/Prob/Simplex.h>

#include "MultiProcessSequenceEvolution.h"

namespace bpp
{
/**
 * @brief Sequence evolution framework based on a mixture of
 * substitution processes
 *
 * @see MultiProcessSequencePhyloLikelihood
 */

class MixtureSequenceEvolution :
  public MultiProcessSequenceEvolution
{
private:
  Simplex simplex_;

public:
  MixtureSequenceEvolution(
    std::shared_ptr<SubstitutionProcessCollection> processColl,
    std::vector<size_t>& nProc);

  MixtureSequenceEvolution(const MixtureSequenceEvolution& mlc) :
    MultiProcessSequenceEvolution(mlc),
    simplex_(mlc.simplex_) {}

  MixtureSequenceEvolution& operator=(const MixtureSequenceEvolution& mlc)
  {
    MultiProcessSequenceEvolution::operator=(mlc);
    simplex_ = mlc.simplex_;
    return *this;
  }

  virtual ~MixtureSequenceEvolution() {}

  MixtureSequenceEvolution* clone() const override { return new MixtureSequenceEvolution(*this); }

public:
  void setNamespace(const std::string& nameSpace) override;

  void fireParameterChanged(const ParameterList& parameters) override;

  ParameterList getNonDerivableParameters() const override;

  const std::vector<double>& getSubProcessProbabilities() const
  {
    return simplex_.getFrequencies();
  }

  Simplex& simplex()
  {
    return simplex_;
  }

  /**
   * @brief return the probability of a  subprocess
   *
   * @param i the index of the subprocess
   */
  double getSubProcessProb(size_t i) const
  {
    return simplex_.prob(i);
  }

  /**
   * @brief Set the probabilities of the subprocess
   */
  void setSubProcessProb(const Simplex& si);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_MIXTURESEQUENCEEVOLUTION_H
