// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_AUTOCORRELATIONSEQUENCEEVOLUTION_H
#define BPP_PHYL_LIKELIHOOD_AUTOCORRELATIONSEQUENCEEVOLUTION_H


#include "HmmProcessAlphabet.h"
#include "MultiProcessSequenceEvolution.h"

// From Numeric
#include <Bpp/Numeric/Hmm/AutoCorrelationTransitionMatrix.h>

#include <memory>

namespace bpp
{
/**
 * @brief Sequence evolution framework based on an auto-correlation of
 * substitution processes.
 *
 */
class AutoCorrelationSequenceEvolution :
  public MultiProcessSequenceEvolution
{
private:
  std::shared_ptr<HmmProcessAlphabet> hmmAlph_;
  std::shared_ptr<AutoCorrelationTransitionMatrix> autoCorrTransMat_;

public:
  AutoCorrelationSequenceEvolution(
    std::shared_ptr<SubstitutionProcessCollection> processColl,
    std::vector<size_t>& nProc);

  AutoCorrelationSequenceEvolution(const AutoCorrelationSequenceEvolution& mlc) :
    MultiProcessSequenceEvolution(mlc),
    hmmAlph_(mlc.hmmAlph_->clone()),
    autoCorrTransMat_(mlc.autoCorrTransMat_->clone()){}

  AutoCorrelationSequenceEvolution& operator=(const AutoCorrelationSequenceEvolution& mlc)
  {
    MultiProcessSequenceEvolution::operator=(mlc);

    hmmAlph_ = std::shared_ptr<HmmProcessAlphabet>(new HmmProcessAlphabet(*mlc.hmmAlph_.get()));
    autoCorrTransMat_ = std::shared_ptr<AutoCorrelationTransitionMatrix>(new AutoCorrelationTransitionMatrix(*mlc.autoCorrTransMat_.get()));

    return *this;
  }

  virtual ~AutoCorrelationSequenceEvolution() {}

  AutoCorrelationSequenceEvolution* clone() const { return new AutoCorrelationSequenceEvolution(*this); }

public:
  void setNamespace(const std::string& nameSpace);

  void fireParameterChanged(const ParameterList& parameters);

  const AutoCorrelationTransitionMatrix& hmmTransitionMatrix() const
  {
    return *autoCorrTransMat_.get();
  }

  AutoCorrelationTransitionMatrix& hmmTransitionMatrix()
  {
    return *autoCorrTransMat_.get();
  }

  std::shared_ptr<AutoCorrelationTransitionMatrix> getHmmTransitionMatrix()
  {
    return autoCorrTransMat_;
  }

  std::shared_ptr<const AutoCorrelationTransitionMatrix> getHmmTransitionMatrix() const
  {
    return autoCorrTransMat_;
  }

  HmmProcessAlphabet& hmmProcessAlphabet()
  {
    return *hmmAlph_.get();
  }

  const HmmProcessAlphabet& hmmProcessAlphabet() const
  {
    return *hmmAlph_.get();
  }

  std::shared_ptr<HmmProcessAlphabet> getHmmProcessAlphabet()
  {
    return hmmAlph_;
  }

  std::shared_ptr<const HmmProcessAlphabet> getHmmProcessAlphabet() const
  {
    return hmmAlph_;
  }

  ParameterList getNonDerivableParameters() const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_AUTOCORRELATIONSEQUENCEEVOLUTION_H
