// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_HMMSEQUENCEEVOLUTION_H
#define BPP_PHYL_LIKELIHOOD_HMMSEQUENCEEVOLUTION_H


#include "HmmProcessAlphabet.h"
#include "MultiProcessSequenceEvolution.h"

// From Numeric
#include <Bpp/Numeric/Hmm/FullHmmTransitionMatrix.h>

#include <memory>

namespace bpp
{
/**
 * @brief Sequence evolution framework based on a hmm
 */
class HmmSequenceEvolution :
  public MultiProcessSequenceEvolution
{
private:
  std::shared_ptr<HmmProcessAlphabet> hmmAlph_;
  std::shared_ptr<FullHmmTransitionMatrix> hmmTransMat_;

public:
  HmmSequenceEvolution(
      std::shared_ptr<SubstitutionProcessCollection> processColl,
      std::vector<size_t>& nProc);

  HmmSequenceEvolution(const HmmSequenceEvolution& mlc) :
    MultiProcessSequenceEvolution(mlc),
    hmmAlph_(mlc.hmmAlph_->clone()),
    hmmTransMat_(mlc.hmmTransMat_->clone()) {}

  HmmSequenceEvolution& operator=(const HmmSequenceEvolution& mlc)
  {
    MultiProcessSequenceEvolution::operator=(mlc);
    hmmAlph_ = std::make_shared<HmmProcessAlphabet>(*mlc.hmmAlph_);
    hmmTransMat_ = std::make_shared<FullHmmTransitionMatrix>(*mlc.hmmTransMat_);

    return *this;
  }

  virtual ~HmmSequenceEvolution() {}

  HmmSequenceEvolution* clone() const { return new HmmSequenceEvolution(*this); }

public:
  void setNamespace(const std::string& nameSpace);

  void fireParameterChanged(const ParameterList& parameters);

  const FullHmmTransitionMatrix& hmmTransitionMatrix() const
  {
    return *hmmTransMat_;
  }

  FullHmmTransitionMatrix& hmmTransitionMatrix()
  {
    return *hmmTransMat_;
  }

  std::shared_ptr<FullHmmTransitionMatrix> getHmmTransitionMatrix()
  {
    return hmmTransMat_;
  }

  const HmmProcessAlphabet& hmmProcessAlphabet() const
  {
    return *hmmAlph_;
  }

  HmmProcessAlphabet& hmmProcessAlphabet()
  {
    return *hmmAlph_;
  }

  std::shared_ptr<HmmProcessAlphabet> getHmmProcessAlphabet()
  {
    return hmmAlph_;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_HMMSEQUENCEEVOLUTION_H
