// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_HMMPROCESSALPHABET_H
#define BPP_PHYL_LIKELIHOOD_HMMPROCESSALPHABET_H


// From Numeric
#include <Bpp/Numeric/Hmm/HmmStateAlphabet.h>

#include "SubstitutionProcessCollection.h"


namespace bpp
{
/**
 * @brief Hidden states alphabet.
 *
 * Implementation of HmmStateAlphabet where Alphabet States are
 * Substitution Process belonging to a collection.
 */
class HmmProcessAlphabet :
  public virtual HmmStateAlphabet,
  public AbstractParametrizable
{
private:
  std::shared_ptr<const SubstitutionProcessCollection> processColl_;

  /**
   * @brief the vector of the substitution process numbers.
   */
  std::vector<size_t> nProc_;

public:
  HmmProcessAlphabet(
      std::shared_ptr<const SubstitutionProcessCollection> pSub,
      std::vector<size_t> nProc) :
    AbstractParametrizable(""),
    processColl_(pSub),
    nProc_(nProc)
  {}

  HmmProcessAlphabet(const HmmProcessAlphabet& hpa) :
    AbstractParametrizable(hpa),
    processColl_(hpa.processColl_),
    nProc_(hpa.nProc_)
  {}

  HmmProcessAlphabet& operator=(const HmmProcessAlphabet& hpa)
  {
    AbstractParametrizable::operator=(*this);
    processColl_ = hpa.processColl_;
    nProc_ = hpa.nProc_;

    return *this;
  }

  HmmProcessAlphabet* clone() const override {return new HmmProcessAlphabet(*this);}

  virtual ~HmmProcessAlphabet() {}

  /**
   * @param stateIndex The index of a hidden state.
   * @return The corresponding hidden state.
   * @see getNumberOfStates
   */
  const Clonable& getState(size_t stateIndex) const override
  {
    return processColl_->substitutionProcess(nProc_[stateIndex]);
  }

  size_t getNumberOfStates() const override
  {
    return nProc_.size();
  }

  /**
   * @brief Tell if this instance can work with the instance of alphabet given as input.
   *
   * In many case, this will return true if the pointer provided as argument refers to this object.
   *
   * @param stateAlphabet The alphabet to check.
   * @return true if the matrix is compatible with the given alphabet.
   */
  bool worksWith(const HmmStateAlphabet& stateAlphabet) const override
  {
    return &stateAlphabet == this;
  }
};
}
#endif // BPP_PHYL_LIKELIHOOD_HMMPROCESSALPHABET_H
