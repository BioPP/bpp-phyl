// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_MULTIPROCESSSEQUENCEEVOLUTION_H
#define BPP_PHYL_LIKELIHOOD_MULTIPROCESSSEQUENCEEVOLUTION_H

#include <Bpp/Numeric/AbstractParametrizable.h>

#include "SequenceEvolution.h"
#include "SubstitutionProcess.h"
#include "SubstitutionProcessCollection.h"

// From SeqLib:
#include <Bpp/Seq/Container/AlignmentData.h>

using namespace std;

namespace bpp
{
/**
 * @brief Partial implementation of multiple processes of sequences.
 *
 */

class MultiProcessSequenceEvolution :
  virtual public SequenceEvolution,
  public AbstractParameterAliasable
{
protected:
  std::shared_ptr<SubstitutionProcessCollection> processColl_;

  /**
   * @brief the vector of the substitution process numbers, as
   * they are used in this order.
   */
  std::vector<size_t> nProc_;

public:
  MultiProcessSequenceEvolution(
      std::shared_ptr<SubstitutionProcessCollection> processColl,
      std::vector<size_t> nProc,
      const std::string& prefix = "");

  MultiProcessSequenceEvolution(const MultiProcessSequenceEvolution& lik) :
    AbstractParameterAliasable(lik),
    processColl_(lik.processColl_),
    nProc_(lik.nProc_)
  {}

  MultiProcessSequenceEvolution& operator=(const MultiProcessSequenceEvolution& lik)
  {
    AbstractParameterAliasable::operator=(lik);

    processColl_ = lik.processColl_;
    nProc_ = lik.nProc_;

    return *this;
  }

public:
  /**
   * @brief The collection
   */
  const SubstitutionProcessCollection& collection() const { return *processColl_; }

  std::shared_ptr<const SubstitutionProcessCollection> getCollection() const { return processColl_; }

  SubstitutionProcessCollection& collection() { return *processColl_; }

  std::shared_ptr<SubstitutionProcessCollection> getCollection() { return processColl_; }

  /**
   * @brief Return the number of process used for computation.
   */
  size_t getNumberOfSubstitutionProcess() const { return nProc_.size(); }

  /**
   * @brief Return the SubstitutionProcess of a given index
   * position (in nProc_ vector).
   */
  const SubstitutionProcessInterface& substitutionProcess(size_t number) const
  {
    return processColl_->substitutionProcess(number);
  }

  std::shared_ptr<const SubstitutionProcessInterface> getSubstitutionProcess(size_t number) const
  {
    return processColl_->getSubstitutionProcess(number);
  }

  const std::vector<size_t>& getSubstitutionProcessNumbers() const
  {
    return nProc_;
  }

  ParameterList getSubstitutionProcessParameters(bool independent) const;

  ParameterList getSubstitutionModelParameters(bool independent) const;

  ParameterList getRateDistributionParameters(bool independent) const;

  ParameterList getRootFrequenciesParameters(bool independent) const;

  ParameterList getBranchLengthParameters(bool independent) const;

  virtual ParameterList getNonDerivableParameters() const;

  virtual void fireParameterChanged(const ParameterList& parameters);

  void setParameters(const ParameterList& parameters);

  /**
   * @brief test if data fits this model
   */
  virtual bool isCompatibleWith(const AlignmentDataInterface& data) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_MULTIPROCESSSEQUENCEEVOLUTION_H
