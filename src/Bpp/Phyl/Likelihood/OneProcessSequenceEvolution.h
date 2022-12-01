//
// File: OneProcessSequenceEvolution.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: mardi 28 avril 2015, ÃÂ  11h 06
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_LIKELIHOOD_ONEPROCESSSEQUENCEEVOLUTION_H
#define BPP_PHYL_LIKELIHOOD_ONEPROCESSSEQUENCEEVOLUTION_H



// From bpp-core:
#include <Bpp/Numeric/ParameterAliasable.h>

// From the STL:
#include <memory>

#include "SequenceEvolution.h"
#include "SubstitutionProcess.h"

namespace bpp
{
/**
 * @brief Evolution of a sequence performed by a unique
 * SubstitutionProcess all along the sequence.
 *
 */

class OneProcessSequenceEvolution :
  public virtual SequenceEvolution,
  public AbstractParameterAliasable
{
protected:
  std::shared_ptr<SubstitutionProcessInterface> subsProc_;

  /**
   * @brief the substitution process number.
   */
  size_t nProc_;

  /**
   * @brief not nice, for inheritance compatibility.
   */
  std::vector<size_t> vProc_;

public:
  OneProcessSequenceEvolution(std::shared_ptr<SubstitutionProcessInterface> process, size_t nProc = 0);

  OneProcessSequenceEvolution(const OneProcessSequenceEvolution& evol);

  OneProcessSequenceEvolution& operator=(const OneProcessSequenceEvolution& evol);

  OneProcessSequenceEvolution* clone() const override
  {
    return new OneProcessSequenceEvolution(*this);
  }

  bool isCompatibleWith(const AlignmentDataInterface& data) const override
  {
    return subsProc_->isCompatibleWith(data);
  }

  const size_t getSubstitutionProcessNumber() const
  {
    return nProc_;
  }

  const SubstitutionProcessInterface& substitutionProcess() const
  {
    return *subsProc_;
  }

  std::shared_ptr<const SubstitutionProcessInterface> getSubstitutionProcess() const
  {
    return subsProc_;
  }

  SubstitutionProcessInterface& substitutionProcess()
  {
    return *subsProc_;
  }

  std::shared_ptr<SubstitutionProcessInterface> getSubstitutionProcess()
  {
    return subsProc_;
  }

  const std::vector<size_t>& getSubstitutionProcessNumbers() const override
  {
    return vProc_;
  }

  const SubstitutionProcessInterface& substitutionProcess(size_t number) const override
  {
    return *subsProc_;
  }

  std::shared_ptr<const SubstitutionProcessInterface> getSubstitutionProcess(size_t number) const override
  {
    return subsProc_;
  }

  /**
   * @brief Get the tree (topology and branch lengths).
   *
   * @return The tree of this OneProcessSequenceEvolution object.
   */

  const ParametrizablePhyloTree& tree() const
  {
    return subsProc_->parametrizablePhyloTree();
  }

  std::shared_ptr<const ParametrizablePhyloTree> getTree() const
  {
    return subsProc_->getParametrizablePhyloTree();
  }

  void fireParameterChanged(const ParameterList& pl) override
  {
    // Updates substitution process:
    subsProc_->matchParametersValues(pl);
  }


  /**
   * @brief Get several categories of parameters
   *
   **/
  ParameterList getBranchLengthParameters(bool independent) const override
  {
    return subsProc_->getBranchLengthParameters(independent);
  }

  ParameterList getRootFrequenciesParameters(bool independent) const override
  {
    return subsProc_->getRootFrequenciesParameters(independent);
  }

  ParameterList getRateDistributionParameters(bool independent) const override
  {
    return subsProc_->getRateDistributionParameters(independent);
  }

  ParameterList getSubstitutionModelParameters(bool independent) const override
  {
    return subsProc_->getSubstitutionModelParameters(independent);
  }

  ParameterList getSubstitutionProcessParameters(bool independent) const override
  {
    if (independent)
      return subsProc_->getIndependentParameters();
    else
      return subsProc_->getParameters();
  }

  const ParameterList& getParameters() const override
  {
    return subsProc_->getParameters();
  }

  const ParameterList& getIndependentParameters() const override
  {
    return subsProc_->getIndependentParameters();
  }

  /**
   * @brief get (Non)Derivable INDEPENDENT parameters
   *
   **/
  ParameterList getNonDerivableParameters() const override
  {
    return subsProc_->getNonDerivableParameters();
  }
};
} // end namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_ONEPROCESSSEQUENCEEVOLUTION_H
