// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_PHYLOMAPPINGS_ONEPROCESSSEQUENCESUBSTITUTIONMAPPING_H
#define BPP_PHYL_MAPPING_PHYLOMAPPINGS_ONEPROCESSSEQUENCESUBSTITUTIONMAPPING_H


#include "../../Likelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.h"
#include "../SubstitutionMappingTools.h"

#include "AbstractSinglePhyloSubstitutionMapping.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

namespace bpp
{
/**
 * @brief The OneProcessSequenceSubstitutionMapping class: substitution
 * mapping linked with a OneProcessSequencePhyloLikelihood
 *
 */

class OneProcessSequenceSubstitutionMapping :
  public AbstractSinglePhyloSubstitutionMapping,
  public std::enable_shared_from_this<OneProcessSequenceSubstitutionMapping>
{
private:
  std::shared_ptr<OneProcessSequencePhyloLikelihood> pOPSP_;

  /**
   * @brief Set the models of the BranchedModelSet to the adhoc
   * branches, for normalization.
   */
  void setBranchedModelSet_();

public:
  OneProcessSequenceSubstitutionMapping(
      std::shared_ptr<OneProcessSequencePhyloLikelihood> spp,
      std::shared_ptr<SubstitutionRegisterInterface> reg,
      std::shared_ptr<const AlphabetIndex2> weights,
      std::shared_ptr<const AlphabetIndex2> distances);

  OneProcessSequenceSubstitutionMapping(const OneProcessSequenceSubstitutionMapping& sppm) :
    AbstractSinglePhyloSubstitutionMapping(sppm),
    pOPSP_(sppm.pOPSP_)
  {}

  OneProcessSequenceSubstitutionMapping& operator=(const OneProcessSequenceSubstitutionMapping& sppm)
  {
    AbstractSinglePhyloSubstitutionMapping::operator=(sppm);
    pOPSP_ = sppm.pOPSP_;

    return *this;
  }

  virtual ~OneProcessSequenceSubstitutionMapping() {}

  OneProcessSequenceSubstitutionMapping* clone() const override
  {
    return new OneProcessSequenceSubstitutionMapping(*this);
  }

  /**
   * @brief ComputeCounts
   */
  void computeCounts(short unresolvedOption = SubstitutionMappingTools::UNRESOLVED_ZERO,
      double threshold = -1, bool verbose = true) override
  {
    counts_ = SubstitutionMappingTools::computeCounts(
          getLikelihoodCalculationSingleProcess(),
          getSubstitutionRegister(),
          getWeights(),
          getDistances(),
          unresolvedOption,
          threshold,
          verbose);
  }

  /*
   * @param nullParams parameters values used for normalization
   * @param unresolvedOption  mgmt of gaps in the counts (default: counted as zeros)
   * @param verbose
   */
  
  void computeNormalizations(const ParameterList& nullParams,
      short unresolvedOption = SubstitutionMappingTools::UNRESOLVED_ZERO,
      bool verbose = true) override;


  size_t getNumberOfModels() const override
  {
    return pOPSP_->substitutionProcess().getNumberOfModels();
  }

  std::vector<size_t> getModelNumbers() const override
  {
    return pOPSP_->substitutionProcess().getModelNumbers();
  }

  LikelihoodCalculationSingleProcess& getLikelihoodCalculationSingleProcess()
  {
    return *pOPSP_->getLikelihoodCalculationSingleProcess();
  }

  const LikelihoodCalculationSingleProcess& getLikelihoodCalculationSingleProcess() const
  {
    return *pOPSP_->getLikelihoodCalculationSingleProcess();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_PHYLOMAPPINGS_ONEPROCESSSEQUENCESUBSTITUTIONMAPPING_H
