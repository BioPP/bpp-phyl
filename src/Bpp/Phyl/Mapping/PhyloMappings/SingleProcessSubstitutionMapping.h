// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_PHYLOMAPPINGS_SINGLEPROCESSSUBSTITUTIONMAPPING_H
#define BPP_PHYL_MAPPING_PHYLOMAPPINGS_SINGLEPROCESSSUBSTITUTIONMAPPING_H


#include "../../Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h"
#include "../ProbabilisticSubstitutionMapping.h"
#include "../SubstitutionMappingTools.h"

#include "AbstractSinglePhyloSubstitutionMapping.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

namespace bpp
{
/**
 * @brief The SingleProcessSubstitutionMapping class: substitution
 * mapping linked with a SingleProcessPhyloLikelihood
 *
 */

class SingleProcessSubstitutionMapping :
  public AbstractSinglePhyloSubstitutionMapping,
  public std::enable_shared_from_this<SingleProcessSubstitutionMapping>
{
private:
  std::shared_ptr<SingleProcessPhyloLikelihood> pSPP_;

  /**
   * @brief Set the models of the BranchedModelSet to the adhoc
   * branches, for normalization.
   */
  void setBranchedModelSet_();

public:
  SingleProcessSubstitutionMapping(
      std::shared_ptr<SingleProcessPhyloLikelihood> spp,
      std::shared_ptr<SubstitutionRegisterInterface> reg,
      std::shared_ptr<const AlphabetIndex2> weights,
      std::shared_ptr<const AlphabetIndex2> distances,
      double threshold = -1,
      bool verbose = true);

  SingleProcessSubstitutionMapping(const SingleProcessSubstitutionMapping& sppm) :
    AbstractSinglePhyloSubstitutionMapping(sppm),
    pSPP_(sppm.pSPP_)
  {}

  SingleProcessSubstitutionMapping& operator=(const SingleProcessSubstitutionMapping& sppm)
  {
    AbstractSinglePhyloSubstitutionMapping::operator=(sppm);
    pSPP_ = sppm.pSPP_;
    return *this;
  }

  virtual ~SingleProcessSubstitutionMapping() {}

  SingleProcessSubstitutionMapping* clone() const { return new SingleProcessSubstitutionMapping(*this); }

  /*
   * @brief ComputeCounts
   *
   */
  void computeCounts(short unresolvedOption = SubstitutionMappingTools::UNRESOLVED_ZERO, double threshold = -1, bool verbose = true)
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

  void computeNormalizations(const ParameterList& nullParams,
      short unresolvedOption = SubstitutionMappingTools::UNRESOLVED_ZERO,
      bool verbose = true);

  /*
   * @brief Return the tree of counts
   *
   */
  size_t getNumberOfModels() const
  {
    return pSPP_->substitutionProcess().getNumberOfModels();
  }

  std::vector<size_t> getModelNumbers() const
  {
    return pSPP_->substitutionProcess().getModelNumbers();
  }

  LikelihoodCalculationSingleProcess& getLikelihoodCalculationSingleProcess()
  {
    return pSPP_->likelihoodCalculationSingleProcess();
  }

  const LikelihoodCalculationSingleProcess& getLikelihoodCalculationSingleProcess() const
  {
    return pSPP_->likelihoodCalculationSingleProcess();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_PHYLOMAPPINGS_SINGLEPROCESSSUBSTITUTIONMAPPING_H
