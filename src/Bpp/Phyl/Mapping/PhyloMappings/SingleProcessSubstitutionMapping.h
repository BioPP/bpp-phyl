//
// File: SingleProcessSubstitutionMapping.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: dimanche 3 dÃÂ©cembre 2017, ÃÂ  13h 56
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

#ifndef BPP_PHYL_MAPPING_PHYLOMAPPINGS_SINGLEPROCESSSUBSTITUTIONMAPPING_H
#define BPP_PHYL_MAPPING_PHYLOMAPPINGS_SINGLEPROCESSSUBSTITUTIONMAPPING_H


#include "../../Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h"
#include "../ProbabilisticSubstitutionMapping.h"
#include "../SubstitutionMappingTools.h"

// #include "../SubstitutionMappingToolsForASite.h"

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
  SingleProcessPhyloLikelihood* pSPP_;

  /**
   * @brief Set the models of the BranchedModelSet to the adhoc
   * branches, for normalization.
   */
  void setBranchedModelSet_();

public:
  SingleProcessSubstitutionMapping(
      SingleProcessPhyloLikelihood& spp, 
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
  void computeCounts(double threshold = -1, bool verbose = true)
  {
    counts_ = SubstitutionMappingTools::computeCounts(
        getLikelihoodCalculationSingleProcess(),
	getSubstitutionRegister(),
	getWeights(),
	getDistances(),
	threshold,
	verbose);
  }

  // void computeCountsForASite(size_t site, double threshold = -1, bool verbose=true)
  // {
  //   counts_.reset(SubstitutionMappingToolsForASite::computeCounts(
  //                   site,
  //                   getLikelihoodCalculationSingleProcess(),
  //                   getRegister(),
  //                   getWeights(),
  //                   getDistances(),
  //                   threshold,
  //                   verbose));
  // }

  void computeNormalizations(const ParameterList& nullParams,
                             bool verbose = true);

  // void computeNormalizationsForASite(size_t site,
  //                                    const ParameterList& nullParams,
  //                                    bool verbose = true);

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
