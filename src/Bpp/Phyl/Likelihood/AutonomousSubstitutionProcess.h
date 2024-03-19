// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_AUTONOMOUS_SUBSTITUTIONPROCESS_H
#define BPP_PHYL_LIKELIHOOD_AUTONOMOUS_SUBSTITUTIONPROCESS_H


#include "SubstitutionProcess.h"

// From the STL:
#include <memory>

namespace bpp
{
/**
 * @brief Interface for SubstitutionProcess objects that own their own
 * ParametrizablePhyloTree & Scenario.
 *
 */
class AutonomousSubstitutionProcessInterface :
  public virtual SubstitutionProcessInterface
{
public:
  /**
   * @brief set the ParametrizablePhyloTree.
   *
   * Will build a unique_ptr<ParametrizablePhyloTree> from the given PhyloTree
   *
   * @param phyloTree Tree, in the form of PhyloTree, to be associated with this instance.
   */
  virtual void setPhyloTree(const PhyloTree& phyloTree) = 0;

  /**
   * @brief set the RootFrequency.
   *
   * @param rootfrequency The root frequencies to be associated with this instance.
   */
  virtual void setRootFrequencySet(std::shared_ptr<FrequencySetInterface> rootfrequency) = 0;

  /**
   * @brief set the ModelScenario.
   *
   * @param modelScenario The scenario to be associated with this instance.
   */
  virtual void setModelScenario(std::shared_ptr<ModelScenario> modelScenario) = 0;
};
} // end namespace bpp

#endif // BPP_PHYL_LIKELIHOOD_AUTONOMOUS_SUBSTITUTIONPROCESS_H
