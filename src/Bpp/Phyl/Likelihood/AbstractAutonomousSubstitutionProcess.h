// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_ABSTRACT_AUTONOMOUS_SUBSTITUTIONPROCESS_H
#define BPP_PHYL_LIKELIHOOD_ABSTRACT_AUTONOMOUS_SUBSTITUTIONPROCESS_H

#include "AbstractSubstitutionProcess.h"
#include "AutonomousSubstitutionProcess.h"

// From the STL:
#include <memory>

namespace bpp
{
/**
 * @brief A partial implementation of the SubstitutionProcess interface.
 *
 * This class OWNS a pointer toward a ParametrizableTree object, as
 * well as convenient arrays for storing previously computed
 * probabilities.
 */
class AbstractAutonomousSubstitutionProcess :
  public virtual AutonomousSubstitutionProcessInterface,
  public virtual AbstractSubstitutionProcess
{
protected:
  std::shared_ptr<ParametrizablePhyloTree> pTree_;

  /**
   * @brief Root frequencies.
   */
  std::shared_ptr<FrequencySetInterface> rootFrequencies_;
  
  std::shared_ptr<ModelScenario> modelScenario_;

protected:
  
  /**
   * @brief Builds using an optional pointer towards a PhyloTree.
   *
   * If the pointer is non-null, a ParametrizablePhyloTree will be
   * built from this PhyloTree, and owned by the
   * AbstractAutonomousSubstitutionProcess.
   *
   */
  AbstractAutonomousSubstitutionProcess(
      std::shared_ptr<const PhyloTree> tree = 0,
      std::shared_ptr<FrequencySetInterface> rootFrequencies = 0,
      const std::string& prefix = "");

  /**
   * @brief Builds using a pointer towards a ParametrizablePhyloTree.
   * This pointer will be owned by the
   * AbstractAutonomousSubstitutionProcess.
   */
  AbstractAutonomousSubstitutionProcess(
      std::shared_ptr<ParametrizablePhyloTree> tree,
      std::shared_ptr<FrequencySetInterface> rootFrequencies = 0,
      const std::string& prefix = "");

  AbstractAutonomousSubstitutionProcess(const AbstractAutonomousSubstitutionProcess& asp);

  AbstractAutonomousSubstitutionProcess& operator=(const AbstractAutonomousSubstitutionProcess& asp);

public:

  /**
   * @brief sets the ParametrizablePhyloTree.
   *
   * Will build a unique_ptr<ParametrizablePhyloTree> from the given PhyloTree
   */
  void setPhyloTree(const PhyloTree& phyloTree);

  const ParametrizablePhyloTree& parametrizablePhyloTree() const
  {
    return *pTree_;
  }

  std::shared_ptr<const ParametrizablePhyloTree> getParametrizablePhyloTree() const
  {
    return pTree_;
  }
  
  /**
   * @return true if the process has parametrized root frequencies (non-stationary model)
   */
  bool hasRootFrequencySet() const
  {
    return rootFrequencies_!=0;
  }

  const FrequencySetInterface& rootFrequencySet() const
  {
    return *rootFrequencies_;
  }

  std::shared_ptr<const FrequencySetInterface> getRootFrequencySet() const
  {
    return rootFrequencies_;
  }

  FrequencySetInterface& rootFrequencySet()
  {
    return *rootFrequencies_;
  }

  std::shared_ptr<FrequencySetInterface> getRootFrequencySet()
  {
    return rootFrequencies_;
  }
  
  /**
   * @brief set the RootFrequency.
   *
   * @param rootfrequency The root frequencies to be associated with this instance.
   */
  void setRootFrequencySet(std::shared_ptr<FrequencySetInterface> rootfrequency)
  {
    if (rootFrequencies_)
      getParameters_().deleteParameters(rootFrequencies_->getParameters().getParameterNames(),false);
    rootFrequencies_ = rootfrequency;
    if (rootFrequencies_)
      addParameters_(rootFrequencies_->getParameters());
  }

  /**
   * @brief Get the parameters corresponding to the root frequencies.
   *
   * @return The parameters corresponding to the root frequencies.
   */
  ParameterList getRootFrequenciesParameters(bool independent) const
  {
    if (!hasRootFrequencySet())
      return ParameterList();
    else
      return rootFrequencies_->getParameters();
  }


  /**
   * @brief AbsractParametrizable interface
   *
   **/

  void fireParameterChanged(const ParameterList& pl);

  /**
   * @brief get the ModelScenario.
   *
   **/
  std::shared_ptr<const ModelScenario> getModelScenario() const
  {
    return modelScenario_;
  }

  /**
   * @brief set the ModelScenario.
   *
   **/

  virtual void setModelScenario(std::shared_ptr<ModelScenario> modelScenario) = 0;
};
} // end namespace bpp

#endif // BPP_PHYL_LIKELIHOOD_ABSTRACT_AUTONOMOUS_SUBSTITUTIONPROCESS
