//
// File: AbstractAutonomousSubstitutionProcess.h
// Authors: Laurent Guéguen
// Created: jeudi 16 décembre 2021, à 21h 48
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
 *
 */
  
class AbstractAutonomousSubstitutionProcess :
  public virtual AutonomousSubstitutionProcess,
  public virtual AbstractSubstitutionProcess
{
protected:
  std::shared_ptr<ParametrizablePhyloTree> pTree_;

  /**
   * @brief Root frequencies.
   */
  std::shared_ptr<FrequencySet> rootFrequencies_;
  
  std::shared_ptr<ModelScenario> modelScenario_;

protected:
  /*
   * @brief Builds using an optional pointer towards a PhyloTree.
   *
   * If the pointer is non-null, a ParametrizablePhyloTree will be
   * built from this PhyloTree, and owned by the
   * AbstractAutonomousSubstitutionProcess.
   *
   */
  AbstractAutonomousSubstitutionProcess(std::shared_ptr<const PhyloTree> tree = 0, std::shared_ptr<FrequencySet> rootFrequencies = 0, const std::string& prefix = "");

  /*
   * @brief Builds using a pointer towards a ParametrizablePhyloTree.
   * This pointer will be owned by the
   * AbstractAutonomousSubstitutionProcess.
   *
   */
  
  AbstractAutonomousSubstitutionProcess(std::shared_ptr<ParametrizablePhyloTree> tree, std::shared_ptr<FrequencySet> rootFrequencies = 0, const std::string& prefix = "");

  AbstractAutonomousSubstitutionProcess(const AbstractAutonomousSubstitutionProcess& asp);

  AbstractAutonomousSubstitutionProcess& operator=(const AbstractAutonomousSubstitutionProcess& asp);

public:

  /**
   * @brief sets the ParametrizablePhyloTree.
   *
   * Will build a unique_ptr<ParametrizablePhyloTree> from the given PhyloTree
   *
   **/

  void setPhyloTree(const PhyloTree& phyloTree);

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

  std::shared_ptr<const FrequencySet> getRootFrequencySet() const
  {
    return rootFrequencies_;
  }

  /**
   * @brief set the RootFrequency.
   *
   * @param rootfrequency The root frequencies to be associated with this instance.
   */

  void setRootFrequencySet(std::shared_ptr<FrequencySet> rootfrequency)
  {
    if (rootFrequencies_)
      getParameters_().deleteParameters(rootFrequencies_->getParameters().getParameterNames(),false);
    rootFrequencies_=rootfrequency;
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
