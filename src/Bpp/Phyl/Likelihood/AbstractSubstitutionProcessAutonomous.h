//
// File: AbstractSubstitutionProcessAutonomous.h
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

#ifndef BPP_PHYL_LIKELIHOOD_ABSTRACTSUBSTITUTIONPROCESSAUTONOMOUS_H
#define BPP_PHYL_LIKELIHOOD_ABSTRACTSUBSTITUTIONPROCESSAUTONOMOUS_H

#include "AbstractSubstitutionProcess.h"
#include "SubstitutionProcessAutonomous.h"

// From the STL:
#include <memory>

namespace bpp
{
/**
 * @brief A partial implementation of the SubstitutionProcess interface.
 *
 * This class handles a pointer toward a ParametrizableTree object, as well
 * as convenient arrays for storing previously computed probabilities.
 */
class AbstractSubstitutionProcessAutonomous :
  public virtual SubstitutionProcessAutonomous,
  public virtual AbstractSubstitutionProcess
{
protected:
  std::unique_ptr<ParametrizablePhyloTree> pTree_;

  std::shared_ptr<ModelScenario> modelScenario_;

protected:
  AbstractSubstitutionProcessAutonomous(ParametrizablePhyloTree* tree, const std::string& prefix = "");

  AbstractSubstitutionProcessAutonomous(const AbstractSubstitutionProcessAutonomous& asp);

  AbstractSubstitutionProcessAutonomous& operator=(const AbstractSubstitutionProcessAutonomous& asp);

public:
  const ParametrizablePhyloTree& getParametrizablePhyloTree() const { return *pTree_; }

  /**
   * @brief AbsractParametrizable interface
   *
   **/

  void fireParameterChanged(const ParameterList& pl);

  /**
   * @brief Return if process has ModelScenario.
   *
   **/
  bool hasModelScenario() const
  {
    return modelScenario_ != 0;
  }

  /**
   * @brief get the ModelScenario.
   *
   **/
  const ModelScenario& getModelScenario() const
  {
    return *modelScenario_;
  }

  /**
   * @brief set the ParametrizablePhyloTree.
   *
   * Will build a unique_ptr<ParametrizablePhyloTree> from the given PhyloTree
   *
   **/

  void setPhyloTree(const PhyloTree& phyloTree);

  /**
   * @brief set the ModelScenario.
   *
   **/

  virtual void setModelScenario(std::shared_ptr<ModelScenario> modelscenario) = 0;
};
} // end namespace bpp

#endif // BPP_PHYL_LIKELIHOOD_ABSTRACTSUBSTITUTIONPROCESSAUTONOMOUS_Hx
