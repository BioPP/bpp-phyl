//
// File: AutonomousSubstitutionProcess.h
// Authors: Laurent Guéguen
// Created: jeudi 16 décembre 2021, à 21h 39
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
  class AutonomousSubstitutionProcess :
    public virtual SubstitutionProcess
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
    virtual void setRootFrequencySet(std::shared_ptr<FrequencySet> rootfrequency) = 0;

    /**
     * @brief set the ModelScenario.
     *
     * @param modelScenario The scenario to be associated with this instance.
     */
    virtual void setModelScenario(std::shared_ptr<ModelScenario> modelScenario) = 0;

  };
} // end namespace bpp

#endif // BPP_PHYL_LIKELIHOOD_AUTONOMOUS_SUBSTITUTIONPROCESS_H
