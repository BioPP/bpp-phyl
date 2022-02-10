//
// File: BppPhylogeneticsApplication.h
// Authors:
//   Laurent GuÃÂ©guen, Julien Dutheil
// Created: 2021-06-15 15:14:00
//

/*
  Copyright or ÃÂ© or Copr. Development Team, (November 17, 2021)
  
  This software is a computer program whose purpose is to provide basal and
  utilitary classes. This file belongs to the Bio++ Project.
  
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

#ifndef BPP_PHYL_APP_BPPPHYLOGENETICSAPPLICATION_H
#define BPP_PHYL_APP_BPPPHYLOGENETICSAPPLICATION_H


#include "../Likelihood/PhyloLikelihoods/PhyloLikelihoodContainer.h"
#include "../Likelihood/SequenceEvolution.h"
#include "../Likelihood/SubstitutionProcessCollection.h"
#include "../Tree/PhyloTree.h"

// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>

using namespace std;

// From bpp-seq:
#include <Bpp/Seq/App/BppSequenceApplication.h>

namespace bpp
{
class BppPhylogeneticsApplication :
  public virtual BppSequenceApplication
{
public:
  BppPhylogeneticsApplication(int argc, char* argv[], const std::string& name) :
    BppApplication(argc, argv, name),
    BppSequenceApplication(argc, argv, name)
  {}

public:
  /**
   * @brief Methods to build objects
   *
   * @{
   */


  /**
   * @brief Get the std::map of initial phylo trees
   */
  virtual std::map<size_t, std::shared_ptr<PhyloTree> > getPhyloTreesMap(
    const std::map<size_t, const AlignedValuesContainer*>& mSites,
    std::map<std::string, std::string>& unparsedParams,
    const std::string& prefix = "input.",
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;

  /**
   * @brief get the collection of objects necessary to build
   * substitution processes.
   */
  virtual SubstitutionProcessCollection* getCollection(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    const std::map<size_t, const AlignedValuesContainer*>& mSites,
    const std::map<size_t, std::shared_ptr<PhyloTree> >& mpTree,
    std::map<std::string, std::string>& unparsedParams,
    const std::string& prefix = "input.",
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;

  virtual SubstitutionProcessCollection* getCollection(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    const std::map<size_t, const AlignedValuesContainer*>& mSites,
    std::map<std::string, std::string>& unparsedParams,
    const std::string& prefix = "input.",
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;

  /**
   * @brief get the substitution processes.
   */
  virtual std::map<size_t, SequenceEvolution*> getProcesses(
    SubstitutionProcessCollection& collection,
    std::map<std::string, std::string>& unparsedParams,
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;

  /**
   * @brief get the phylolikelihoods.
   */
  virtual std::shared_ptr<PhyloLikelihoodContainer> getPhyloLikelihoods(
    Context& context,
    std::map<size_t, SequenceEvolution*> mSeqEvol,
    SubstitutionProcessCollection& collection,
    const std::map<size_t,
                   const AlignedValuesContainer*>& mSites,
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;


  /**
   * @}
   */

  /**
   * @brief Method to have a clean likelihood (ie not saturated, nor infinite).
   */
  virtual void fixLikelihood(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    PhyloLikelihood* phylolik,
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;

  /**
   * @brief Display parameter values.
   *
   * @param tl the phylolikelihood
   * @param displaylL if log-likelihood is displayed (default: true)
   */
  virtual void displayParameters(
    const PhyloLikelihood& tl,
    bool displaylL = true) const;
};
} // end of namespace bpp;
#endif // BPP_PHYL_APP_BPPPHYLOGENETICSAPPLICATION_H
