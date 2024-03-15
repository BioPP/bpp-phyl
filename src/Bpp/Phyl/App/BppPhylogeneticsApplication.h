// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface> >& mSites,
    std::map<std::string, std::string>& unparsedParams,
    const std::string& prefix = "input.",
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;

  /**
   * @brief get the collection of objects necessary to build
   * substitution processes.
   */
  virtual std::unique_ptr<SubstitutionProcessCollection> getCollection(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface> >& mSites,
    const std::map<size_t, std::shared_ptr<PhyloTree> >& mpTree,
    std::map<std::string, std::string>& unparsedParams,
    const std::string& prefix = "input.",
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;

  virtual std::unique_ptr<SubstitutionProcessCollection> getCollection(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface> >& mSites,
    std::map<std::string, std::string>& unparsedParams,
    const std::string& prefix = "input.",
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;

  /**
   * @brief get the substitution processes.
   */
  virtual std::map<size_t, std::unique_ptr<SequenceEvolution> > getProcesses(
    std::shared_ptr<SubstitutionProcessCollection> collection,
    std::map<std::string, std::string>& unparsedParams,
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;

  /**
   * @brief get the phylolikelihoods.
   */
  virtual std::shared_ptr<PhyloLikelihoodContainer> getPhyloLikelihoods(
    Context& context,
    std::map<size_t, std::shared_ptr<SequenceEvolution> > mSeqEvol,
    std::shared_ptr<SubstitutionProcessCollection> collection,
    const std::map<size_t,
                   std::shared_ptr<const AlignmentDataInterface> >& mSites,
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;


  /**
   * @}
   */

  /**
   * @brief Method to have a clean likelihood (ie not saturated, nor infinite).
   */
  virtual void fixLikelihood(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    std::shared_ptr<PhyloLikelihoodInterface> phylolik,
    const std::string& suffix = "",
    bool suffixIsOptional = true) const;

  /**
   * @brief Display parameter values.
   *
   * @param tl the phylolikelihood
   * @param displaylL if log-likelihood is displayed (default: true)
   */
  virtual void displayParameters(
    const PhyloLikelihoodInterface& tl,
    bool displaylL = true) const;
};
} // end of namespace bpp;
#endif // BPP_PHYL_APP_BPPPHYLOGENETICSAPPLICATION_H
