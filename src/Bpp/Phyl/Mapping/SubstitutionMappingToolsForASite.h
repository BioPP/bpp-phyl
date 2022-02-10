//
// File: SubstitutionMappingToolsForASite.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: mercredi 26 septembre 2018, ÃÂ  21h 11
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)
  
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

#ifndef BPP_PHYL_MAPPING_SUBSTITUTIONMAPPINGTOOLSFORASITE_H
#define BPP_PHYL_MAPPING_SUBSTITUTIONMAPPINGTOOLSFORASITE_H

#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

#include "../Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h"
#include "BranchedModelSet.h"
#include "OneJumpSubstitutionCount.h"
#include "ProbabilisticSubstitutionMapping.h"
#include "SubstitutionCount.h"

namespace bpp
{
/**
 * @brief Provide methods to compute substitution mappings.
 *
 * For now, only 4 methods are implemented, and provide probabilistic
 * substitution mappings.
 *
 * See:
 * Dutheil J, Pupko T, Jean-Marie A, Galtier N.
 * A model-based approach for detecting coevolving positions in a molecule.
 * Mol Biol Evol. 2005 Sep;22(9):1919-28. Epub 2005 Jun 8.
 *
 * @author Julien Dutheil
 */
class SubstitutionMappingToolsForASite
{
public:
  typedef std::map<const SubstitutionModel*, std::shared_ptr<SubstitutionCount > > t_Sm_Sc;
  typedef std::map<const SubstitutionRegister*, t_Sm_Sc > t_Sr_Sm_Sc;

private:
  static t_Sr_Sm_Sc m_Sr_Sm_Sc;

public:
  SubstitutionMappingToolsForASite() {}
  virtual ~SubstitutionMappingToolsForASite() {}

public:
  /**
   * @brief Methods to compute mapping Trees
   *
   * @{
   */

  /**
   * @brief Compute the substitutions tree for a particular dataset
   *
   * @param site              The site
   * @param rltc              A LikelihoodCalculationSingleProcess object.
   * @param reg               The SubstitutionRegister to use.
   * @param weights           Pointer to AlphabetIndex2 for weights
   *                          for all substitutions (default: null
   *                          means no weight),
   * @param distances         Pointer to AlphabetIndex2 for distances
   *                          for all substitutions (default: null
   *                          means each distance = 1),
   * @param threshold         value above which counts are considered
   *                          saturated (default: -1 means no threshold).
   * @param verbose           Print info to screen.
   * @return A tree <PhyloNode, PhyloBranchMapping>
   */
  static ProbabilisticSubstitutionMapping* computeCounts(
    size_t site,
    LikelihoodCalculationSingleProcess& rltc,
    const SubstitutionRegister& reg,
    std::shared_ptr<const AlphabetIndex2> weights = 0,
    std::shared_ptr<const AlphabetIndex2> distances = 0,
    double threshold = -1,
    bool verbose = true)
  {
    std::vector<uint> edgeIds = rltc.getSubstitutionProcess().getParametrizablePhyloTree()->getAllEdgesIndexes();
    return computeCounts(site, rltc, edgeIds, reg, weights, distances, threshold, verbose);
  }

  /**
   * @brief Compute the substitutions tree for a particular dataset
   *
   * @param site              The site
   * @param rltc              A LikelihoodCalculationSingleProcess object.
   * @param m_Sm_Sc           a map from a pointer to SubstitutionModel to a
   *                          shared_ptr to a SubstitutionCount
   * @param edgeIds           The Ids of the nodes the substitutions
   *                          are counted on. If empty, count substitutions
   *                          on all nodes.
   * @param threshold         value above which counts are considered
   *                          saturated (default: -1 means no threshold).
   * @param verbose           Print info to screen.
   * @return A tree <PhyloNode, PhyloBranchMapping>
   */

  static ProbabilisticSubstitutionMapping* computeCounts(
    size_t site,
    LikelihoodCalculationSingleProcess& rltc,
    t_Sm_Sc& m_Sm_Sc,
    const std::vector<uint>& edgeIds,
    double threshold = -1,
    bool verbose = true);

  /**
   * @brief Compute the substitutions tree for a particular dataset
   *
   * @param site              The site
   * @param rltc              A LikelihoodCalculationSingleProcess object.
   * @param edgeIds           The Ids of the nodes the substitutions
   *                          are counted on. If empty, count substitutions
   *                          on all nodes.
   * @param reg               The SubstitutionRegister to use.
   * @param weights           Pointer to AlphabetIndex2 for weights
   *                          for all substitutions (default: null
   *                          means no weight),
   * @param distances         Pointer to AlphabetIndex2 for distances
   *                          for all substitutions (default: null
   *                          means each distance = 1),
   * @param threshold         value above which counts are considered
   *                          saturated (default: -1 means no threshold).
   * @param verbose           Print info to screen.
   * @return A tree <PhyloNode, PhyloBranchMapping>
   */

  static ProbabilisticSubstitutionMapping* computeCounts(
    size_t site,
    LikelihoodCalculationSingleProcess& rltc,
    const std::vector<uint>& edgeIds,
    const SubstitutionRegister& reg,
    std::shared_ptr<const AlphabetIndex2> weights = 0,
    std::shared_ptr<const AlphabetIndex2> distances = 0,
    double threshold = -1,
    bool verbose = true);

  /**
   * @brief Compute the normalizations tree due to the models of "null"
   * process on each branch, for each register.
   *
   * @param site              The site
   * @param rltc              A LikelihoodCalculationSingleProcess object.
   * @param edgeIds           The Ids of the nodes the substitutions
   *                          are counted on. If empty, count substitutions
   *                          on all nodes.
   * @param nullModels        The "null" models used for normalization
   * @param reg               the Substitution Register
   * @param distances         Pointer to AlphabetIndex2 for distances
   *                          for all substitutions (default: null
   *                          means each distance = 1),
   * @param verbose           Display progress messages.
   * @return A tree <PhyloNode, PhyloBranchMapping> of normalization factors.
   */

  static ProbabilisticSubstitutionMapping* computeNormalizations(
    size_t site,
    LikelihoodCalculationSingleProcess& rltc,
    const std::vector<uint>& edgeIds,
    const BranchedModelSet* nullModels,
    const SubstitutionRegister& reg,
    std::shared_ptr<const AlphabetIndex2> distances = 0,
    bool verbose = true);

  /**
   * @brief Compute the normalizations tree due to the models of "null"
   * process on each branch, for each register.
   *
   * @param site              The site
   * @param rltc              A LikelihoodCalculationSingleProcess object.
   * @param nullModels        The "null" models used for normalization
   * @param reg               the Substitution Register
   * @param distances         Pointer to AlphabetIndex2 for distances
   *                          for all substitutions (default: null
   *                          means each distance = 1),
   * @param verbose           Display progress messages.
   * @return A tree <PhyloNode, PhyloBranchMapping> of normalization factors.
   */
  static ProbabilisticSubstitutionMapping* computeNormalizations(
    size_t site,
    LikelihoodCalculationSingleProcess& rltc,
    const BranchedModelSet* nullModels,
    const SubstitutionRegister& reg,
    std::shared_ptr<const AlphabetIndex2> distances = 0,
    bool verbose = true)
  {
    std::vector<uint> edgeIds = rltc.getSubstitutionProcess().getParametrizablePhyloTree()->getAllEdgesIndexes();
    return computeNormalizations(site, rltc, edgeIds, nullModels, reg, distances, verbose);
  }

  /**
   * @brief Compute the substitution counts tree, normalized by the
   * models of "null" process on each branch, for each register.
   *
   * @param site              The site
   * @param rltc              A LikelihoodCalculationSingleProcess object.
   * @param edgeIds           The Ids of the nodes the substitutions
   *                          are counted on. If empty, count substitutions
   *                          on all nodes.
   * @param nullModels        The "null" models used for normalization
   * @param reg               the Substitution Register
   * @param weights           Pointer to AlphabetIndex2 for weights
   *                          for all substitutions (default: null
   *                          means no weight),
   * @param distances         Pointer to AlphabetIndex2 for distances
   *                          for all substitutions (default: null
   *                          means each distance = 1),
   * @param perTimeUnit       If true, normalized counts are per unit of
   *                          time (otherwise they are multiplied by
   *                          the length of the branches) (default:
   *                          false).
   * @param siteSize          The length of a site, as considered as
   *                          a counting unit (default = 1)
   * @param threshold         value above which non-normalized counts are
   *                          considered saturated (default: -1
   *                          means no threshold).
   * @param verbose           Display progress messages.
   * @return A tree <PhyloNode, PhyloBranchMapping> of normalized counts..
   */

  static ProbabilisticSubstitutionMapping* computeNormalizedCounts(
    size_t site,
    LikelihoodCalculationSingleProcess& rltc,
    const std::vector<uint>& edgeIds,
    const BranchedModelSet* nullModels,
    const SubstitutionRegister& reg,
    std::shared_ptr<const AlphabetIndex2> weights = 0,
    std::shared_ptr<const AlphabetIndex2> distances = 0,
    bool perTimeUnit = false,
    uint siteSize = 1,
    double threshold = -1,
    bool verbose = true);

  static ProbabilisticSubstitutionMapping* computeNormalizedCounts(
    const ProbabilisticSubstitutionMapping* counts,
    const ProbabilisticSubstitutionMapping* factors,
    const std::vector<uint>& edgeIds,
    bool perTimeUnit = false,
    uint siteSize = 1);

  static ProbabilisticSubstitutionMapping* computeNormalizedCounts(
    size_t site,
    LikelihoodCalculationSingleProcess& rltc,
    const BranchedModelSet* nullModels,
    const SubstitutionRegister& reg,
    std::shared_ptr<const AlphabetIndex2> weights = 0,
    std::shared_ptr<const AlphabetIndex2> distances = 0,
    bool perTimeUnit = false,
    uint siteSize = 1,
    double threshold = -1,
    bool verbose = true)
  {
    std::vector<uint> edgeIds = rltc.getSubstitutionProcess().getParametrizablePhyloTree()->getAllEdgesIndexes();
    return computeNormalizedCounts(site, rltc, edgeIds, nullModels, reg, weights, distances, perTimeUnit, siteSize, threshold, verbose);
  }

  static ProbabilisticSubstitutionMapping* computeNormalizedCounts(
    const ProbabilisticSubstitutionMapping* counts,
    const ProbabilisticSubstitutionMapping* factors,
    bool perTimeUnit,
    uint siteSize = 1)
  {
    std::vector<uint> edgeIds = factors->getAllEdgesIndexes();
    return computeNormalizedCounts(counts, factors, edgeIds, perTimeUnit, siteSize);
  }

  /**
   *@}
   *
   */
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_SUBSTITUTIONMAPPINGTOOLSFORASITE_H
