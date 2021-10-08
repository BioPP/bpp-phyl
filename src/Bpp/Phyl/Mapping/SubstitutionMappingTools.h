//
// File: SubstitutionMappingTools.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 13:04 2006
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

   This software is a computer program whose purpose is to provide classes
   for phylogenetic data analysis.

   This software is governed by the CeCILL  license under French law and
   abiding by the rules of distribution of free software.  You can  use,
   modify and/ or redistribute the software under the terms of the CeCILL
   license as circulated by CEA, CNRS and INRIA at the following URL
   "http://www.cecill.info".

   As a counterpart to the access to the source code and  rights to copy,
   modify and redistribute granted by the license, users are provided only
   with a limited warranty  and the software's author,  the holder of the
   economic rights,  and the successive licensors  have only  limited
   liability.

   In this respect, the user's attention is drawn to the risks associated
   with loading,  using,  modifying and/or developing or reproducing the
   software by the user in light of its specific status of free software,
   that may mean  that it is complicated to manipulate,  and  that  also
   therefore means  that it is reserved for developers  and  experienced
   professionals having in-depth computer knowledge. Users are therefore
   encouraged to load and test the software's suitability as regards their
   requirements in conditions enabling the security of their systems and/or
   data to be ensured and,  more generally, to use and operate it in the
   same conditions as regards security.

   The fact that you are presently reading this means that you have had
   knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef _SUBSTITUTIONMAPPINGTOOLS_H_
#define _SUBSTITUTIONMAPPINGTOOLS_H_

#include "ProbabilisticSubstitutionMapping.h"
#include "SubstitutionCount.h"
#include "OneJumpSubstitutionCount.h"
#include "../NewLikelihood/DataFlow/LikelihoodCalculationSingleProcess.h"
#include "BranchedModelSet.h"
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

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
class SubstitutionMappingTools
{
public:
  SubstitutionMappingTools() {}
  virtual ~SubstitutionMappingTools() {}

public:
  /**
   * @brief Methods to compute mapping Trees
   *
   * @{
   */

  /**
   * @brief Compute the substitutions tree for a particular dataset
   *
   * @param rltc              A LikelihoodCalculationSingleProcess object.
   * @param substitutionCount The SubstitutionCount to use.
   * @param threshold         value above which counts are considered
   *                          saturated (default: -1 means no threshold).
   * @param verbose           Print info to screen.
   * @return A tree <PhyloNode, PhyloBranchMapping>
   */
  static ProbabilisticSubstitutionMapping* computeCounts(
    LikelihoodCalculationSingleProcess& rltc,
    SubstitutionCount& substitutionCount,
    double threshold = -1,
    bool verbose = true)
  {
    std::vector<uint> edgeIds = rltc.getSubstitutionProcess().getParametrizablePhyloTree().getAllEdgesIndexes();
    return computeCounts(rltc, edgeIds, substitutionCount, threshold, verbose);
  }

  /**
   * @brief Compute the substitutions tree for a particular dataset
   *
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
   *
   * @return A tree <PhyloNode, PhyloBranchMapping>
   */
  static ProbabilisticSubstitutionMapping* computeCounts(
    LikelihoodCalculationSingleProcess& rltc,
    const SubstitutionRegister& reg,
    std::shared_ptr<const AlphabetIndex2> weights = 0,
    std::shared_ptr<const AlphabetIndex2> distances = 0,
    double threshold = -1,
    bool verbose = true)
  {
    std::vector<uint> edgeIds = rltc.getSubstitutionProcess().getParametrizablePhyloTree().getAllEdgesIndexes();
    return computeCounts(rltc, edgeIds, reg, weights, distances, threshold, verbose);
  }

  /**
   * @brief Compute the substitutions tree for a particular dataset
   *
   * @param rltc              A LikelihoodCalculationSingleProcess object.
   * @param speciesIds        The Species Ids of the edges the substitutions
   *                          are counted on.
   * @param substitutionCount The SubstitutionCount to use.
   * @param threshold         value above which counts are considered
   *                          saturated (default: -1 means no threshold).
   * @param verbose           Print info to screen.
   *
   * @return A tree <PhyloNode, PhyloBranchMapping>
   */

  static ProbabilisticSubstitutionMapping* computeCounts(
    LikelihoodCalculationSingleProcess& rltc,
    const std::vector<uint>& speciesIds,
    SubstitutionCount& substitutionCount,
    double threshold = -1,
    bool verbose = true);

  /**
   * @brief Compute the substitutions tree for a particular dataset
   *
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
   *
   * @return A tree <PhyloNode, PhyloBranchMapping>
   */

  static ProbabilisticSubstitutionMapping* computeCounts(
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
   *
   * @return A tree <PhyloNode, PhyloBranchMapping> of normalization factors.
   */

  static ProbabilisticSubstitutionMapping* computeNormalizations(
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
   * @param rltc              A LikelihoodCalculationSingleProcess object.
   * @param nullModels        The "null" models used for normalization
   * @param reg               the Substitution Register
   * @param distances         Pointer to AlphabetIndex2 for distances
   *                          for all substitutions (default: null
   *                          means each distance = 1),
   * @param verbose           Display progress messages.
   *
   * @return A tree <PhyloNode, PhyloBranchMapping> of normalization factors.
   */
  static ProbabilisticSubstitutionMapping* computeNormalizations(
    LikelihoodCalculationSingleProcess& rltc,
    const BranchedModelSet* nullModels,
    const SubstitutionRegister& reg,
    std::shared_ptr<const AlphabetIndex2> distances = 0,
    bool verbose = true)
  {
    std::vector<uint> edgeIds = rltc.getSubstitutionProcess().getParametrizablePhyloTree().getAllEdgesIndexes();
    return computeNormalizations(rltc, edgeIds, nullModels, reg, distances, verbose);
  }

  /**
   * @brief Compute the substitution counts tree, normalized by the
   * models of "null" process on each branch, for each register.
   *
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
   *
   * @return A tree <PhyloNode, PhyloBranchMapping> of normalized counts..
   */

  static ProbabilisticSubstitutionMapping* computeNormalizedCounts(
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
    std::vector<uint> edgeIds = rltc.getSubstitutionProcess().getParametrizablePhyloTree().getAllEdgesIndexes();
    return computeNormalizedCounts(rltc, edgeIds, nullModels, reg, weights, distances, perTimeUnit, siteSize, threshold, verbose);
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
   * @brief This method computes for each site and for each branch
   * the probability that at least one jump occurred.
   *
   * Here 'jump' refer to a change in the model state. Depending on
   * the model, this might not be the same as a substitution (an
   * alphabet state change).
   *
   */
  static ProbabilisticSubstitutionMapping* computeOneJumpCounts(
    LikelihoodCalculationSingleProcess& rltc,
    bool verbose = true)
  {
    OneJumpSubstitutionCount ojsm(0);
    return computeCounts(rltc, ojsm, 0);
  }


  /**
   * @}
   *
   **/

  /**
   *
   * @brief Methods to get trees of counts. These methods can use
   * raw or normalized counting trees (ie computed with
   * computeCounts and computeNormalizedCounts methods).
   *
   * @{
   *
   **/

  /**
   * @brief Sum all sites substitutions a given type.
   *
   * @param counts The substitution map to use.
   * @param type  The number of the type
   *
   * @return A PhyloTree where branch lengths carry the branch counts.
   */

  static PhyloTree* getTreeForType(const ProbabilisticSubstitutionMapping& counts,
                                   size_t type);

  /**
   * @brief Sum all sites substitutions a given type.
   *
   * @param counts The substitution map to use.
   * @param factors The substitution normalization to use.
   * @param type  The number of the type
   *
   * @return A PhyloTree where branch lengths carry the branch counts.
   */
  static PhyloTree* getTreeForType(const ProbabilisticSubstitutionMapping& counts,
                                   const ProbabilisticSubstitutionMapping& factors,
                                   size_t type);

  /**
   *
   * @brief Methods to get std::vectors of counts
   *
   * @{
   **/

  /**
   * @brief PER BRANCH methods can use raw or normalized counting
   * trees (ie computed with computeCounts and
   * computeNormalizedCounts methods).
   *
   * @{
   *
   **/

  /**
   * @brief Sum all type of substitutions for each site for each
   * branch for each type.
   *
   * @param counts The substitution map to use.
   * @param ids  The ids of the branches where the the substitutions
   * are counted (default : all ids)
   *
   * @return A std::vector will all counts for all types of substitutions.
   */

  static VVVdouble getCountsPerSitePerBranchPerType(
    const ProbabilisticSubstitutionMapping& counts, const std::vector<uint>& ids = Vuint(0));


  /**
   * @brief Sum all type of substitutions for each branch of a given
   * position.
   *
   * @param counts The substitution map to use.
   * @param site The site  for which the counts should be computed.
   *
   * @return A std::vector will all counts for all types of substitutions summed.
   */

  static Vdouble getCountsForSitePerBranch(
    const ProbabilisticSubstitutionMapping& counts, size_t site);

  static Vdouble getCountsForSitePerBranch(
    const ProbabilisticSubstitutionMapping& counts,
    const ProbabilisticSubstitutionMapping& factors,
    size_t site);

  /**
   * @brief Sum all type of substitutions for each site for each branch.
   *
   * @param counts The substitution map to use.
   * @param ids  The ids of the branches where the the substitutions
   * are counted (default : all ids)
   *
   * @return A std::vector will all counts for all types of substitutions summed.
   */

  static VVdouble getCountsPerSitePerBranch(
    const ProbabilisticSubstitutionMapping& counts, const std::vector<uint>& ids = Vuint(0));

  static VVdouble getCountsPerSitePerBranch(
    const ProbabilisticSubstitutionMapping& counts,
    const ProbabilisticSubstitutionMapping& factors,
    const std::vector<uint>& ids = Vuint(0));

  /**
   * @brief Sum all sites substitutions for each type of a given branch.
   *
   * @param counts The substitution map to use.
   * @param branchId The id of the branch of the substitution tree for which the counts should be computed.
   * @return A std::vector will all counts summed for each type of substitutions.
   */

  static Vdouble getCountsForBranchPerType(
    const ProbabilisticSubstitutionMapping& counts, uint branchId);

  /**
   * @brief Sum all sites substitutions for each branch for each type.
   *
   * @param counts The substitution map to use.
   * @param ids  The ids of the branches where the the substitutions
   * are counted (default : all ids)
   *
   * @return A std::vector will all counts for all branches for all types
   * of substitutions summed.
   */

  static VVdouble getCountsPerBranchPerType(const ProbabilisticSubstitutionMapping& counts, const std::vector<uint>& ids = Vuint(0));

  static VVdouble getCountsPerTypePerBranch(const ProbabilisticSubstitutionMapping& counts, const std::vector<uint>& ids = Vuint(0));

  /**
   * @brief Compute the sum over all branches of the counts per type
   * per branch.
   *
   * @param rltc              A LikelihoodCalculationSingleProcess object.
   * @param ids               The vector of numbers of the nodes of the tree
   * @param reg               the Substitution Register
   *
   *           If the SubstitutionRegister is a non-stationary
   *           CategorySubstitutionRegister, a correction is made.
   *
   * @param weights           Pointer to AlphabetIndex2 for weights
   *                          for all substitutions (default: null
   *                          means no weight),
   * @param distances         Pointer to AlphabetIndex2 for distances
   *                          for all substitutions (default: null
   *                          means each distance = 1),
   * @param threshold         value above which counts are considered
   *                          saturated (default: -1 means no threshold).
   * @param verbose Display progress messages.
   */

  static VVdouble computeCountsPerTypePerBranch(
    LikelihoodCalculationSingleProcess& rltc,
    const std::vector<uint>& ids,
    const SubstitutionRegister& reg,
    std::shared_ptr<const AlphabetIndex2> weights = 0,
    std::shared_ptr<const AlphabetIndex2> distances = 0,
    double threshold = -1,
    bool verbose = true);

  /**
   * @}
   *
   */

  /**
   * @brief Methods that sum on branchs need raw mapping trees, or
   * raw mapping trees and normalization trees (ie computed with
   * computeCounts and computeNormalizations).
   *
   * @{
   *
   */

  /**
   * @brief Sum all type of substitutions for each type of a given
   * position.
   *
   * @param counts The substitution map to use.
   * @param site The site  for which the counts should be computed.
   * @param ids  The ids of the branches where the the substitutions
   *        are counted (default : all ids)
   *
   * @return A std::vector will all counts for all types of substitutions summed.
   */

  static Vdouble getCountsForSitePerType(
    const ProbabilisticSubstitutionMapping& counts,
    size_t site,
    const std::vector<uint>& ids = Vuint(0));

  /**
   * @brief Sum and normalize all type of substitutions for each
   * type of a given position on all nodes.
   *
   * @param counts The substitution map to use.
   * @param factors The substitution normalization to use.
   * @param site The site  for which the counts should be computed.
   * @param perTimeUnit           If true, normalized counts are per unit of
   *                          time (otherwise they are multiplied by
   *                          the length of the branches).
   * @param siteSize          The length of a site, as considered as
   *                          a counting unit (default = 1)
   * @return A std::vector will all counts for all types of substitutions summed.
   */

  static Vdouble getCountsForSitePerType(
    const ProbabilisticSubstitutionMapping& counts,
    const ProbabilisticSubstitutionMapping& factors,
    size_t site,
    bool perTimeUnit,
    uint siteSize = 1);

  /**
   * @brief Sum and normalize all type of substitutions for each
   * type of a given position on a set of nodes.
   *
   * @param counts The substitution map to use.
   * @param factors The substitution normalization to use.
   * @param site The site  for which the counts should be computed.
   * @param ids  The ids of the branches where the the substitutions
   *        are counted (default : all ids)
   * @param perTimeUnit           If true, normalized counts are per unit of
   *                          time (otherwise they are multiplied by
   *                          the length of the branches).
   * @param siteSize          The length of a site, as considered as
   *                          a counting unit (default = 1)
   * @return A std::vector will all counts for all types of substitutions summed.
   */

  static Vdouble getCountsForSitePerType(
    const ProbabilisticSubstitutionMapping& counts,
    const ProbabilisticSubstitutionMapping& factors,
    size_t site,
    const std::vector<uint>& ids,
    bool perTimeUnit,
    uint siteSize = 1);

  /**
   * @brief Sum all type of substitutions for each site for each type.
   *
   * @param counts The substitution map to use.
   * @param ids  The ids of the branches where the the substitutions
   *        are counted (default: all ids)
   * @return A std::vector will all counts for all types of substitutions summed.
   */

  static VVdouble getCountsPerSitePerType(
    const ProbabilisticSubstitutionMapping& counts,
    const std::vector<uint>& ids = Vuint(0));

  /**
   * @brief Sum and normalize all type of substitutions for each
   * site for each type on all nodes.
   *
   * @param counts The substitution map to use.
   * @param factors The substitution normalization to use.
   * @param perTimeUnit       If true, normalized counts are per unit of
   *                          time (otherwise they are multiplied by
   *                          the length of the branches).
   * @param siteSize          The length of a site, as considered as
   *                          a counting unit (default = 1)
   * @return A std::vector will all counts for all types of substitutions summed.
   */

  static VVdouble getCountsPerSitePerType(
    const ProbabilisticSubstitutionMapping& counts,
    const ProbabilisticSubstitutionMapping& factors,
    bool perTimeUnit,
    uint siteSize = 1);

  /**
   * @brief Sum and normalize all type of substitutions for each
   * site for each type on a set of nodes.
   *
   * @param counts The substitution map to use.
   * @param factors The substitution normalization to use.
   * @param ids  The ids of the branches where the the substitutions
   * are counted (default: all ids)
   * @param perTimeUnit       If true, normalized counts are per unit of
   *                          time (otherwise they are multiplied by
   *                          the length of the branches).
   * @param siteSize          The length of a site, as considered as
   *                          a counting unit (default = 1)
   * @return A std::vector will all counts for all types of substitutions summed.
   */

  static VVdouble getCountsPerSitePerType(
    const ProbabilisticSubstitutionMapping& counts,
    const ProbabilisticSubstitutionMapping& factors,
    const std::vector<uint>& ids,
    bool perTimeUnit,
    uint siteSize = 1);

  /**
   * @brief Compute the norm of a substitution std::vector for a given position.
   *
   * The norm is computed as:
   * @f$ N_i = \sqrt{\left(\sum_l {\left(\sum_t n_{l, i, t}\right)}^2\right)}@f$,
   * where @f$n_{l, i, t}@f$ is the number of substitutions of type t on site i on branch l, obtained using the () operator for the SubstitutionMapping object.
   *
   * @param counts The substitution map to use.
   * @param site The position for which the norm should be computed.
   * @return The norm of the substitution std::vector.
   */

  static double getNormForSite(
    const ProbabilisticSubstitutionMapping& counts,
    size_t site);

  /*
   *
   * @}
   */

  /*
   *
   * @}
   *
   */

  /*
   *
   * @brief Outputs of counts
   *
   * @{
   */

  /**
   * @brief Output Per Site Per Branch
   */
  static void outputPerSitePerBranch(const std::string& filename,
                                     const std::vector<uint>& ids,
                                     const VVdouble& counts);

  /**
   * @brief Output Per Site Per Type
   */
  static void outputPerSitePerType(const std::string& filename,
                                   const SubstitutionRegister& reg,
                                   const VVdouble& counts);

  /**
   * @brief Output Per Site Per Branch Per Type
   */
  static void outputPerSitePerBranchPerType(const std::string& filenamePrefix,
                                            const std::vector<uint>& ids,
                                            const SubstitutionRegister& reg,
                                            const VVVdouble& counts);


  /**
   * @brief Write the substitutions std::vectors to a stream.
   *
   * @param substitutions The substitutions std::vectors to write.
   * @param sites         The dataset associated to the std::vectors
   * (needed to know the position of each site in the dataset).
   * @param type          The type of substitutions to be output. See SubstitutionCount class.
   * Only one type of substitution can be output at a time.
   * @param out           The output stream where to write the std::vectors.
   * @throw IOException If an output error happens.
   */

  static void writeToStream(
    const ProbabilisticSubstitutionMapping& substitutions,
    const AlignedValuesContainer& sites,
    size_t type,
    std::ostream& out);

  /**
   * @brief Read the substitutions std::vectors from a stream.
   *
   * @param in            The input stream where to read the std::vectors.
   * @param substitutions The mapping object to fill.
   * @param type          The type of substitutions that are read. Should be in supported by the substittuion count obect assiciated to the mapping, if any.
   * @throw IOException If an input error happens.
   */

  static void readFromStream(std::istream& in, ProbabilisticSubstitutionMapping & substitutions, size_t type);

  /*
   *
   *@}
   */


  /**
   *@}
   *
   */
};
} // end of namespace bpp.

#endif// _SUBSTITUTIONMAPPINGTOOLS_H_
