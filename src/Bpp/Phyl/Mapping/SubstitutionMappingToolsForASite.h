//
// File: SubstitutionMappingToolsForASite.h
// Created by: Laurent  Guéguen
// Created on: mercredi 26 septembre 2018, à 21h 11
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

#ifndef _SUBSTITUTIONMAPPINGTOOLS_FOR_A_SITE_H_
#define _SUBSTITUTIONMAPPINGTOOLS_FOR_A_SITE_H_

#include "ProbabilisticSubstitutionMapping.h"
#include "SubstitutionCount.h"
#include "OneJumpSubstitutionCount.h"
#include "../NewLikelihood/DataFlow/LikelihoodCalculationSingleProcess.h"
#include "../Model/BranchedModelSet.h"
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
  class SubstitutionMappingToolsForASite
  {
  public:
    
    typedef std::map<const SubstitutionModel*, std::unique_ptr<SubstitutionCount > > t_Sm_Sc;
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
     * and a particular site.
     *
     *
     * @param site              The site
     * @param rltc              A RecursiveLikelihoodTreeCalculation object.
     * @param substitutionCount The SubstitutionCount to use.
     * @param threshold         value above which counts are considered
     *                          saturated (default: -1 means no threshold).
     * @param verbose           Print info to screen.
     * @return A tree <PhyloNode, PhyloBranchMapping>
     */

    // static ProbabilisticSubstitutionMapping* computeCounts(
    //   size_t site,
    //   RecursiveLikelihoodTreeCalculation& rltc,
    //   SubstitutionCount& substitutionCount,
    //   double threshold = -1,
    //   bool verbose = true)
    // {
    //   std::vector<uint> nodeIds=rltc.getSubstitutionProcess().getParametrizablePhyloTree().getAllEdgesIndexes();
    //   return computeCounts(site, rltc, nodeIds, substitutionCount, threshold, verbose);
    // }

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
      std::vector<uint> nodeIds=rltc.getSubstitutionProcess().getParametrizablePhyloTree().getAllEdgesIndexes();
      return computeCounts(site, rltc, nodeIds, reg, weights, distances, threshold, verbose);
    }

    /**
     * @brief Compute the substitutions tree for a particular dataset
     *
     * @param site              The site
     * @param rltc              A LikelihoodCalculationSingleProcess object.
     * @param nodeIds           The Ids of the nodes the substitutions
     *                          are counted on. If empty, count substitutions
     *                          on all nodes.
     * @param substitutionCount The SubstitutionCount to use.
     * @param threshold         value above which counts are considered
     *                          saturated (default: -1 means no threshold).
     * @param verbose           Print info to screen.
     * @return A tree <PhyloNode, PhyloBranchMapping>
     */

    static ProbabilisticSubstitutionMapping* computeCounts(
      size_t site,
      LikelihoodCalculationSingleProcess& rltc,
      t_Sm_Sc& m_Sm_Sc,
      const std::vector<uint>& nodeIds,
      double threshold = -1,
      bool verbose = true);

    /**
     * @brief Compute the substitutions tree for a particular dataset
     *
     * @param site              The site
     * @param rltc              A LikelihoodCalculationSingleProcess object.
     * @param nodeIds           The Ids of the nodes the substitutions
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
      const std::vector<uint>& nodeIds,
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
     * @param nodeIds           The Ids of the nodes the substitutions
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
      const std::vector<uint>& nodeIds,
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
      std::vector<uint> nodeIds=rltc.getSubstitutionProcess().getParametrizablePhyloTree().getAllEdgesIndexes();
      return computeNormalizations(site, rltc, nodeIds, nullModels, reg, distances, verbose);
    }

    /**
     * @brief Compute the substitution counts tree, normalized by the
     * models of "null" process on each branch, for each register.
     *
     * @param site              The site
     * @param rltc              A LikelihoodCalculationSingleProcess object.
     * @param nodeIds           The Ids of the nodes the substitutions
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
      const std::vector<uint>& nodeIds,
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
      const std::vector<uint>& nodeIds,
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
      std::vector<uint> nodeIds=rltc.getSubstitutionProcess().getParametrizablePhyloTree().getAllEdgesIndexes();
      return computeNormalizedCounts(site, rltc, nodeIds, nullModels, reg, weights, distances, perTimeUnit, siteSize, threshold, verbose);
    }

    static ProbabilisticSubstitutionMapping* computeNormalizedCounts(
      const ProbabilisticSubstitutionMapping* counts,
      const ProbabilisticSubstitutionMapping* factors,
      bool perTimeUnit,
      uint siteSize = 1)
    {
      std::vector<uint> nodeIds=factors->getAllEdgesIndexes();
      return computeNormalizedCounts(counts, factors, nodeIds, perTimeUnit, siteSize);
    }

    /**
     * @brief This method computes for each branch the probability
     * that at least one jump occurred.
     *
     * Here 'jump' refer to a change in the model state. Depending on
     * the model, this might not be the same as a substitution (an
     * alphabet state change).
     *
     */

    // static ProbabilisticSubstitutionMapping* computeOneJumpCounts(
    //   size_t site,
    //   LikelihoodCalculationSingleProcess& rltc,
    //   bool verbose = true)
    // {
    //   OneJumpSubstitutionCount ojsm(0);
    //   return computeCounts(site, rltc, ojsm, 0);
    // }


    /**
     * @}
     *
     **/

    // /**
    //  *
    //  * @brief Methods to get trees of counts. These methods can use
    //  * raw or normalized counting trees (ie computed with
    //  * computeCounts and computeNormalizedCounts methods).
    //  *
    //  * @{
    //  *
    //  **/

    // /**
    //  * @brief Sum all sites substitutions a given type.
    //  *
    //  * @param counts The substitution map to use.
    //  *
    //  * @return A PhyloTree where branch lengths carry the branch counts.
    //  */

    // static PhyloTree* getTreeForType(const ProbabilisticSubstitutionMapping& counts,
    //                           size_t type);

    // static PhyloTree* getTreeForType(const ProbabilisticSubstitutionMapping& counts,
    //                                  const ProbabilisticSubstitutionMapping& factors,
    //                                  size_t type);

    // /**
    //  *
    //  * @brief Methods to get std::vectors of counts
    //  *
    //  * @{
    //  **/

    // /**
    //  * @brief PER BRANCH methods can use raw or normalized counting
    //  * trees (ie computed with computeCounts and
    //  * computeNormalizedCounts methods).
    //  *
    //  * @{
    //  *
    //  **/

    // /**
    //  * @brief Sum all type of substitutions for each branch of a given
    //  * position. 
    //  *
    //  * @param counts The substitution map to use.
    //  * @param site The site  for which the counts should be computed.
    //  *
    //  * @return A std::vector will all counts for all types of substitutions summed.
    //  */
    
    // static Vdouble getCountsPerBranch(
    //   const ProbabilisticSubstitutionMapping& counts,
    //   const std::vector<uint>& ids = Vuint(0));

    // static Vdouble getCountsPerBranch(
    //   const ProbabilisticSubstitutionMapping& counts,
    //   const ProbabilisticSubstitutionMapping& factors,
    //   const std::vector<uint>& ids = Vuint(0));

    // /**
    //  * @brief Return substitutions for each type of a given branch.
    //  *
    //  * @param counts The substitution map to use.
    //  * @param branchInd The id of the branch of the substitution tree for which the counts should be computed.
    //  * @return A std::vector will all counts summed for each type of substitutions.
    //  */
    
    // static Vdouble getCountsForBranchPerType(
    //   const ProbabilisticSubstitutionMapping& counts, uint branchId);

    // /**
    //  * @brief Return substitutions for each branch for each type. 
    //  *
    //  * @param counts The substitution map to use.
    //  * @param ids  The ids of the branches where the the substitutions
    //  * are counted (default : all ids)
    //  *
    //  * @return A std::vector will all counts for all branches for all types
    //  * of substitutions summed.
    //  */
    
    // static VVdouble getCountsPerBranchPerType(const ProbabilisticSubstitutionMappingForASite& counts, const std::vector<uint>& ids = Vuint(0));

    // static VVdouble getCountsPerTypePerBranch(const ProbabilisticSubstitutionMappingForASite& counts, const std::vector<uint>& ids = Vuint(0));

    // /**
    //  * @}
    //  *
    //  */

    // /**
    //  * @brief Methods that sum on branchs need raw mapping trees, or
    //  * raw mapping trees and normalization trees (ie computed with
    //  * computeCounts and computeNormalizations).
    //  *
    //  * @{
    //  *
    //  */
    
    // /**
    //  * @brief Sum all type of substitutions for each type of a given
    //  * position. 
    //  *
    //  * @param counts The substitution map to use.
    //  * @param site The site  for which the counts should be computed.
    //  * @param ids  The ids of the branches where the the substitutions
    //  *        are counted (default : all ids)
    //  *
    //  * @return A std::vector will all counts for all types of substitutions summed.
    //  */
    
    // static Vdouble getCountsPerType(
    //   const ProbabilisticSubstitutionMappingForASite& counts,
    //   const std::vector<uint>& ids = Vuint(0));

    // /**
    //  * @brief Sum and normalize all type of substitutions for each
    //  * type of a given position.
    //  *
    //  * @param site              The site
    //  * @param counts The substitution map to use.
    //  * @param factors The substitution normalization to use.
    //  * @param site The site  for which the counts should be computed.
    //  * @param ids  The ids of the branches where the the substitutions
    //  *        are counted (default : all ids)
    //  * @param perTimeUnit           If true, normalized counts are per unit of
    //  *                          time (otherwise they are multiplied by
    //  *                          the length of the branches).
    //  * @param siteSize          The length of a site, as considered as
    //  *                          a counting unit (default = 1)
    //  * @return A std::vector will all counts for all types of substitutions summed.
    //  */
    
    // static Vdouble getCountsPerType(
    //   const ProbabilisticSubstitutionMappingForASite& counts,
    //   const ProbabilisticSubstitutionMappingForASite& factors,
    //   bool perTimeUnit,
    //   uint siteSize = 1)
    // {
    //   return getCountsPerType(counts, factors, Vuint(0), perTimeUnit, siteSize);
    // }

    // static Vdouble getCountsPerType(
    //   const ProbabilisticSubstitutionMappingForASite& counts,
    //   const ProbabilisticSubstitutionMappingForASite& factors,
    //   const std::vector<uint>& ids,
    //   bool perTimeUnit,
    //   uint siteSize = 1);

    // /*
    //  *
    //  * @}
    //  */
    
    // /*
    //  *
    //  * @}
    //  *
    //  */
    
    // /*
    //  *
    //  * @brief Outputs of counts
    //  *
    //  * @{
    //  */

    // /**
    //  * @brief Write the substitutions std::vectors to a stream.
    //  *
    //  * @param substitutions The substitutions std::vectors to write.
    //  * (needed to know the position of each site in the dataset).
    //  * @param type          The type of substitutions to be output. See SubstitutionCount class.
    //  * Only one type of substitution can be output at a time.
    //  * @param out           The output stream where to write the std::vectors.
    //  * @throw IOException If an output error happens.
    //  */

    // static void writeToStream(
    //   const ProbabilisticSubstitutionMappingForASite& substitutions,
    //   size_t type,
    //   std::ostream& out);

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

#endif // _SUBSTITUTIONMAPPINGTOOLS_FOR_A_SITE_H_
