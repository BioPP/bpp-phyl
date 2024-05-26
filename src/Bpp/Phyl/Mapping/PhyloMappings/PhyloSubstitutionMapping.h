// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_PHYLOMAPPINGS_PHYLOSUBSTITUTIONMAPPING_H
#define BPP_PHYL_MAPPING_PHYLOMAPPINGS_PHYLOSUBSTITUTIONMAPPING_H


// From bpp-seq:
#include <Bpp/Numeric/ParameterAliasable.h>
#include <Bpp/Numeric/ParametrizableCollection.h>
#include "../SubstitutionRegister.h"
#include "../../Model/SubstitutionModel.h"
#include "../SubstitutionMappingTools.h"

#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
class PhyloSubstitutionMapping :
  public Clonable
{
public:
  PhyloSubstitutionMapping() {}
  virtual ~PhyloSubstitutionMapping() {}

  PhyloSubstitutionMapping* clone() const = 0;

public:
  /**
   * @brief Get the substitution model corresponding to a certain
   * branch and model class.
   *
   * @param branchId The id of the branch.
   * @param classIndex The model class index.
   */
  virtual std::shared_ptr<const TransitionModelInterface> getModel(unsigned int branchId, size_t classIndex) const = 0;

  /**
   * @brief Checks and sets the models with given parameters.
   */
  virtual bool matchParametersValues(const ParameterList& nullParams) = 0;

  /**
   * @brief Gets the parameters.
   */
  virtual const ParameterList& getParameters() const = 0;

  /**
   * @brief compute Normalizations
   * @param nullParams a list of null parameters
   * @param unresolvedOption  mgmt of gaps in the counts (default: counted as zeros)
   * @param verbose
   */

  virtual void computeNormalizations(const ParameterList& nullParams,
      short unresolvedOption = SubstitutionMappingTools::UNRESOLVED_ZERO,
      bool verbose = true) = 0;

  /**
   * @brief return if normalizations have been performed.
   */
  virtual bool normalizationsPerformed() const = 0;

  /**
   * @brief return if counts have been performed.
   */
  virtual bool countsPerformed() const = 0;

  /**
   * @brief ComputeCounts
   * @param unresolvedOption  mgmt of gaps in the counts (default: counted as zeros)
   * @param threshold 
   * @param verbose 
   */

  virtual void computeCounts(short unresolvedOption = SubstitutionMappingTools::UNRESOLVED_ZERO,
      double threshold = -1, bool verbose = true) = 0;

  /**
   * @brief Return the tree of counts
   */
  virtual ProbabilisticSubstitutionMapping& counts() = 0;

  virtual const ProbabilisticSubstitutionMapping& counts() const = 0;

  /**
   * @brief Return the tree of factors
   */
  virtual ProbabilisticSubstitutionMapping& normalizations() = 0;

  virtual const ProbabilisticSubstitutionMapping& normalizations() const = 0;

  /**
   * @brief change Distances
   */
  virtual void setDistances(const AlphabetIndex2& ndist) = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_PHYLOMAPPINGS_PHYLOSUBSTITUTIONMAPPING_H
