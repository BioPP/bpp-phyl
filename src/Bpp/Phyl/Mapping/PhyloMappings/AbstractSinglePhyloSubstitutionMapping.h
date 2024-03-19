// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_PHYLOMAPPINGS_ABSTRACTSINGLEPHYLOSUBSTITUTIONMAPPING_H
#define BPP_PHYL_MAPPING_PHYLOMAPPINGS_ABSTRACTSINGLEPHYLOSUBSTITUTIONMAPPING_H

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/ParametrizableCollection.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

#include "../../Tree/PhyloTree.h"
#include "../BranchedModelSet.h"
#include "../ProbabilisticSubstitutionMapping.h"
#include "PhyloSubstitutionMapping.h"

namespace bpp
{
struct ModelBranch
{
  std::shared_ptr<TransitionModelInterface> pMod_;
};


/**
 * @brief The AbstractSinglePhyloSubstitutionMapping class: substitution
 * mapping linked with a Single Process PhyloLikelihood
 */
class AbstractSinglePhyloSubstitutionMapping :
  virtual public BranchedModelSet,
  virtual public PhyloSubstitutionMapping,
  public AssociationTreeGlobalGraphObserver<uint, ModelBranch>
{
public:
  typedef AssociationTreeGlobalGraphObserver<uint, ModelBranch> modelTree;

private:
  std::shared_ptr<const SubstitutionRegisterInterface> pReg_;

  /**
   * @brief weights of the substitutions. If null, no weights are
   * used.
   */
  std::shared_ptr<const AlphabetIndex2> weights_;

  /**
   * @brief distances of the substitutions. If null, no distances are
   * used.
   */
  std::shared_ptr<const AlphabetIndex2> distances_;

protected:
  std::unique_ptr<ProbabilisticSubstitutionMapping> counts_;
  std::unique_ptr<ProbabilisticSubstitutionMapping> factors_;

private:
  /**
   * @brief A collection of Transition Models
   */
  ParametrizableCollection<TransitionModelInterface> modelColl_;

  /**
   *
   * @brief a map <model index, vector of branch ids>
   *
   */
  std::map<size_t, std::vector<uint>> mModBrid_;

public:
  AbstractSinglePhyloSubstitutionMapping(
      std::shared_ptr<TreeGlobalGraph> graph,
      std::shared_ptr<const SubstitutionRegisterInterface> reg,
      std::shared_ptr<const AlphabetIndex2> weights,
      std::shared_ptr<const AlphabetIndex2> distances) :
    modelTree(graph), pReg_(reg), weights_(weights), distances_(distances), counts_(), factors_(), modelColl_(), mModBrid_()
  {}

  AbstractSinglePhyloSubstitutionMapping(const AbstractSinglePhyloSubstitutionMapping& sppm);

  AbstractSinglePhyloSubstitutionMapping& operator=(const AbstractSinglePhyloSubstitutionMapping& sppm);

  virtual ~AbstractSinglePhyloSubstitutionMapping() {}

  /**
   * @brief From BranchedModelSet
   *
   * @{
   */
  std::shared_ptr<TransitionModelInterface> getModelForBranch(uint branchId) override
  {
    return (*getEdge(branchId)).pMod_;
  }

  std::shared_ptr<const TransitionModelInterface> getModelForBranch(uint branchId) const override
  {
    return (*getEdge(branchId)).pMod_;
  }

  std::shared_ptr<const TransitionModelInterface> getModel(unsigned int branchId, size_t classIndex) const override
  {
    return (*getEdge(branchId)).pMod_;
  }

  std::shared_ptr<TransitionModelInterface> getModel(size_t index)
  {
    return modelColl_[index];
  }

  std::shared_ptr<const TransitionModelInterface> getModel(size_t index) const override
  {
    return modelColl_[index];
  }

  std::vector<uint> getBranchesWithModel(size_t index) const override
  {
    return mModBrid_.at(index);
  }

  /**
   * @}
   */

  /**
   * @brief From PhyloSubstitutionMapping.
   *
   * @{
   */

  /**
   * @brief Return the tree of factors
   */
  bool normalizationsPerformed() const override
  {
    return factors_ != 0;
  }

  ProbabilisticSubstitutionMapping& normalizations() override
  {
    return *factors_;
  }

  const ProbabilisticSubstitutionMapping& normalizations() const override
  {
    return *factors_;
  }

  bool countsPerformed() const override
  {
    return counts_ != 0;
  }

  ProbabilisticSubstitutionMapping& counts() override
  {
    return *counts_;
  }

  const ProbabilisticSubstitutionMapping& counts() const override
  {
    return *counts_;
  }


  /**
   * @brief For registers
   */
  void setSubstitutionRegister(std::shared_ptr<const SubstitutionRegisterInterface> reg)
  {
    pReg_ = reg;
  }

  const SubstitutionRegisterInterface& substitutionRegister() const
  {
    return *pReg_;
  }

  std::shared_ptr<const SubstitutionRegisterInterface> getSubstitutionRegister() const
  {
    return pReg_;
  }

  bool matchParametersValues(const ParameterList& nullParams) override
  {
    return modelColl_.matchParametersValues(nullParams);
  }

  const ParameterList& getParameters() const override
  {
    return modelColl_.getParameters();
  }

  std::shared_ptr<const AlphabetIndex2> getWeights() const
  {
    return weights_;
  }

  std::shared_ptr<const AlphabetIndex2> getDistances() const
  {
    return distances_;
  }

  /**
   * @brief add a Substitition Model in the map, and on all branches
   * with given Ids.
   *
   * @param index the index of the model
   * @param model the model that will be COPIED.
   * @param brIds the Ids of the branches that will carry this model.
   */
  void addModel(size_t index, const TransitionModelInterface& model, Vuint brIds);

  /*
   * @brief change Distances
   *
   *  BEWARE: counts are not updated automatically
   */
  void setDistances(const AlphabetIndex2& ndist) override
  {
    distances_.reset(ndist.clone());
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_PHYLOMAPPINGS_ABSTRACTSINGLEPHYLOSUBSTITUTIONMAPPING_H
