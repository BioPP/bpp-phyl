// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_PROCESSCOMPUTATIONTREE_H
#define BPP_PHYL_LIKELIHOOD_PROCESSCOMPUTATIONTREE_H


#include "../Tree/PhyloBranchParam.h"
#include "../Tree/PhyloTree.h"
#include "ParametrizablePhyloTree.h"
#include "SubstitutionProcess.h"

namespace bpp
{
/**
 * @brief Tree Organization of Computing Nodes
 *
 * Stores computation tools for all nodes, for all classes.
 *
 * This object has the parameters of the Tree and the rates
 * distribution, since it manages the SpeciationComputingNodes.
 *
 */

// Node specific DataFlow objects
class ProcessComputationNode :
  public PhyloNode
{
  /**
   * @brief the index of the species in the phyloTree matching this node.
   */
  const unsigned int speciesIndex_;

public:
  /**
   * @brief Build from a node in the phylo tree, with a specific
   * speciesIndex (because indexes are not the same as in the
   * ParametrizablePhyloTree.
   */
  ProcessComputationNode(const PhyloNode& node, unsigned int speciesIndex) :
    PhyloNode(node),
    speciesIndex_(speciesIndex) {}

  ProcessComputationNode(const ProcessComputationNode& node) :
    PhyloNode(node),
    speciesIndex_(node.speciesIndex_) {}

  unsigned int getSpeciesIndex() const
  {
    return speciesIndex_;
  }

  bool isSpeciation() const
  {
    auto prop = dynamic_cast<const NodeEvent*>(getProperty("event"));
    if (!prop)
      throw Exception("ProcessNode::isSpeciation : Node has no event associated: Node id " + TextTools::toString(getSpeciesIndex()));
    return prop->isSpeciation();
  }

  bool isMixture() const
  {
    auto prop = dynamic_cast<const NodeEvent*>(getProperty("event"));
    if (!prop)
      throw Exception("ProcessNode::isMixture : Node has no event associated: Node id " + TextTools::toString(getSpeciesIndex()));
    return prop->isMixture();
  }
};

using ProcessComputationNodeRef = std::shared_ptr<ProcessComputationNode>;

// Class for the edges: In case of mixture nodes, the probability is
// set through submodel probability
class ProcessComputationEdge
{
private:
  /**
   * @brief Model carried by the branch
   */
  std::shared_ptr<const BranchModelInterface> model_;

  /**
   * @brief Number of the model carried by the branch
   */
  unsigned int nmodel_;

  /**
   * @brief numbers of the submodels used, if any.
   *
   * In practice, only vectors of size <=1 are implemented (ie
   * individual submodels), to be fixed later.
   */
  std::vector<unsigned int> vSubNb_;

  /**
   * @brief use the probability associated to the edge in case of
   * mixture. If true, the probability is used and not the
   * transition probabilities.
   */
  bool useProb_;

  /**
   * @brief the index of the species in the phyloTree matching this node.
   */
  const unsigned int speciesIndex_;

public:
  ProcessComputationEdge(
      std::shared_ptr<const BranchModelInterface> model,
      unsigned int nmodel,
      unsigned int speciesIndex,
      bool useProb = false,
      const std::vector<unsigned int>& vNb = Vuint(0)) :
    model_(model),
    nmodel_(nmodel),
    vSubNb_(vNb),
    useProb_(useProb),
    speciesIndex_(speciesIndex)
  {}

  ProcessComputationEdge(const ProcessComputationEdge& edge) :
    model_(edge.model_),
    nmodel_(edge.nmodel_),
    vSubNb_(edge.vSubNb_),
    useProb_(edge.useProb_),
    speciesIndex_(edge.speciesIndex_)
  {}

  std::shared_ptr<const BranchModelInterface> getModel() const
  {
    return model_;
  }

  unsigned int getModelNumber() const
  {
    return nmodel_;
  }

  unsigned int getSpeciesIndex() const
  {
    return speciesIndex_;
  }

  const std::vector<unsigned int>& subModelNumbers() const
  {
    return vSubNb_;
  }

  bool useProb() const
  {
    return useProb_;
  }
};

using BaseTree = AssociationTreeGlobalGraphObserver<ProcessComputationNode, ProcessComputationEdge>;

class ProcessComputationTree :
  public BaseTree
{
private:
  std::shared_ptr<const SubstitutionProcessInterface> process_;

  /**
   * @brief Build rest of the tree under given father (event on father will
   * be set in this method)
   */
  void buildFollowingPath_(
      std::shared_ptr<ProcessComputationNode> father,
      const ModelPath& path);

  void buildFollowingScenario_(
      std::shared_ptr<ProcessComputationNode> father,
      const ModelScenario& scenario,
      std::map<std::shared_ptr<MixedTransitionModelInterface>,
      unsigned int>& mMrca);

public:
  /**
   * @brief construction of a ProcessComputationTree from a SubstitutionProcess
   *
   * @param process the SubstitutionProcess
   */
  ProcessComputationTree(std::shared_ptr<const SubstitutionProcessInterface> process);
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PROCESSCOMPUTATIONTREE_H
