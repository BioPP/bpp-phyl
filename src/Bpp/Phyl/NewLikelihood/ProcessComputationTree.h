//
// File: ProcessComputationTree.h
// Created by: Laurent Guéguen
// Created on: Sat Dec 30 12:48 2006
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _PROCESS_COMPUTATION_TREE_H_
#define _PROCESS_COMPUTATION_TREE_H_

#include "../Tree/PhyloTree.h"
#include "../Tree/PhyloBranchParam.h"

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
  /*
   * @brief the index of the species in the phyloTree matching this node.
   *
   */

  const uint speciesIndex_;

public:
  /*
   * @brief Build from a node in the phylo tree, with a specific
   * speciesIndex (because indexes are not the same as in the
   * ParametrizablePhyloTree.
   *
   */

  ProcessComputationNode(const PhyloNode& node, uint speciesIndex) :
    PhyloNode(node),
    speciesIndex_(speciesIndex) {}

  ProcessComputationNode(const ProcessComputationNode& node) :
    PhyloNode(node),
    speciesIndex_(node.speciesIndex_) {}

  uint getSpeciesIndex() const
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
  /*
   * Model carried by the branch
   *
   */

  const BranchModel* model_;

  /*
   * Number of the model carried by the branch
   *
   */

  uint nmodel_;

  /*
   * @brief numbers of the submodels used, if any.
   *
   * In practice, only vectors of size <=1 are implemented (ie
   * individual submodels), to be fixed later.
   */

  std::vector<uint> vSubNb_;

  /*
   * @brief use the probability associated to the edge in case of
   * mixture. If true, the probability is used and not the
   * transition probabilities.
   *
   */

  bool useProb_;

  /*
   * @brief the index of the species in the phyloTree matching this node.
   *
   */

  const uint speciesIndex_;

public:
  ProcessComputationEdge(const BranchModel* model, uint nmodel, uint speciesIndex, bool useProb = false, const std::vector<uint>& vNb = Vuint(0)) :
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

  const BranchModel* getModel() const
  {
    return model_;
  }

  uint getModelNumber() const
  {
    return nmodel_;
  }

  uint getSpeciesIndex() const
  {
    return speciesIndex_;
  }

  const std::vector<uint>& subModelNumbers() const
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
  const SubstitutionProcess& process_;

  /*
   * Build rest of the tree under given father (event on father will
   * be set in this method)
   *
   */

  void _build_following_path(std::shared_ptr<ProcessComputationNode> father, const ModelPath& path);

  void _build_following_scenario(std::shared_ptr<ProcessComputationNode> father, const ModelScenario& scenario,  std::map<std::shared_ptr<MixedTransitionModel>, uint>& mMrca);

public:
  /*
   * @brief construction of a ProcessComputationTree from a SubstitutionProcess
   *
   * @param process the SubstitutionProcess
   *
   */

  ProcessComputationTree(const SubstitutionProcess& process);
};
} // end of namespace bpp.

#endif//_COMPUTINGTREE_H_
