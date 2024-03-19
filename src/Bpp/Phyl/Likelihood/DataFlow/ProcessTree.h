// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_PROCESSTREE_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_PROCESSTREE_H

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>
#include <Bpp/Numeric/ParametrizableCollection.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Parametrizable.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/ProcessComputationTree.h>
#include <Bpp/Phyl/Likelihood/SubstitutionProcess.h>

#include "Definitions.h"
#include "Model.h"
#include "Parameter.h"

// From the stl:
#include <string>

namespace bpp
{
class CollectionNodes;

/**
 * @brief Helper: create a map with mutable dataflow nodes for each
 *  branch of the tree.
 *  The map is indexed by branch ids.
 */

using DAGindexes = std::vector<uint>;
using Speciesindex = uint;


// Branch specific DataFlow objects
class ProcessEdge
{
private:
  /**
   * @brief the index of the species in the phyloTree matching this node.
   */
  const Speciesindex speciesIndex_;

  /**
   * @brief Model & BrLen, = 0 if not supporting a model
   */
  std::shared_ptr<ConfiguredParameter> brlen_;

  std::shared_ptr<ConfiguredModel> model_;

  /**
   * @brief Optional number of submodels, in case model_ is mixed
   * and a submodel is used.
   */

  // Not just "size_t nMod_" because dataflow dependency is needed
  // for createMatrix for TransitionMatrixFromModel

  std::shared_ptr<NumericConstant<size_t>> nMod_;

  ValueRef<Eigen::MatrixXd> transitionMatrix_;

  /**
   * @brief Probablity of the edge, used in case of mixture models.
   */
  ValueRef<double> brprob_;

public:
  /**
   * @brief Construction with model and brlen.
   */
  ProcessEdge(uint speciesIndex,
      std::shared_ptr<ConfiguredParameter> brlen,
      std::shared_ptr<ConfiguredModel> model,
      std::shared_ptr<NumericConstant<size_t>> nMod = 0) : speciesIndex_(speciesIndex), brlen_(brlen), model_(model), nMod_(nMod), transitionMatrix_(0), brprob_(0){}

  /**
   * @brief Construction with probability ref from Mixture model
   */
  ProcessEdge(uint speciesIndex,
      ValueRef<double> brprob) : speciesIndex_(speciesIndex), brlen_(0), model_(0), nMod_(0), transitionMatrix_(0), brprob_(brprob){}

  /**
   * @brief Copy construction
   */
  ProcessEdge(const ProcessEdge& edge) : speciesIndex_(edge.speciesIndex_), brlen_(edge.brlen_), model_(edge.model_), nMod_(edge.nMod_), transitionMatrix_(edge.transitionMatrix_), brprob_(edge.brprob_){}

  std::shared_ptr<ConfiguredModel> getModel()
  {
    return model_;
  }

  ConfiguredModel& model()
  {
    return *model_;
  }

  // void setModel(std::shared_ptr<ConfiguredModel> model)
  // {
  //   model_=model;
  // }

  std::shared_ptr<ConfiguredParameter> getBrLen()
  {
    return brlen_;
  }

  void setBrLen(std::shared_ptr<ConfiguredParameter> brlen)
  {
    brlen_ = brlen;
  }

  void setTransitionMatrix(ValueRef<Eigen::MatrixXd> transitionMatrix)
  {
    transitionMatrix_ = transitionMatrix;
  }

  ValueRef<Eigen::MatrixXd> getTransitionMatrix()
  {
    return transitionMatrix_;
  }

  ValueRef<double> getProba()
  {
    return brprob_;
  }

  // void setNMod(std::shared_ptr<NumericConstant<size_t>> nMod)
  // {
  //   nMod_=nMod;
  // }

  std::shared_ptr<NumericConstant<size_t>> getNMod()
  {
    return nMod_;
  }

  uint getSpeciesIndex() const
  {
    return speciesIndex_;
  }
};

typedef ProcessComputationNode ProcessNode;

using ProcessEdgePtr = std::shared_ptr<ProcessEdge>;
using ProcessNodePtr = std::shared_ptr<ProcessNode>;

class ProcessTree : public AssociationTreeGlobalGraphObserver<ProcessNode, ProcessEdge>
{
  Context context_;

  // CollectionNodes* nodes_;

public:
  /**
   * @brief Build a ProcessTree with same topology as a given
   * ParametrizablePhyloTree, and new ConfiguredParameter BrLen.
   * are created
   */
  ProcessTree(Context& context,
      const ParametrizablePhyloTree& tree);

  /**
   * @brief Copy a ProcessTree with all BrLen multiplied by a given rate DF double    *
   */
  ProcessTree(const ProcessTree& tree,
      ValueRef<double> rate);

  /**
   * @brief Build a ProcessTree with same topology as a given
   * ParametrizablePhyloTree. ConfiguredParameter objects are
   * linked to those of a ParameterList, with names BrLen_"suff".
   */
  ProcessTree(Context& context,
      const ParametrizablePhyloTree& tree,
      const ParameterList& parList,
      const std::string& suff);

  /**
   * @brief Build a ProcessTree following a given
   * ProcessComputationTree. So the resulting topology may be
   * different from the given ParametrizablePhyloTree.
   *
   * ConfiguredModels are in a ParametrizableCollection, and BrLen
   * Parameters are in a ProcessTree with same topology as in the
   * SubstitutionProcess.
   */
  ProcessTree(const ProcessComputationTree& tree,
      ParametrizableCollection<ConfiguredModel>& modelColl,
      const ProcessTree& phyloTree);


  ProcessTree* clone() const override
  {
    throw Exception("ProcessTree::clone should not be called.");
  }

  ProcessTree(const ProcessTree& pTree) :
    AssociationTreeGlobalGraphObserver<ProcessNode, ProcessEdge>(pTree.getGraph())
  {
    throw Exception("ProcessTree::ProcessTree should not be called.");
  }

  ProcessTree& operator=(const ProcessTree& pTree)
  {
    throw Exception("ProcessTree::operator= should not be called.");
    // AssociationTreeGlobalGraphObserver<ProcessNode,Value<T>>::operator=(pTree);
    return *this;
  }

  /*
   * @brief Get the edges indexes of the DAG that correspond to
   * the species Index (of the Process tree).
   */

  DAGindexes getDAGEdgesIndexes(const Speciesindex speciesIndex) const;

  /*
   * For inclusion in ParametrizableCollection. Not used
   *
   */
  const ParameterList getParameters() {return ParameterList();}

  bool matchParametersValues(ParameterList&) {return true;}

  /*
   * Static construction methods.
   *
   */

  /**
   * @briefCreate a Process Tree following a Substitution Process. Tree
   * Node parameters are got from ConfiguredParameters PREVIOUSLY
   * built and stored in a ParameterList.
   */
  static std::shared_ptr<ProcessTree> makeProcessTree(
      Context& context,
      std::shared_ptr<const SubstitutionProcessInterface> process,
      ParameterList& parList,
      const std::string& suff = "");

  /**
   * @brief Create a Process Tree following a Substitution Process in a Collection. Tree
   * Node parameters are got from ConfiguredParameters PREVIOUSLY
   * built and stored in a ParameterList.
   */
  static std::shared_ptr<ProcessTree> makeProcessTree(CollectionNodes& collection, size_t pNum);
};

/**
 * @brief Make a Collection of ConfiguredModel, from the models
 * described in the SubstitutionProcess, and sharing
 * ConfiguredParameters from a ParameterList.
 *
 * Paremeter names from the ParameterList may have "_suff"
 * terminations, with suff matching the numbers of the models in the
 * process.
 *
 */
inline ParametrizableCollection<ConfiguredModel> makeConfiguredModelCollection(
    Context& context,
    const SubstitutionProcessInterface& process,
    ParameterList& parList)
{
  ParametrizableCollection<ConfiguredModel> modelColl;

  // Build the ConfiguredModels from the BranchModels
  auto vnMod = process.getModelNumbers();

  for (auto nMod:vnMod)
  {
    auto mod = process.getModel(nMod);

    modelColl.addObject(ConfiguredParametrizable::createConfigured<BranchModelInterface, ConfiguredModel>(context, *mod, parList, (nMod == 1 ? "" : "_" + TextTools::toString(nMod))), nMod); // suffix "_1" will be added if necessary
  }

  return modelColl;
}
} // end of namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_PROCESSTREE_H
