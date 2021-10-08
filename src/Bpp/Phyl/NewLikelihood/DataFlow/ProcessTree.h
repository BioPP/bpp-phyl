//
// File: ProcessTree.h
// Authors:
//   Laurent GuÃ©guen
// Created: jeudi 15 septembre 2016, Ã  06h 40
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_PROCESSTREE_H
#define BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_PROCESSTREE_H

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>
#include <Bpp/Numeric/ParametrizableCollection.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/Parametrizable.h>
#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/NewLikelihood/ProcessComputationTree.h>
#include <Bpp/Phyl/NewLikelihood/SubstitutionProcess.h>

#include "Definitions.h"
#include "Model.h"
#include "Parameter.h"

// From the stl:
#include <string>

namespace bpp
{
class CollectionNodes;

/** Helper: create a map with mutable dataflow nodes for each
 *  branch of the tree.
 *  The map is indexed by branch ids.
 */

// Branch specific DataFlow objects
class ProcessEdge
{
private:
  /*
   * @brief the index of the species in the phyloTree matching this node.
   *
   */

  const uint speciesIndex_;

  /*
   * @brief Model & BrLen, = 0 if not supporting a model
   *
   *
   */

  std::shared_ptr<ConfiguredParameter> brlen_;

  std::shared_ptr<ConfiguredModel> model_;

  /*
   *@brief Optional number of submodels, in case model_ is mixed
   * and a submodel is used.
   *
   */

  // Not just "size_t nMod_" because dataflow dependency is needed
  // for createMatrix for TransitionMatrixFromModel

  std::shared_ptr<NumericConstant<size_t> > nMod_;

  ValueRef<Eigen::MatrixXd> transitionMatrix_;
  /*
   *@ brief Probablity of the edge, used in case of mixture models.
   *
   */

  ValueRef<double> brprob_;

public:
  /*
   * @brief Construction with model and brlen.
   *
   */

  ProcessEdge(uint speciesIndex,
              std::shared_ptr<ConfiguredParameter> brlen,
              std::shared_ptr<ConfiguredModel> model,
              std::shared_ptr<NumericConstant<size_t> > nMod = 0) : speciesIndex_(speciesIndex), brlen_(brlen), model_(model), nMod_(nMod), transitionMatrix_(0), brprob_(0){}

  /*
   * @brief Construction with probability ref from Mixture model
   *
   */

  ProcessEdge(uint speciesIndex,
              ValueRef<double> brprob) : speciesIndex_(speciesIndex), brlen_(0), model_(0), nMod_(0), transitionMatrix_(0), brprob_(brprob){}

  /*
   * @brief Copy construction
   *
   */

  ProcessEdge(const ProcessEdge& edge) : speciesIndex_(edge.speciesIndex_), brlen_(edge.brlen_), model_(edge.model_), nMod_(edge.nMod_), transitionMatrix_(edge.transitionMatrix_), brprob_(edge.brprob_){}

  std::shared_ptr<ConfiguredModel> getModel()
  {
    return model_;
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

  std::shared_ptr<NumericConstant<size_t> > getNMod()
  {
    return nMod_;
  }

  uint getSpeciesIndex() const
  {
    return speciesIndex_;
  }
};

typedef ProcessComputationNode ProcessNode;

using ProcessEdgeRef = std::shared_ptr<ProcessEdge>;
using ProcessNodeRef = std::shared_ptr<ProcessNode>;

class ProcessTree : public AssociationTreeGlobalGraphObserver<ProcessNode, ProcessEdge>
{
  Context context_;

  // CollectionNodes* nodes_;

public:
  /*
   * @brief Build a ProcessTree with same topology as a given
   * ParametrizablePhyloTree, and new ConfiguredParameter BrLen.
   * are created
   *
   */

  ProcessTree(Context& context,
              const ParametrizablePhyloTree& tree);

  /*
   * @brief Copy a ProcessTree with all BrLen multiplied by a given rate DF double    *
   */

  ProcessTree(const ProcessTree& tree,
              ValueRef<double> rate);

  /*
   * @brief Build a ProcessTree with same topology as a given
   * ParametrizablePhyloTree. ConfiguredParameter objects are
   * linked to those of a ParameterList, with names BrLen_"suff".
   *
   */

  ProcessTree(Context& context,
              const ParametrizablePhyloTree& tree,
              const ParameterList& parList,
              const std::string& suff);

  /*
   * @brief Build a ProcessTree following a given
   * ProcessComputationTree. So the resulting topology may be
   * different from the given ParametrizablePhyloTree.
   *
   * ConfiguredModels are in a ParametrizableCollection, and BrLen
   * Parameters are in a ProcessTree with same topology as in the
   * SubstitutionProcess.
   *
   */

  ProcessTree(const ProcessComputationTree& tree,
              ParametrizableCollection<ConfiguredModel>& modelColl,
              const ProcessTree& phyloTree);


  ProcessTree* clone() const
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
   * For inclusion in ParametrizableCollection. Not used
   *
   */
  const ParameterList getParameters() {return ParameterList();}

  bool matchParametersValues(ParameterList&) {return true;}

  /*
   * Static construction methods.
   *
   */

  /***************************************************/
  /** Create a Process Tree following a Substitution Process. Tree
   * Node parameters are got from ConfiguredParameters PREVIOUSLY
   * built and stored in a ParameterList.
   *
   */

  static std::shared_ptr<ProcessTree> makeProcessTree(Context& context, const SubstitutionProcess& process, ParameterList& parList, const std::string& suff = "");

  /***************************************************/
  /** Create a Process Tree following a Substitution Process in a Collection. Tree
   * Node parameters are got from ConfiguredParameters PREVIOUSLY
   * built and stored in a ParameterList.
   *
   */

  static std::shared_ptr<ProcessTree> makeProcessTree(CollectionNodes& collection, size_t pNum);
};

/*
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
  const SubstitutionProcess& process,
  ParameterList& parList)
{
  ParametrizableCollection<ConfiguredModel> modelColl;

  // Build the ConfiguredModels from the BranchModels
  auto vnMod = process.getModelNumbers();

  for (auto nMod:vnMod)
  {
    auto mod = process.getModel(nMod);

    modelColl.addObject(ConfiguredParametrizable::createConfigured<BranchModel, ConfiguredModel>(context, *mod, parList, (nMod == 1 ? "" : "_" + TextTools::toString(nMod))), nMod); // suffix "_1" will be added if necessary
  }

  return modelColl;
}
} // end of namespace bpp
#endif // BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_PROCESSTREE_H
