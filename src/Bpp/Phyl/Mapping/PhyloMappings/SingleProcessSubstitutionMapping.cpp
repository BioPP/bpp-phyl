// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SingleProcessSubstitutionMapping.h"

using namespace bpp;
using namespace std;

SingleProcessSubstitutionMapping::SingleProcessSubstitutionMapping(
    shared_ptr<SingleProcessPhyloLikelihood> spp,
    shared_ptr<SubstitutionRegisterInterface> reg,
    std::shared_ptr<const AlphabetIndex2> weights,
    std::shared_ptr<const AlphabetIndex2> distances,
    double threshold,
    bool verbose) :
  AbstractSinglePhyloSubstitutionMapping(spp->tree()->getGraph(), reg, weights, distances),
  pSPP_(spp)
{
  setBranchedModelSet_();

  // assigns edge indexes
  const auto tree = spp->tree();

  unique_ptr<modelTree::EdgeIterator> eIT = allEdgesIterator();

  for ( ; !eIT->end(); eIT->next())
  {
    auto edge1 = tree->getEdgeFromGraphid(getEdgeGraphid(**eIT));
    if (tree->hasEdgeIndex(edge1))
      setEdgeIndex(**eIT, tree->getEdgeIndex(edge1));
  }

  // assigns node indexes
  unique_ptr<modelTree::NodeIterator> nIT = allNodesIterator();

  for ( ; !nIT->end(); nIT->next())
  {
    auto node1 = tree->getNodeFromGraphid(getNodeGraphid(**nIT));
    if (tree->hasNodeIndex(node1))
      setNodeIndex(**nIT, tree->getNodeIndex(node1));
  }
}


void SingleProcessSubstitutionMapping::computeNormalizations(const ParameterList& nullParams,
                                                             short unresolvedOption,
                                                             bool verbose)
{
  matchParametersValues(nullParams);

  factors_ = SubstitutionMappingTools::computeNormalizations(
      getLikelihoodCalculationSingleProcess(),
      shared_from_this(),
      getSubstitutionRegister(),
      getDistances(),
      unresolvedOption,
      verbose);
}

void SingleProcessSubstitutionMapping::setBranchedModelSet_()
{
  const SubstitutionProcessInterface& sp = pSPP_->substitutionProcess();

  vector<size_t> vId = sp.getModelNumbers();

  for (auto id : vId)
  {
    addModel(id, dynamic_cast<const TransitionModelInterface&>(sp.model(id)), sp.getNodesWithModel(id));
  }
}
