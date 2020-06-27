//
// File: SingleProcessSubstitutionMapping.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 7 décembre 2017, à 22h 59
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

#include "SingleProcessSubstitutionMapping.h"

using namespace bpp;
using namespace std;

SingleProcessSubstitutionMapping::SingleProcessSubstitutionMapping(SingleProcessPhyloLikelihood& spp, SubstitutionRegister& reg, std::shared_ptr<const AlphabetIndex2> weights, std::shared_ptr<const AlphabetIndex2> distances, double threshold, bool verbose) :
  AbstractSinglePhyloSubstitutionMapping(spp.getTree().getGraph(), reg, weights, distances),
  pSPP_(&spp)
{
  setBranchedModelSet_();

  // assigns edge indexes
  const auto& tree = spp.getTree();

  unique_ptr<modelTree::EdgeIterator> eIT=allEdgesIterator();

  for (;!eIT->end(); eIT->next())
  {
    auto edge1 = tree.getEdgeFromGraphid(getEdgeGraphid(**eIT));
    if (tree.hasEdgeIndex(edge1))
      setEdgeIndex(**eIT,tree.getEdgeIndex(edge1));
  }
  
  // assigns node indexes
  unique_ptr<modelTree::NodeIterator> nIT=allNodesIterator();

  for (;!nIT->end(); nIT->next())
  {
    auto node1 = tree.getNodeFromGraphid(getNodeGraphid(**nIT));
    if (tree.hasNodeIndex(node1))
      setNodeIndex(**nIT,tree.getNodeIndex(node1));
  }
}



void SingleProcessSubstitutionMapping::computeNormalizations(const ParameterList& nullParams,
                                                             bool verbose)
{
  matchParametersValues(nullParams);

  factors_.reset(SubstitutionMappingTools::computeNormalizations(getLikelihoodCalculationSingleProcess(),
                                                                 this,
                                                                 getRegister(),
                                                                 getDistances(),
                                                                 verbose));
}

// void SingleProcessSubstitutionMapping::computeNormalizationsForASite(
//   size_t site,
//   const ParameterList& nullParams,
//   bool verbose)
// {
//   matchParametersValues(nullParams);

//   factors_.reset(SubstitutionMappingToolsForASite::computeNormalizations(
//                    site,
//                    getLikelihoodCalculationSingleProcess(),
//                    this,
//                    getRegister(),
//                    getDistances(),
//                    verbose));
// }

void SingleProcessSubstitutionMapping::setBranchedModelSet_()
{
  const SubstitutionProcess& sp=pSPP_->getSubstitutionProcess();
  
  vector<size_t> vId=sp.getModelNumbers();
  
  for (auto id:vId)
    addModel(id, dynamic_cast<const TransitionModel&>(*sp.getModel(id)),sp.getNodesWithModel(id));
}

  
