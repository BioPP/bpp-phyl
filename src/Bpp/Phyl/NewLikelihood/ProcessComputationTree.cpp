//
// File: ProcessComputationTree.cpp
// Created by: Laurent Guéguen
// Created on: mardi 9 juillet 2013, à 15h 37
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

#include "ProcessComputationTree.h"

using namespace bpp;
using namespace std;

ProcessComputationTree::ProcessComputationTree(const SubstitutionProcess& process) :
  BaseTree(0)
{
  const ParametrizablePhyloTree& ptree = process.getParametrizablePhyloTree();
  // if no model scenario, copy the basic tree
  if (true)//!process.hasModelScenario())
  {
    setGraph(ptree.getGraph());
    
    auto itN = ptree.allNodesIterator();
    for (itN->start();!itN->end();itN->next())
    {
      auto index=ptree.getNodeIndex(*(*itN));
      auto nnode=std::make_shared<ProcessComputationNode>(*(*(*itN)), index);
      nnode->setProperty("event",NodeEvent::speciationEvent);
      associateNode(nnode,ptree.getNodeGraphid(*(*itN)));
      setNodeIndex(nnode,index);
    }

    auto itE = ptree.allEdgesIterator();
    for (itE->start();!itE->end();itE->next())
    {
      auto index=ptree.getEdgeIndex(*(*itE));
      auto model=process.getModelForNode(index);
      
      auto nedge=std::make_shared<ProcessComputationEdge>(model, index);
      associateEdge(nedge,ptree.getEdgeGraphid(*(*itE)));
      setEdgeIndex(nedge,index);
    }

  }
  
  // const auto scenariModelScenario& scenario,
                                 
  // for (size_t i=0;i<modelSet->getNumberOfHyperNodes();i++){
  //   mvTreeLikelihoods_[tree.getRootId()].push_back(new RNonHomogeneousMixedTreeLikelihood(tree, modelSet, modelSet->getHyperNode(i), upperNode_, rDist, false, usePatterns));
  // }

  // std::vector<int> vDesc; // vector of the explorated descendents

  // int desc;
  // vector<int> vn;
  
  // size_t nbpath = scenario.getNumberOfModelPaths();

  // vDesc.push_back(ptree.getRootId()); //

  // while (vDesc.size() != 0)  {

  //   desc = vDesc.back();
  //   vDesc.pop_back();

  //   vector<int> vExpMod; // vector of the ids of the MixedModels which
  //                        // nodes are not in only one subtree under desc

  //   vector<int> vson = ptree.getSonsId(desc);
  //   std::map<int, vector<int> > mdesc; // map of the subtree nodes for
  //                                      // each son of desc
  //   for (size_t i = 0; i < vson.size(); i++)
  //   {
  //     std::vector<int> vi;
  //     TreeTools::getNodesId(tree, vson[i], vi);
  //     mdesc[vson[i]] = vi;
  //   }

  //   for (size_t i = 0; i < nbmodels; i++)
  //   {
  //     const MixedSubstitutionModelSet::HyperNode::Node& node = hyperNode_.getNode(i);
      
  //     if (node.size()>1)
  //     {
  //       vn = modelSet_->getNodesWithModel(i); // tree nodes associated to model

  //       /* Check if the vn members are in the same subtree */
  //       size_t flag = 0; // count of the subtrees that have vn members
  //       std::map<int, vector<int> >::iterator it;
  //       for (it = mdesc.begin(); it != mdesc.end(); it++)
  //       {
  //         for (size_t j = 0; j < it->second.size(); j++)
  //         {
  //           if (it->second[j] != it->first)
  //           {
  //             if (find(vn.begin(), vn.end(), it->second[j]) != vn.end())
  //             {
  //               flag += (find(vn.begin(), vn.end(), it->first) != vn.end()) ? 2 : 1; // check if the son
  //               // has this model too
  //               break;
  //             }
  //           }
  //           else if (find(vn.begin(), vn.end(), it->first) != vn.end())
  //             flag++;
  //         }
  //         if (flag >= 2)
  //           break;
  //       }
  //       if (flag >= 2)
  //         vExpMod.push_back(static_cast<int>(i));  // mixed model that must be expanded
  //     }
  //   }

  //   if (vExpMod.size() != 0)
  //   {
  //     std::map<int, int> mapmodels;
  //     size_t ttmodels = 1;
  //     for (vector<int>::iterator it = vExpMod.begin(); it != vExpMod.end(); it++)
  //     {
  //       mapmodels[*it] = static_cast<int>(hyperNode_.getNode(static_cast<size_t>(*it)).size());
  //       ttmodels *= static_cast<size_t>(mapmodels[*it]);
  //     }

  //     for (size_t i = 0; i < ttmodels; i++)
  //     {
  //       int s = static_cast<int>(i);
  //       MixedSubstitutionModelSet::HyperNode hn(hyperNode_);
        
  //       for (size_t j = 0; j < nbmodels; j++)
  //       {
  //         if ((hyperNode_.getNode(j).size() >= 1) && find(vExpMod.begin(), vExpMod.end(), static_cast<int>(j)) != vExpMod.end())
  //         {
  //           hn.setModel(j, Vuint(1, hyperNode_.getNode(j)[static_cast<size_t>(s % mapmodels[static_cast<int>(j)])]));
  //           s /= mapmodels[static_cast<int>(j)];
  //         }
  //       }
  //       hn.setProbability((dynamic_cast<MixedSubstitutionModelSet*>(modelSet_))->getHyperNodeProbability(hn));
  //       RNonHomogeneousMixedTreeLikelihood* pr;

  //       if (pdata)
  //         pr = new RNonHomogeneousMixedTreeLikelihood(tree, *pdata, dynamic_cast<MixedSubstitutionModelSet*>(modelSet_), hn, desc, rateDistribution_, false, usePatterns);
  //       else
  //         pr = new RNonHomogeneousMixedTreeLikelihood(tree, dynamic_cast<MixedSubstitutionModelSet*>(modelSet_), hn, desc, rateDistribution_, false, usePatterns);
  //       pr->resetParameters_();
  //       mvTreeLikelihoods_[desc].push_back(pr);
  //     }
  //   }
  //   else
  //     for (size_t i = 0; i < vson.size(); i++)
  //     {
  //       vDesc.push_back(vson[i]);
  //     }
  // }

}

