//
// File: PhyloTreeTools.cpp
// Created by: Julien Dutheil
// Created on: Wed Aug  6 13:45:28 2003
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

#include "PhyloTreeTools.h"

#include "BipartitionTools.h"
#include "../Model/Nucleotide/JCnuc.h"
#include "../Distance/DistanceEstimation.h"
#include "../Distance/BioNJ.h"
#include "../Parsimony/DRTreeParsimonyScore.h"
#include "../OptimizationTools.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/BppString.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <sstream>

using namespace std;

/******************************************************************************/

// int PhyloTreeTools::getLeafId(const PhyloTree& tree, PhyloTree::NodeIndex nodeIndex, const std::string& name)
// throw (NodeNotFoundException)
// {
//   int* id = NULL;
//   searchLeaf(tree, nodeId, name, id);
//   if (id == NULL)
//     throw PhyloNodeNotFoundException("PhyloTreeTools::getLeafId().", name);
//   else
//   {
//     int i = *id;
//     delete id;
//     return i;
//   }
// }
// 
// void PhyloTreeTools::searchLeaf(const PhyloTree& tree, PhyloTree::NodeIndex nodeIndex, const string& name, int*& id)
// throw (NodeNotFoundException)
// {
//   if (tree.isLeaf(nodeId))
//   {
//     if (tree.getNodeName(nodeId) == name)
//     {
//       id = new int(nodeId);
//       return;
//     }
//   }
//   vector<int> sons;
//   for (size_t i = 0; i < sons.size(); i++)
//   {
//     searchLeaf(tree, nodeId, name, id);
//   }
// }
// 


// size_t PhyloTreeTools::getDepth(const PhyloTree& tree, PhyloTree::NodeIndex nodeIndex) throw (NodeNotFoundException)
// {
//   //TODO must check if node exists
// //   if (!tree.hasNode(nodeId))
// //     throw PhyloNodeNotFoundException("PhyloTreeTools::getDepth", nodeId);
//   size_t d = 0;
//   vector<int> sons = tree.getSons(nodeId);
//   for (size_t i = 0; i < sons.size(); i++)
//   {
//     size_t c = getDepth(tree, sons[i]) + 1;
//     if (c > d)
//       d = c;
//   }
//   return d;
// }



// size_t PhyloTreeTools::getDepths(const PhyloTree& tree, PhyloTree::NodeIndex nodeIndex, map<int, size_t>& depths) throw (NodeNotFoundException)
// {
// //   if (!tree.hasNode(nodeId))
// //     throw PhyloNodeNotFoundException("PhyloTreeTools::getDepth", nodeId);
//   size_t d = 0;
//   vector<int> sons = tree.getSonsId(nodeId);
//   for (size_t i = 0; i < sons.size(); i++)
//   {
//     size_t c = getDepths(tree, sons[i], depths) + 1;
//     if (c > d)
//       d = c;
//   }
//   depths[nodeId] = d;
//   return d;
// }



double PhyloTreeTools::getHeight(const PhyloTree& tree, const shared_ptr<PhyloNode> node) 
{
  double d = 0;

  vector<shared_ptr<PhyloBranch> > edges = tree.getOutgoingEdges(node);
  for (size_t i = 0; i < edges.size(); i++)
  {
    double dist = 0;
    if (edges[i]->hasLength())
      dist = edges[i]->getLength();
    else
      throw PhyloBranchPException("Branch without length.", edges[i].get());

    double c = getHeight(tree, get<1>(tree.getNodes(edges[i]))) + dist;
    if (c > d)
      d = c;
  }

  return d;
}



// double PhyloTreeTools::getHeights(const PhyloTree& tree, PhyloTree::NodeIndex nodeIndex, map<int, double>& heights) throw (NodeNotFoundException, PhyloNodeException)
// {
//   if (!tree.hasNode(nodeId))
//     throw PhyloNodeNotFoundException("PhyloTreeTools::getHeight", nodeId);
//   double d = 0;
//   vector<int> sons = tree.getSonsId(nodeId);
//   for (size_t i = 0; i < sons.size(); i++)
//   {
//     double dist = 0;
//     if (tree.hasDistanceToFather(sons[i]))
//       dist = tree.getDistanceToFather(sons[i]);
//     else
//       throw PhyloNodeException("Node without branch length.", sons[i]);
//     double c = getHeights(tree, sons[i], heights) + dist;
//     if (c > d)
//       d = c;
//   }
//   heights[nodeId] = d;
//   return d;
// }



// // double PhyloTreeTools::getDistanceBetweenAnyTwoNodes(const PhyloTree& tree, PhyloTree::NodeIndex nodeIndex1, PhyloTree::NodeIndex nodeIndex2)
// // {
// //   std::vector<PhyloEdge*> edgePath = tree->getEdgePathBetweenTwoNodes(nodeId1,nodeId2);
// //   double d = 0;
// //   for (std::vector<PhyloEdge*>::iterator currBranch = edgePath.begin(); currBranch != edgePath.end(); currBranch++)
// //   {
// //     d += (*currBranch)->getLength();
// //   }
// //   return d;
// // }



// Vdouble PhyloTreeTools::getBranchLengths(const PhyloTree& tree, PhyloTree::NodeIndex nodeIndex) throw (NodeNotFoundException, PhyloNodeException)
// {
//   if (!tree.hasNode(nodeId))
//     throw PhyloNodeNotFoundException("PhyloTreeTools::getBranchLengths", nodeId);
//   Vdouble brLen(1);
//   if (tree.hasDistanceToFather(nodeId))
//     brLen[0] = tree.getDistanceToFather(nodeId);
//   else
//     throw PhyloNodeException("PhyloTreeTools::getbranchLengths(). No branch length.", nodeId);
//   vector<int> sons = tree.getSonsId(nodeId);
//   for (size_t i = 0; i < sons.size(); i++)
//   {
//     Vdouble sonBrLen = getBranchLengths(tree, sons[i]);
//     for (size_t j = 0; j < sonBrLen.size(); j++)
//     {
//       brLen.push_back(sonBrLen[j]);
//     }
//   }
//   return brLen;
// }



// double PhyloTreeTools::getTotalLength(const PhyloTree& tree, PhyloTree::NodeIndex nodeIndex, bool includeAncestor) throw (NodeNotFoundException, PhyloNodeException)
// {
//   if (!tree.hasNode(nodeId))
//     throw PhyloNodeNotFoundException("PhyloTreeTools::getTotalLength", nodeId);
//   if (includeAncestor && !tree.hasDistanceToFather(nodeId))
//     throw PhyloNodeException("PhyloTreeTools::getTotalLength(). No branch length.", nodeId);
//   double length = includeAncestor ? tree.getDistanceToFather(nodeId) : 0;
//   vector<int> sons = tree.getSonsId(nodeId);
//   for (size_t i = 0; i < sons.size(); i++)
//   {
//     length += getTotalLength(tree, sons[i], true);
//   }
//   return length;
// }



// void PhyloTreeTools::setBranchLengths(PhyloTree& tree, shared_ptr<PhyloNode>  node, double brLen) throw (NodeNotFoundException)
// {
//   vector<PhyloBranch*> branches = tree->getSubtreeEdges(node);
//   for (vector<PhyloBranch*>::iterator currBranch = branches.begin(); currBranch != branches.end(); currBranch++)
//   {
//     (*currBranch)->setLenght(brLen);
//   }
// }



// void PhyloTreeTools::setVoidBranchLengths(PhyloTree& tree, shared_ptr<PhyloNode>  node, double brLen) throw (NodeNotFoundException)
// {
//  vector<PhyloBranch*> branches = tree->getSubtreeEdges(node);
//   for (vector<PhyloBranch*>::iterator currBranch = branches.begin(); currBranch != branches.end(); currBranch++)
//   {
//     if(!(*currBranch)->hasLength())
//       (*currBranch)->setLenght(brLen);
//   }
// }



size_t PhyloTreeTools::initBranchLengthsGrafen(PhyloTree& tree, shared_ptr<PhyloNode> node)
{
  vector<shared_ptr<PhyloNode> > sons = tree.getSons(node);
  vector<size_t> h(sons.size());
  for (size_t i = 0; i < sons.size(); i++)
  {
    h[i] = initBranchLengthsGrafen(tree, sons[i]);
  }
  size_t thish = sons.size() == 0 ? 0 : VectorTools::sum<size_t>(h) + sons.size() - 1;
  for (size_t i = 0; i < sons.size(); i++)
  {
    tree.getEdgeToFather(sons[i])->setLength((double)(thish - h[i]));
  }
  return thish;
}


void PhyloTreeTools::initBranchLengthsGrafen(PhyloTree& tree)
{
  initBranchLengthsGrafen(tree, tree.getRoot());
}

void PhyloTreeTools::computeBranchLengthsGrafen(
  PhyloTree& tree,
  shared_ptr<PhyloNode> node,
  double power,
  double total,
  double& height,
  double& heightRaised)
{
  vector<shared_ptr<PhyloNode> > sons = tree.getSons(node);
  vector<double> hr(sons.size());
  height = 0;
  for (size_t i = 0; i < sons.size(); i++)
  {
    shared_ptr<PhyloBranch> branch=tree.getEdgeToFather(sons[i]);
    
    if (branch->hasLength())
    {
      double h;
      computeBranchLengthsGrafen(tree, sons[i], power, total, h, hr[i]);
      double d = h + branch->getLength();
      if (d > height)
        height = d;
    }
    else
      throw PhyloBranchPException ("PhyloTreeTools::computeBranchLengthsGrafen. Branch length lacking.", branch.get());
  }
  heightRaised = pow(height / total, power) * total;
  for (size_t i = 0; i < sons.size(); i++)
  {
    tree.getEdgeToFather(sons[i])->setLength(heightRaised - hr[i]);
  }
}


void PhyloTreeTools::computeBranchLengthsGrafen(PhyloTree& tree, double power, bool init)
{
  shared_ptr<PhyloNode>  root = tree.getRoot();
  if (init)
  {
    initBranchLengthsGrafen(tree);
  }
  // Scale by total heigth:
  double totalHeight = getHeight(tree, root);
  double h, hr;
  computeBranchLengthsGrafen(tree, root, power, totalHeight, h, hr);
}

double PhyloTreeTools::convertToClockTree(PhyloTree& tree, shared_ptr<PhyloNode> node)
{
  vector<shared_ptr<PhyloNode> > sons = tree.getSons(node);

  vector<double> h(sons.size());
  // We compute the mean height:
  double l = 0;
  double maxh = -1.;
  for (size_t i = 0; i < sons.size(); i++)
  {
    shared_ptr<PhyloBranch> branch=tree.getEdgeToFather(sons[i]);

    if (branch->hasLength())
    {
      h[i] = convertToClockTree(tree, sons[i]);
      if (h[i] > maxh)
        maxh = h[i];
      l += h[i] + branch->getLength();
    }
    else
      throw PhyloBranchPException ("PhyloTreeTools::convertToClockTree. Branch length lacking.", branch.get());
  }
  if (sons.size() > 0)
    l /= (double)sons.size();
  if (l < maxh)
    l = maxh;
  for (size_t i = 0; i < sons.size(); i++)
  {
    tree.getEdgeToFather(sons[i])->setLength(l - h[i]);
  }
  return l;
}




double PhyloTreeTools::convertToClockTree2(PhyloTree& tree, shared_ptr<PhyloNode> node)
{
  vector<shared_ptr<PhyloNode> > sons = tree.getSons(node);
  vector<double> h(sons.size());
  // We compute the mean height:
  double l = 0;
  double maxh = -1.;
  for (size_t i = 0; i < sons.size(); i++)
  {
    shared_ptr<PhyloBranch> branch=tree.getEdgeToFather(sons[i]);

    if (branch->hasLength())
    {
      h[i] = convertToClockTree2(tree, sons[i]);
      if (h[i] > maxh)
        maxh = h[i];
      l += h[i] + branch->getLength();
    }
    else
      throw PhyloBranchPException("PhyloTreeTools::convertToClockTree2. Branch length lacking.", branch.get());
  }
  if (sons.size() > 0)
    l /= (double)sons.size();
  for (size_t i = 0; i < sons.size(); i++)
  {
    tree.scaleTree(sons[i], h[i] > 0 ? l / h[i] : 0);
  }
  return l;
}

// // 
// // 
// // 
// // DistanceMatrix* PhyloTreeTools::getDistanceMatrix(const PhyloTree& tree)
// // {
// //   vector<string> names = tree.getLeavesNames();
// //   DistanceMatrix* mat = new DistanceMatrix(names);
// //   for (size_t i = 0; i < names.size(); i++)
// //   {
// //     (*mat)(i, i) = 0;
// //     for (size_t j = 0; j < i; j++)
// //     {
// //       (*mat)(i, j) = (*mat)(j, i) = getDistanceBetweenAnyTwoNodes(tree, tree.getLeafId(names[i]), tree.getLeafId(names[j]));
// //     }
// //   }
// //   return mat;
// // }
// // 
// // 
// // 
// // int PhyloTreeTools::getMPNUId(const PhyloTree& tree, int id)
// // {
// //   vector<int> ids = getNodesId(tree, id);
// //   sort(ids.begin(), ids.end());
// //   // Check if some id is "missing" in the subtree:
// //   for (size_t i = 0; i < ids.size(); i++)
// //   {
// //     if (ids[i] != (int)i)
// //       return (int)i;
// //   }
// //   // Well, all ids are from 0 to n, so send n+1:
// //   return (int)ids.size();
// // }
// // 
// // VectorSiteContainer* PhyloTreeTools::MRPEncode(const vector<Tree*>& vecTr)
// // {
// //   vector<BipartitionList*> vecBipL;
// //   for (size_t i = 0; i < vecTr.size(); i++)
// //   {
// //     vecBipL.push_back(new BipartitionList(*vecTr[i]));
// //   }
// // 
// //   VectorSiteContainer* cont = BipartitionTools::MRPEncode(vecBipL);
// // 
// //   for (size_t i = 0; i < vecTr.size(); i++)
// //   {
// //     delete vecBipL[i];
// //   }
// // 
// //   return cont;
// // }
// // 
// // 
// // 
// // VectorSiteContainer* PhyloTreeTools::MRPEncodeMultilabel(const vector<Tree*>& vecTr)
// // {
// //     vector<BipartitionList*> vecBipL;
// //     for (size_t i = 0; i < vecTr.size(); i++)
// //     {
// //         vecBipL.push_back(new BipartitionList(*vecTr[i]));
// //     }
// //     
// //     VectorSiteContainer* cont = BipartitionTools::MRPEncodeMultilabel(vecBipL);
// //     
// //     for (size_t i = 0; i < vecTr.size(); i++)
// //     {
// //         delete vecBipL[i];
// //     }
// //     
// //     return cont;
// // }
// // 
// // 
// // 
// // bool PhyloTreeTools::haveSameTopology(const PhyloTree& tr1, const PhyloTree& tr2)
// // {
// //   size_t jj, nbbip;
// //   BipartitionList* bipL1, * bipL2;
// //   vector<size_t> size1, size2;
// // 
// //   /* compare sets of leaves */
// //   if (!VectorTools::haveSameElements(tr1.getLeavesNames(), tr2.getLeavesNames()))
// //     return false;
// // 
// //   /* construct bipartitions */
// //   bipL1 = new BipartitionList(tr1, true);
// //   bipL1->removeTrivialBipartitions();
// //   bipL1->removeRedundantBipartitions();
// //   bipL1->sortByPartitionSize();
// //   bipL2 = new BipartitionList(tr2, true);
// //   bipL2->removeTrivialBipartitions();
// //   bipL2->removeRedundantBipartitions();
// //   bipL2->sortByPartitionSize();
// // 
// //   /* compare numbers of bipartitions */
// //   if (bipL1->getNumberOfBipartitions() != bipL2->getNumberOfBipartitions())
// //     return false;
// //   nbbip = bipL1->getNumberOfBipartitions();
// // 
// //   /* compare partition sizes */
// //   for (size_t i = 0; i < nbbip; i++)
// //   {
// //     size1.push_back(bipL1->getPartitionSize(i));
// //     size2.push_back(bipL1->getPartitionSize(i));
// //     if (size1[i] != size2[i])
// //       return false;
// //   }
// // 
// //   /* compare bipartitions */
// //   for (size_t i = 0; i < nbbip; i++)
// //   {
// //     for (jj = 0; jj < nbbip; jj++)
// //     {
// //       if (size1[i] == size2[jj] && BipartitionTools::areIdentical(*bipL1, i, *bipL2, jj))
// //         break;
// //     }
// //     if (jj == nbbip)
// //       return false;
// //   }
// // 
// //   return true;
// // }
// // 
// // 
// // 
// // int PhyloTreeTools::robinsonFouldsDistance(const PhyloTree& tr1, const PhyloTree& tr2, bool checkNames, int* missing_in_tr2, int* missing_in_tr1) throw (Exception)
// // {
// //   BipartitionList* bipL1, * bipL2;
// //   size_t i, j;
// //   vector<size_t> size1, size2;
// //   vector<bool> bipOK2;
// // 
// // 
// //   if (checkNames && !VectorTools::haveSameElements(tr1.getLeavesNames(), tr2.getLeavesNames()))
// //     throw Exception("Distinct leaf sets between trees ");
// // 
// //   /* prepare things */
// //   int missing1 = 0;
// //   int missing2 = 0;
// // 
// //   bipL1 = new BipartitionList(tr1, true);
// //   bipL1->removeTrivialBipartitions();
// //   bipL1->sortByPartitionSize();
// //   bipL2 = new BipartitionList(tr2, true);
// //   bipL2->removeTrivialBipartitions();
// //   bipL2->sortByPartitionSize();
// // 
// // 
// //   for (i = 0; i < bipL1->getNumberOfBipartitions(); i++)
// //   {
// //     size1.push_back(bipL1->getPartitionSize(i));
// //   }
// //   for (i = 0; i < bipL2->getNumberOfBipartitions(); i++)
// //   {
// //     size2.push_back(bipL2->getPartitionSize(i));
// //   }
// // 
// //   for (i = 0; i < bipL2->getNumberOfBipartitions(); i++)
// //   {
// //     bipOK2.push_back(false);
// //   }
// // 
// //   /* main loops */
// // 
// //   for (i = 0; i < bipL1->getNumberOfBipartitions(); i++)
// //   {
// //     for (j = 0; j < bipL2->getNumberOfBipartitions(); j++)
// //     {
// //       if (bipOK2[j])
// //         continue;
// //       if (size1[i] == size2[j] && BipartitionTools::areIdentical(*bipL1, i, *bipL2, j))
// //       {
// //         bipOK2[j] = true;
// //         break;
// //       }
// //     }
// //     if (j == bipL2->getNumberOfBipartitions())
// //       missing2++;
// //   }
// // 
// //   missing1 = static_cast<int>(bipL2->getNumberOfBipartitions()) - static_cast<int>(bipL1->getNumberOfBipartitions()) + missing2;
// // 
// //   if (missing_in_tr1)
// //     *missing_in_tr1 = missing1;
// //   if (missing_in_tr2)
// //     *missing_in_tr2 = missing2;
// //   return missing1 + missing2;
// // }
// // 
// // 
// // 
// // BipartitionList* PhyloTreeTools::bipartitionOccurrences(const vector<Tree*>& vecTr, vector<size_t>& bipScore)
// // {
// //   vector<BipartitionList*> vecBipL;
// //   BipartitionList* mergedBipL;
// //   vector<size_t> bipSize;
// //   size_t nbBip;
// // 
// //   /*  build and merge bipartitions */
// //   for (size_t i = 0; i < vecTr.size(); i++)
// //   {
// //     vecBipL.push_back(new BipartitionList(*vecTr[i]));
// //   }
// //   mergedBipL = BipartitionTools::mergeBipartitionLists(vecBipL);
// //   for (size_t i = 0; i < vecTr.size(); i++)
// //   {
// //     delete vecBipL[i];
// //   }
// // 
// //   mergedBipL->removeTrivialBipartitions();
// //   nbBip = mergedBipL->getNumberOfBipartitions();
// //   bipScore.clear();
// //   for (size_t i = 0; i < nbBip; i++)
// //   {
// //     bipSize.push_back(mergedBipL->getPartitionSize(i));
// //     bipScore.push_back(1);
// //   }
// // 
// //   /* compare bipartitions */
// //   for (size_t i = nbBip; i > 0; i--)
// //   {
// //     if (bipScore[i - 1] == 0)
// //       continue;
// //     for (size_t j = i - 1; j > 0; j--)
// //     {
// //       if (bipScore[j - 1] && bipSize[i - 1] == bipSize[j - 1] && mergedBipL->areIdentical(i - 1, j - 1))
// //       {
// //         bipScore[i - 1]++;
// //         bipScore[j - 1] = 0;
// //       }
// //     }
// //   }
// // 
// //   /* keep only distinct bipartitions */
// //   for (size_t i = nbBip; i > 0; i--)
// //   {
// //     if (bipScore[i - 1] == 0)
// //     {
// //       bipScore.erase(bipScore.begin() + static_cast<ptrdiff_t>(i - 1));
// //       mergedBipL->deleteBipartition(i - 1);
// //     }
// //   }
// // 
// //   /* add terminal branches */
// //   mergedBipL->addTrivialBipartitions(false);
// //   for (size_t i = 0; i < mergedBipL->getNumberOfElements(); i++)
// //   {
// //     bipScore.push_back(vecTr.size());
// //   }
// // 
// //   return mergedBipL;
// // }
// // 
// // 
// // 
// // PhyloTree<Node>* PhyloTreeTools::thresholdConsensus(const vector<Tree*>& vecTr, double threshold, bool checkNames) throw (Exception)
// // {
// //   vector<size_t> bipScore;
// //   vector<string> tr0leaves;
// //   BipartitionList* bipL;
// //   double score;
// // 
// //   if (vecTr.size() == 0)
// //     throw Exception("PhyloTreeTools::thresholdConsensus. Empty vector passed");
// // 
// //   /* check names */
// //   if (checkNames)
// //   {
// //     tr0leaves = vecTr[0]->getLeavesNames();
// //     for (size_t i = 1; i < vecTr.size(); i++)
// //     {
// //       if (!VectorTools::haveSameElements(vecTr[i]->getLeavesNames(), tr0leaves))
// //         throw Exception("PhyloTreeTools::thresholdConsensus. Distinct leaf sets between trees");
// //     }
// //   }
// // 
// //   bipL = bipartitionOccurrences(vecTr, bipScore);
// // 
// //   for (size_t i = bipL->getNumberOfBipartitions(); i > 0; i--)
// //   {
// //     if (bipL->getPartitionSize(i - 1) == 1)
// //       continue;
// //     score = static_cast<int>(bipScore[i - 1]) / static_cast<double>(vecTr.size());
// //     if (score <= threshold && score != 1.)
// //     {
// //       bipL->deleteBipartition(i - 1);
// //       continue;
// //     }
// //     if (score > 0.5)
// //       continue;
// //     for (size_t j = bipL->getNumberOfBipartitions(); j > i; j--)
// //     {
// //       if (!bipL->areCompatible(i - 1, j - 1))
// //       {
// //         bipL->deleteBipartition(i - 1);
// //         break;
// //       }
// //     }
// //   }
// // 
// //   PhyloTree<Node>* tr = bipL->toTree();
// //   delete bipL;
// //   return tr;
// // }
// // 
// // 
// // 
// // PhyloTree<Node>* PhyloTreeTools::fullyResolvedConsensus(const vector<Tree*>& vecTr, bool checkNames)
// // {
// //   return thresholdConsensus(vecTr, 0., checkNames);
// // }
// // 
// // 
// // 
// // PhyloTree<Node>* PhyloTreeTools::majorityConsensus(const vector<Tree*>& vecTr, bool checkNames)
// // {
// //   return thresholdConsensus(vecTr, 0.5, checkNames);
// // }
// // 
// // 
// // 
// // PhyloTree<Node>* PhyloTreeTools::strictConsensus(const vector<Tree*>& vecTr, bool checkNames)
// // {
// //   return thresholdConsensus(vecTr, 1., checkNames);
// // }
// // 
// // 
// // 
// // Tree* PhyloTreeTools::MRP(const vector<Tree*>& vecTr)
// // {
// //   // matrix representation
// //   VectorSiteContainer* sites = PhyloTreeTools::MRPEncode(vecTr);
// // 
// //   // starting bioNJ tree
// //   const DNA* alphabet = dynamic_cast<const DNA*>(sites->getAlphabet());
// //   JCnuc* jc = new JCnuc(alphabet);
// //   ConstantDistribution* constRate = new ConstantDistribution(1.);
// //   DistanceEstimation distFunc(jc, constRate, sites, 0, true);
// //   BioNJ bionjTreeBuilder(false, false);
// //   bionjTreeBuilder.setDistanceMatrix(*(distFunc.getMatrix()));
// //   bionjTreeBuilder.computeTree();
// //   if (ApplicationTools::message)
// //     ApplicationTools::message->endLine();
// //   PhyloTree<Node>* startTree = new PhyloTree<Node>(*bionjTreeBuilder.getTree());
// // 
// //   // MP optimization
// //   DRTreeParsimonyScore* MPScore = new DRTreeParsimonyScore(*startTree, *sites, false);
// //   MPScore = OptimizationTools::optimizeTreeNNI(MPScore, 0);
// //   delete startTree;
// //   Tree* retTree = new PhyloTree<Node>(MPScore->getTree());
// //   delete MPScore;
// // 
// //   return retTree;
// // }
// // 
// // 
// // 
// // void PhyloTreeTools::computeBootstrapValues(PhyloTree& tree, const vector<Tree*>& vecTr, bool verbose, int format)
// // {
// //   vector<int> index;
// //   BipartitionList bpTree(tree, true, &index);
// //   vector<size_t> occurences;
// //   BipartitionList* bpList = bipartitionOccurrences(vecTr, occurences);
// // 
// //   vector< Number<double> > bootstrapValues(bpTree.getNumberOfBipartitions());
// // 
// //   for (size_t i = 0; i < bpTree.getNumberOfBipartitions(); i++)
// //   {
// //     if (verbose)
// //       ApplicationTools::displayGauge(i, bpTree.getNumberOfBipartitions() - 1, '=');
// //     for (size_t j = 0; j < bpList->getNumberOfBipartitions(); j++)
// //     {
// //       if (BipartitionTools::areIdentical(bpTree, i, *bpList, j))
// //       {
// //         bootstrapValues[i] = format >= 0 ? round(static_cast<double>(occurences[j]) * pow(10., 2 + format) / static_cast<double>(vecTr.size())) / pow(10., format) : static_cast<double>(occurences[j]);
// //         break;
// //       }
// //     }
// //   }
// // 
// //   for (size_t i = 0; i < index.size(); i++)
// //   {
// //     if (!tree.isLeaf(index[i]))
// //       tree.setBranchProperty(index[i], BOOTSTRAP, bootstrapValues[i]);
// //   }
// // 
// //   delete bpList;
// }
// 
// 
// 
// vector<int> PhyloTreeTools::getAncestors(const PhyloTree& tree, PhyloTree::NodeIndex nodeIndex) throw (NodeNotFoundException)
// {
//   vector<int> ids;
//   int currentId = nodeId;
//   while (tree.hasFather(currentId))
//   {
//     currentId = tree.getFatherId(currentId);
//     ids.push_back(currentId);
//   }
//   return ids;
// }
// 
// 
// 
// int PhyloTreeTools::getLastCommonAncestor(const PhyloTree& tree, const vector<int>& nodeIds) throw (NodeNotFoundException, Exception)
// {
//   if (nodeIds.size() == 0)
//     throw Exception("PhyloTreeTools::getLastCommonAncestor(). You must provide at least one node id.");
//   vector< vector<int> > ancestors(nodeIds.size());
//   for (size_t i = 0; i < nodeIds.size(); i++)
//   {
//     ancestors[i] = getAncestors(tree, nodeIds[i]);
//     ancestors[i].insert(ancestors[i].begin(), nodeIds[i]);
//   }
//   int lca = tree.getRootId();
//   size_t count = 1;
//   for ( ; ; )
//   {
//     if (ancestors[0].size() <= count)
//       return lca;
//     int current = ancestors[0][ancestors[0].size() - count - 1];
//     for (size_t i = 1; i < nodeIds.size(); i++)
//     {
//       if (ancestors[i].size() <= count)
//         return lca;
//       if (ancestors[i][ancestors[i].size() - count - 1] != current)
//         return lca;
//     }
//     lca = current;
//     count++;
//   }
//   // This line is never reached!
//   return lca;
// }
// 
// 
// 

void PhyloTreeTools::constrainedMidPointRooting(PhyloTree& tree)
{
  // is the tree rooted?
  if (!tree.isRooted())
    throw Exception("The tree has to be rooted on the branch of interest to determine the midpoint position of the root");

  vector<shared_ptr<PhyloNode> > sons = tree.getSons(tree.getRoot());

  if (sons.size()>2)
    throw Exception("The tree is multifurcated at the root, which is not allowed.");

  double length = 0.;

  // Length of the branch containing the root:
  shared_ptr<PhyloBranch> branch0=tree.getEdgeToFather(sons[0]);
  shared_ptr<PhyloBranch> branch1=tree.getEdgeToFather(sons[1]);

  length = branch0->getLength() + branch1->getLength();
  
  // The fraction of the original branch allowing to split its length and to place the root:
  double x = bestRootPosition_(tree, sons[0], sons[1], length);
  // The new branch lengths are then computed:
  branch0->setLength(length * x);
  branch1->setLength(length * (1 - x));
}



double PhyloTreeTools::bestRootPosition_(const PhyloTree& tree, const shared_ptr<PhyloNode>  node1, const shared_ptr<PhyloNode> node2, double length)
{
  double x;
  Moments_ m1, m2;
  double A, B; // C;
  // The variance is expressed as a degree 2 polynomial : variance(x) = A * x * x + B * x + C
  // The fraction x is then obtained by differentiating this equation.
  m1 = statFromNode_(tree, node1);
  m2 = statFromNode_(tree, node2);
  A = 4 * m1.N * (m2.N * length) * length;
  B = 4 * length * (m2.N * m1.sum - m1.N * m2.sum - length * m1.N * m2.N);
//   C = (m1.N + m2.N) * (m1.squaredSum + m2.squaredSum) + m1.N * length * m2.N * length +
//     2 * m1.N * length * m2.sum - 2 * m2.N * length * m1.sum -
//     (m1.sum + m2.sum) * (m1.sum + m2.sum);

  if (A < 1e-20)
    x = 0.5;
  else
    x = -B / (2 * A);
  if (x < 0)
    x = 0;
  else if (x > 1)
    x = 1;

  return x;
}



PhyloTreeTools::Moments_ PhyloTreeTools::statFromNode_(const PhyloTree& tree, const shared_ptr<PhyloNode> root)
{
  // This function recursively calculates both the sum of the branch lengths and the sum of the squared branch lengths down the node whose ID is rootId.
  // If below a particular node there are N leaves, the branch between this node and its father is taken into account N times in the calculation.
  Moments_ m;
  static Moments_ mtmp;

  if (tree.isLeaf(root))
  {
    m.N = 1;
    m.sum = 0.;
    m.squaredSum = 0.;
  }
  else
  {
    vector<shared_ptr<PhyloNode> > sons = tree.getSons(root);
    for (size_t i = 0; i < sons.size(); i++)
    {
      mtmp = statFromNode_(tree, sons[i]);
      shared_ptr<PhyloBranch> branch=tree.getEdgeToFather(sons[i]);

      double bLength = branch->getLength();
      m.N += mtmp.N;
      m.sum += mtmp.sum + bLength * mtmp.N;
      m.squaredSum += mtmp.squaredSum + 2 * bLength * mtmp.sum + mtmp.N * bLength * bLength;
    }
  }

  return m;

}
// 
// 
// 
// Tree* PhyloTreeTools::MRPMultilabel(const vector<Tree*>& vecTr)
// {
//     // matrix representation
//     VectorSiteContainer* sites = PhyloTreeTools::MRPEncode(vecTr);
//     
//     // starting bioNJ tree
//     const DNA* alphabet = dynamic_cast<const DNA*>(sites->getAlphabet());
//     JCnuc* jc = new JCnuc(alphabet);
//     ConstantDistribution* constRate = new ConstantDistribution(1.);
//     DistanceEstimation distFunc(jc, constRate, sites, 0, true);
//     BioNJ bionjTreeBuilder(false, false);
//     bionjTreeBuilder.setDistanceMatrix(*(distFunc.getMatrix()));
//     bionjTreeBuilder.computeTree();
//     if (ApplicationTools::message)
//         ApplicationTools::message->endLine();
//     PhyloTree<Node>* startTree = new PhyloTree<Node>(*bionjTreeBuilder.getTree());
//     
//     // MP optimization
//     DRTreeParsimonyScore* MPScore = new DRTreeParsimonyScore(*startTree, *sites, false);
//     MPScore = OptimizationTools::optimizeTreeNNI(MPScore, 0);
//     delete startTree;
//     Tree* retTree = new PhyloTree<Node>(MPScore->getTree());
//     delete MPScore;
//     
//     return retTree;
// }*/


