// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/BppString.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>

#include "PhyloTreeTools.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <sstream>

using namespace std;

/******************************************************************************/

const string PhyloTreeTools::BOOTSTRAP = "bootstrap";

/******************************************************************************/

std::shared_ptr<PhyloTree> PhyloTreeTools::buildFromTreeTemplate(const TreeTemplate<Node>& treetemp)
{
  auto phyloT = std::make_shared<PhyloTree>(true);
  const Node& root = *treetemp.getRootNode();

  auto rooti = std::make_shared<PhyloNode>(root.hasName() ? root.getName() : "");
  phyloT->createNode(rooti);
  phyloT->setRoot(rooti);
  phyloT->setNodeIndex(rooti, (unsigned int)root.getId());

  auto propi = root.getNodePropertyNames();
  for (const auto& prop:propi)
  {
    rooti->setProperty(prop, *root.getNodeProperty(prop));
  }

  phyloT->addSubTree(rooti, root);

  return phyloT;
}

double PhyloTreeTools::getHeight(const PhyloTree& tree, const std::shared_ptr<PhyloNode> node)
{
  double d = 0;

  vector<shared_ptr<PhyloBranch>> edges = tree.getOutgoingEdges(node);
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


size_t PhyloTreeTools::initBranchLengthsGrafen(PhyloTree& tree, std::shared_ptr<PhyloNode> node)
{
  vector<shared_ptr<PhyloNode>> sons = tree.getSons(node);
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
    std::shared_ptr<PhyloNode> node,
    double power,
    double total,
    double& height,
    double& heightRaised)
{
  vector<shared_ptr<PhyloNode>> sons = tree.getSons(node);
  vector<double> hr(sons.size());
  height = 0;
  for (size_t i = 0; i < sons.size(); i++)
  {
    shared_ptr<PhyloBranch> branch = tree.getEdgeToFather(sons[i]);

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
  heightRaised = std::pow(height / total, power) * total;
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
  // Scale by total height:
  double totalHeight = getHeight(tree, root);
  double h, hr;
  computeBranchLengthsGrafen(tree, root, power, totalHeight, h, hr);
}

double PhyloTreeTools::convertToClockTree(PhyloTree& tree, std::shared_ptr<PhyloNode> node)
{
  vector<shared_ptr<PhyloNode>> sons = tree.getSons(node);

  vector<double> h(sons.size());
  // We compute the mean height:
  double l = 0;
  double maxh = -1.;
  for (size_t i = 0; i < sons.size(); i++)
  {
    shared_ptr<PhyloBranch> branch = tree.getEdgeToFather(sons[i]);

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


double PhyloTreeTools::convertToClockTree2(PhyloTree& tree, std::shared_ptr<PhyloNode> node)
{
  vector<shared_ptr<PhyloNode>> sons = tree.getSons(node);
  vector<double> h(sons.size());
  // We compute the mean height:
  double l = 0;
  double maxh = -1.;
  for (size_t i = 0; i < sons.size(); i++)
  {
    shared_ptr<PhyloBranch> branch = tree.getEdgeToFather(sons[i]);

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


void PhyloTreeTools::constrainedMidPointRooting(PhyloTree& tree)
{
  // is the tree rooted?
  if (!tree.isRooted())
    throw Exception("The tree has to be rooted on the branch of interest to determine the midpoint position of the root");

  vector<shared_ptr<PhyloNode>> sons = tree.getSons(tree.getRoot());

  if (sons.size() > 2)
    throw Exception("The tree is multifurcated at the root, which is not allowed.");

  double length = 0.;

  // Length of the branch containing the root:
  shared_ptr<PhyloBranch> branch0 = tree.getEdgeToFather(sons[0]);
  shared_ptr<PhyloBranch> branch1 = tree.getEdgeToFather(sons[1]);

  length = branch0->getLength() + branch1->getLength();

  // The fraction of the original branch allowing to split its length and to place the root:
  double x = bestRootPosition_(tree, sons[0], sons[1], length);
  // The new branch lengths are then computed:
  branch0->setLength(length * x);
  branch1->setLength(length * (1 - x));
}


double PhyloTreeTools::bestRootPosition_(const PhyloTree& tree, const std::shared_ptr<PhyloNode>  node1, const std::shared_ptr<PhyloNode> node2, double length)
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


PhyloTreeTools::Moments_ PhyloTreeTools::statFromNode_(const PhyloTree& tree, const std::shared_ptr<PhyloNode> root)
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
    vector<shared_ptr<PhyloNode>> sons = tree.getSons(root);
    for (size_t i = 0; i < sons.size(); i++)
    {
      mtmp = statFromNode_(tree, sons[i]);
      shared_ptr<PhyloBranch> branch = tree.getEdgeToFather(sons[i]);

      double bLength = branch->getLength();
      m.N += mtmp.N;
      m.sum += mtmp.sum + bLength * mtmp.N;
      m.squaredSum += mtmp.squaredSum + 2 * bLength * mtmp.sum + mtmp.N * bLength * bLength;
    }
  }

  return m;
}
