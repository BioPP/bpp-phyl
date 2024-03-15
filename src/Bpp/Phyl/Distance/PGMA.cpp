// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../Tree/NodeTemplate.h"
#include "../Tree/Tree.h"
#include "../Tree/TreeTemplate.h"
#include "../Tree/TreeTemplateTools.h"
#include "PGMA.h"

using namespace bpp;

// From the STL:
#include <cmath>
#include <iostream>

using namespace std;

vector<size_t> PGMA::getBestPair()
{
  vector<size_t> bestPair(2);
  double distMin = -std::log(0.);
  for (map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
  {
    size_t id = i->first;
    map<size_t, Node*>::iterator j = i;
    j++;
    for ( ; j != currentNodes_.end(); j++)
    {
      size_t jd = j->first;
      double dist = matrix_(id, jd);
      if (dist < distMin)
      {
        distMin = dist;
        bestPair[0] = id;
        bestPair[1] = jd;
      }
    }
  }

  if (distMin == -std::log(0.))
  {
    throw Exception("Unexpected error: no minimum found in the distance matrix.");
  }

  return bestPair;
}

vector<double> PGMA::computeBranchLengthsForPair(const vector<size_t>& pair)
{
  vector<double> d(2);
  double dist = matrix_(pair[0], pair[1]) / 2.;
  d[0] = dist - dynamic_cast<NodeTemplate<PGMAInfos>*>(currentNodes_[pair[0]])->getInfos().time;
  d[1] = dist - dynamic_cast<NodeTemplate<PGMAInfos>*>(currentNodes_[pair[1]])->getInfos().time;
  return d;
}

double PGMA::computeDistancesFromPair(const vector<size_t>& pair, const vector<double>& branchLengths, size_t pos)
{
  double w1, w2;
  if (weighted_)
  {
    w1 = 1;
    w2 = 1;
  }
  else
  {
    w1 = static_cast<double>(dynamic_cast<NodeTemplate<PGMAInfos>*>(currentNodes_[pair[0]])->getInfos().numberOfLeaves);
    w2 = static_cast<double>(dynamic_cast<NodeTemplate<PGMAInfos>*>(currentNodes_[pair[1]])->getInfos().numberOfLeaves);
  }
  return (w1 * matrix_(pair[0], pos) + w2 * matrix_(pair[1], pos)) / (w1 + w2);
}

void PGMA::finalStep(int idRoot)
{
  NodeTemplate<PGMAInfos>* root = new NodeTemplate<PGMAInfos>(idRoot);
  map<size_t, Node*>::iterator it = currentNodes_.begin();
  size_t i1 = it->first;
  Node* n1        = it->second;
  it++;
  size_t i2 = it->first;
  Node* n2        = it->second;
  double d = matrix_(i1, i2) / 2;
  root->addSon(n1);
  root->addSon(n2);
  n1->setDistanceToFather(d - dynamic_cast<NodeTemplate<PGMAInfos>*>(n1)->getInfos().time);
  n2->setDistanceToFather(d - dynamic_cast<NodeTemplate<PGMAInfos>*>(n2)->getInfos().time);
  tree_.reset(new TreeTemplate<NodeTemplate<PGMAInfos>>(root));
}

Node* PGMA::getLeafNode(int id, const string& name)
{
  PGMAInfos infos;
  infos.numberOfLeaves = 1;
  infos.time = 0.;
  NodeTemplate<PGMAInfos>* leaf = new NodeTemplate<PGMAInfos>(id, name);
  leaf->setInfos(infos);
  return leaf;
}

Node* PGMA::getParentNode(int id, Node* son1, Node* son2)
{
  PGMAInfos infos;
  infos.numberOfLeaves =
    dynamic_cast<NodeTemplate<PGMAInfos>*>(son1)->getInfos().numberOfLeaves
    + dynamic_cast<NodeTemplate<PGMAInfos>*>(son2)->getInfos().numberOfLeaves;
  infos.time = dynamic_cast<NodeTemplate<PGMAInfos>*>(son1)->getInfos().time + son1->getDistanceToFather();
  Node* parent = new NodeTemplate<PGMAInfos>(id);
  dynamic_cast<NodeTemplate<PGMAInfos>*>(parent)->setInfos(infos);
  parent->addSon(son1);
  parent->addSon(son2);
  return parent;
}
