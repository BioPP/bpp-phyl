// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../Tree/NodeTemplate.h"
#include "HierarchicalClustering.h"

using namespace bpp;
using namespace std;

const string HierarchicalClustering::COMPLETE = "Complete";
const string HierarchicalClustering::SINGLE   = "Single";
const string HierarchicalClustering::AVERAGE  = "Average";
const string HierarchicalClustering::MEDIAN   = "Median";
const string HierarchicalClustering::WARD     = "Ward";
const string HierarchicalClustering::CENTROID = "Centroid";

vector<size_t> HierarchicalClustering::getBestPair()
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
    cout << "---------------------------------------------------------------------------------" << endl;
    for (map<size_t, Node*>::iterator i = currentNodes_.begin(); i != currentNodes_.end(); i++)
    {
      size_t id = i->first;
      map<size_t, Node*>::iterator j = i;
      j++;
      for ( ; j != currentNodes_.end(); j++)
      {
        size_t jd = j->first;
        double dist = matrix_(id, jd);
        cout << dist << "\t";
      }
      cout << endl;
    }
    cout << "---------------------------------------------------------------------------------" << endl;

    throw Exception("Unexpected error: no minimum found in the distance matrix.");
  }

  return bestPair;
}
vector<double> HierarchicalClustering::computeBranchLengthsForPair(const vector<size_t>& pair)
{
  vector<double> d(2);
  double dist = matrix_(pair[0], pair[1]) / 2.;
  d[0] = dist - dynamic_cast<NodeTemplate<ClusterInfos>*>(currentNodes_[pair[0]])->getInfos().length;
  d[1] = dist - dynamic_cast<NodeTemplate<ClusterInfos>*>(currentNodes_[pair[1]])->getInfos().length;
  return d;
}

double HierarchicalClustering::computeDistancesFromPair(const vector<size_t>& pair, const vector<double>& branchLengths, size_t pos)
{
  double w1, w2, w3, w4;
  if (method_ == "Single")
  {
    w1 = .5;
    w2 = .5;
    w3 = 0.;
    w4 = -.5;
  }
  else if (method_ == "Complete")
  {
    w1 = .5;
    w2 = .5;
    w3 = 0.;
    w4 = .5;
  }
  else if (method_ == "Median")
  {
    w1 = .5;
    w2 = .5;
    w3 = -0.25;
    w4 = 0.;
  }
  else if (method_ == "Average")
  {
    double n1 = static_cast<double>(dynamic_cast<NodeTemplate<ClusterInfos>*>(currentNodes_[pair[0]])->getInfos().numberOfLeaves);
    double n2 = static_cast<double>(dynamic_cast<NodeTemplate<ClusterInfos>*>(currentNodes_[pair[1]])->getInfos().numberOfLeaves);
    w1 = n1 / (n1 + n2);
    w2 = n2 / (n1 + n2);
    w3 = 0.;
    w4 = 0.;
  }
  else if (method_ == "Ward")
  {
    double n1 = static_cast<double>(dynamic_cast<NodeTemplate<ClusterInfos>*>(currentNodes_[pair[0]])->getInfos().numberOfLeaves);
    double n2 = static_cast<double>(dynamic_cast<NodeTemplate<ClusterInfos>*>(currentNodes_[pair[1]])->getInfos().numberOfLeaves);
    double n3 = static_cast<double>(dynamic_cast<NodeTemplate<ClusterInfos>*>(currentNodes_[pos])->getInfos().numberOfLeaves);
    w1 = (n1 + n3) / (n1 + n2 + n3);
    w2 = (n2 + n3) / (n1 + n2 + n3);
    w3 = -n3 / (n1 + n2 + n3);
    w4 = 0.;
  }
  else if (method_ == "Centroid")
  {
    double n1 = static_cast<double>(dynamic_cast<NodeTemplate<ClusterInfos>*>(currentNodes_[pair[0]])->getInfos().numberOfLeaves);
    double n2 = static_cast<double>(dynamic_cast<NodeTemplate<ClusterInfos>*>(currentNodes_[pair[1]])->getInfos().numberOfLeaves);
    w1 = n1 / (n1 + n2);
    w2 = n2 / (n1 + n2);
    w3 = -n1 * n2 / pow(n1 + n2, 2.);
    w4 = 0.;
  }
  else
    throw Exception("HierarchicalClustering::computeBranchLengthsForPair. unknown method '" + method_ + "'.");
  double d1 = matrix_(pair[0], pos);
  double d2 = matrix_(pair[1], pos);
  double d3 = matrix_(pair[0], pair[1]);
  return w1 * d1 + w2 * d2 + w3 * d3 + w4 * std::abs(d1 - d2);
}

void HierarchicalClustering::finalStep(int idRoot)
{
  NodeTemplate<ClusterInfos>* root = new NodeTemplate<ClusterInfos>(idRoot);
  map<size_t, Node*>::iterator it = currentNodes_.begin();
  size_t i1 = it->first;
  Node* n1        = it->second;
  it++;
  size_t i2 = it->first;
  Node* n2        = it->second;
  double d = matrix_(i1, i2) / 2;
  root->addSon(n1);
  root->addSon(n2);
  n1->setDistanceToFather(d - dynamic_cast<NodeTemplate<ClusterInfos>*>(n1)->getInfos().length);
  n2->setDistanceToFather(d - dynamic_cast<NodeTemplate<ClusterInfos>*>(n2)->getInfos().length);
  tree_.reset(new TreeTemplate<NodeTemplate<ClusterInfos>>(root));
}

Node* HierarchicalClustering::getLeafNode(int id, const string& name)
{
  ClusterInfos infos;
  infos.numberOfLeaves = 1;
  infos.length = 0.;
  NodeTemplate<ClusterInfos>* leaf = new NodeTemplate<ClusterInfos>(id, name);
  leaf->setInfos(infos);
  return leaf;
}

Node* HierarchicalClustering::getParentNode(int id, Node* son1, Node* son2)
{
  ClusterInfos infos;
  infos.numberOfLeaves =
      dynamic_cast<NodeTemplate<ClusterInfos>*>(son1)->getInfos().numberOfLeaves
      + dynamic_cast<NodeTemplate<ClusterInfos>*>(son2)->getInfos().numberOfLeaves;
  infos.length = dynamic_cast<NodeTemplate<ClusterInfos>*>(son1)->getInfos().length + son1->getDistanceToFather();
  Node* parent = new NodeTemplate<ClusterInfos>(id);
  dynamic_cast<NodeTemplate<ClusterInfos>*>(parent)->setInfos(infos);
  parent->addSon(son1);
  parent->addSon(son2);
  return parent;
}
