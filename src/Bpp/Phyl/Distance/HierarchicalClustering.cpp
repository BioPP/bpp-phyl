//
// File: HierarchicalClustering.cpp
// From file Cluster.cpp in CoMap package.
// Created by: Julien Dutheil
// Created on: Tue Aug 30 17:19 2005
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

   This software is a computer program whose purpose is to map substitutions
   on a tree and to detect co-evolving positions in a dataset.

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

#include "HierarchicalClustering.h"
#include "../Tree/NodeTemplate.h"

using namespace bpp;
using namespace std;

const string HierarchicalClustering::COMPLETE = "Complete";
const string HierarchicalClustering::SINGLE   = "Single";
const string HierarchicalClustering::AVERAGE  = "Average";
const string HierarchicalClustering::MEDIAN   = "Median";
const string HierarchicalClustering::WARD     = "Ward";
const string HierarchicalClustering::CENTROID = "Centroid";

TreeTemplate<Node>* HierarchicalClustering::getTree() const
{
  Node* root = TreeTemplateTools::cloneSubtree<Node>(*dynamic_cast<TreeTemplate<NodeTemplate<ClusterInfos> >*>(tree_)->getRootNode());
  return new TreeTemplate<Node>(root);
}

vector<size_t> HierarchicalClustering::getBestPair() throw (Exception)
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
  return w1 * d1 + w2 * d2 + w3 * d3 + w4* std::abs(d1 - d2);
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
  tree_ = new TreeTemplate<NodeTemplate<ClusterInfos> >(root);
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

