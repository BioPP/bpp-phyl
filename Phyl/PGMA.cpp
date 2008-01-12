//
// File: PGMA.cpp
// Created by: Julien Dutheil
// Created on: Mon jul 11 11:41 2005
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "PGMA.h"
#include "NodeTemplate.h"
#include "Tree.h"
#include "TreeTemplate.h"
#include "TreeTemplateTools.h"

using namespace bpp;

// From the STL:
#include <cmath>
#include <iostream>

using namespace std;

TreeTemplate<Node> * PGMA::getTree() const
{
	Node * root = TreeTemplateTools::cloneSubtree<Node>(* dynamic_cast<TreeTemplate<NodeTemplate<PGMAInfos> > *>(_tree) -> getRootNode());
	return new TreeTemplate<Node>(* root);
}
	
vector<unsigned int> PGMA::getBestPair() throw (Exception)
{
	vector<unsigned int> bestPair(2);
	double distMin = -std::log(0.);
	for(map<unsigned int, Node *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++) {
		unsigned int id = i -> first;
		map<unsigned int, Node *>::iterator j = i;
		j++;
		for(; j != _currentNodes.end(); j++) {
			unsigned int jd = j -> first;
			double dist = _matrix(id, jd);
			if(dist < distMin) {
				distMin = dist;
				bestPair[0] = id;
				bestPair[1] = jd;
			}
		}
	}

	if(distMin == -std::log(0.)) {
		throw Exception("Unexpected error: no minimum found in the distance matrix.");
	}

	return bestPair;	
}

vector<double> PGMA::computeBranchLengthsForPair(const vector<unsigned int> & pair)
{
	vector<double> d(2);
	double dist = _matrix(pair[0], pair[1]) / 2.;
	d[0] = dist - dynamic_cast<NodeTemplate<PGMAInfos> *>(_currentNodes[pair[0]]) -> getInfos().time; 
	d[1] = dist - dynamic_cast<NodeTemplate<PGMAInfos> *>(_currentNodes[pair[1]]) -> getInfos().time; 
	return d;
}

double PGMA::computeDistancesFromPair(const vector<unsigned int> & pair, const vector<double> & branchLengths, unsigned int pos)
{
	double w1, w2;
	if(_weighted) {
		w1 = 1;
		w2 = 1;
	} else {
		w1 = dynamic_cast<NodeTemplate<PGMAInfos> *>(_currentNodes[pair[0]]) -> getInfos().numberOfLeaves;
		w2 = dynamic_cast<NodeTemplate<PGMAInfos> *>(_currentNodes[pair[1]]) -> getInfos().numberOfLeaves;
	}
	return (w1 * _matrix(pair[0], pos) + w2 * _matrix(pair[1], pos)) / (w1 + w2); 
}

void PGMA::finalStep(int idRoot)
{
	NodeTemplate<PGMAInfos> * root = new NodeTemplate<PGMAInfos>(idRoot);
	map<unsigned int, Node * >::iterator it = _currentNodes.begin();
	unsigned int i1 = it -> first;
	Node * n1       = it -> second;
	it++;
	unsigned int i2 = it -> first;
	Node * n2       = it -> second;
	double d = _matrix(i1, i2) / 2;
	root -> addSon(*n1);
	root -> addSon(*n2);
	n1 -> setDistanceToFather(d - dynamic_cast<NodeTemplate<PGMAInfos> *>(n1) -> getInfos().time); 
	n2 -> setDistanceToFather(d - dynamic_cast<NodeTemplate<PGMAInfos> *>(n2) -> getInfos().time); 
	_tree = new TreeTemplate<NodeTemplate<PGMAInfos> >(*root);
}

Node * PGMA::getLeafNode(int id, const string & name)
{
	PGMAInfos infos;
	infos.numberOfLeaves = 1;
	infos.time = 0.;
	NodeTemplate<PGMAInfos> * leaf = new NodeTemplate<PGMAInfos>(id, name);
	leaf -> setInfos(infos);
	return leaf;
}

Node * PGMA::getParentNode(int id, Node * son1, Node * son2)
{
	PGMAInfos infos;
	infos.numberOfLeaves = 
		dynamic_cast<NodeTemplate<PGMAInfos> *>(son1) -> getInfos().numberOfLeaves
	+ dynamic_cast<NodeTemplate<PGMAInfos> *>(son2) -> getInfos().numberOfLeaves;
	infos.time = dynamic_cast<NodeTemplate<PGMAInfos> *>(son1) -> getInfos().time + son1 -> getDistanceToFather();
	Node * parent = new NodeTemplate<PGMAInfos>(id);
	dynamic_cast<NodeTemplate<PGMAInfos> *>(parent) -> setInfos(infos);
	parent -> addSon(* son1);
	parent -> addSon(* son2);
	return parent;
}

