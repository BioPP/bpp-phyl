//
// File: NeighborJoining.cpp
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Thu jun 23 10:39 2005
//

/*
Copyright ou © ou Copr. CNRS, (16 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant à fournir des classes
pour l'analyse de données phylogénétiques.

Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.
*/

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

Julien.Dutheil@univ-montp2.fr

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

#include "NeighborJoining.h"
#include "Tree.h"
#include <cmath>
#include <iostream>
using namespace std;

vector<unsigned int> NeighborJoining::getBestPair()
{
	for(map<unsigned int, Node *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++) {
		unsigned int id = i -> first;
		_sumDist[id] = 0;
		for(map<unsigned int, Node *>::iterator j = _currentNodes.begin(); j != _currentNodes.end(); j++) {
			unsigned int jd = j -> first;
			_sumDist[id] += _matrix(id, jd);
		}
	}

	vector<unsigned int> bestPair(2);
	double critMax = std::log(0.);
	for(map<unsigned int, Node *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++) {
		unsigned int id = i -> first;
		map<unsigned int, Node *>::iterator j = i;
		j++;
		for(; j != _currentNodes.end(); j++) {
			unsigned int jd = j -> first;
			double crit = _sumDist[id] + _sumDist[jd] - (_currentNodes.size() - 2) * _matrix(id, jd);
			//cout << "\t" << id << "\t" << jd << "\t" << crit << endl;
			if(crit > critMax) {
				critMax = crit;
				bestPair[0] = id;
				bestPair[1] = jd;
			}
		}
	}

	return bestPair;	
}

vector<double> NeighborJoining::computeBranchLengthsForPair(const vector<unsigned int> & pair)
{
	double ratio = (_sumDist[pair[0]] - _sumDist[pair[1]]) / (_currentNodes.size() - 2);
	vector<double> d(2);
	d[0] = std::max(.5 * (_matrix(pair[0], pair[1]) + ratio), 0.); 
	d[1] = std::max(.5 * (_matrix(pair[0], pair[1]) - ratio), 0.); 
	return d;
}

double NeighborJoining::computeDistancesFromPair(const vector<unsigned int> & pair, const vector<double> & branchLengths, unsigned int pos)
{
	return std::max(.5 * (_matrix(pair[0], pos) - branchLengths[0] + _matrix(pair[1], pos) - branchLengths[1]), 0.); 
}

void NeighborJoining::finalStep(int idRoot)
{
	Node * root = new Node(idRoot);
	map<unsigned int, Node* >::iterator it = _currentNodes.begin();
	unsigned int i1 = it -> first;
	Node * n1       = it -> second;
	it++;
	unsigned int i2 = it -> first;
	Node * n2       = it -> second;
	if(_currentNodes.size() == 2) { //Rooted
		double d = _matrix(i1, i2) / 2;
		root -> addSon(*n1);
		root -> addSon(*n2);
		n1 -> setDistanceToFather(d);
		n2 -> setDistanceToFather(d);
	} else { //Unrooted
		it++;
		unsigned int i3 = it -> first;
		Node * n3       = it -> second;
		double d1 = std::max(_matrix(i1, i2) + _matrix(i1, i3) - _matrix(i2, i3), 0.);
		double d2 = std::max(_matrix(i2, i1) + _matrix(i2, i3) - _matrix(i1, i3), 0.);
		double d3 = std::max(_matrix(i3, i1) + _matrix(i3, i2) - _matrix(i1, i2), 0.);
		root -> addSon(*n1);
		root -> addSon(*n2);
		root -> addSon(*n3);
		n1 -> setDistanceToFather(d1);
		n2 -> setDistanceToFather(d2);
		n3 -> setDistanceToFather(d3);
	}
	_tree = new TreeTemplate<Node>(*root);
}

