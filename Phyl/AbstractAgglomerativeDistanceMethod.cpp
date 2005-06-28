//
// File: AbstractAgglomerativeDistanceMethod.cpp
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Wed jun 22 10:00 2005
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

#include "AbstractAgglomerativeDistanceMethod.h"
// From the STL:
#include <iostream>
using namespace std;

AbstractAgglomerativeDistanceMethod::~AbstractAgglomerativeDistanceMethod()
{
	delete _tree;
}

void AbstractAgglomerativeDistanceMethod::setDistanceMatrix(const DistanceMatrix & matrix)
{
	_matrix = matrix;
	if(_tree != NULL) delete _tree;
}
	
Tree<Node> * AbstractAgglomerativeDistanceMethod::getTree() const
{
	Node * root = TreeTools::cloneSubtree<Node>(* _tree -> getRootNode());
	return new Tree<Node>(* root);
}
			
void AbstractAgglomerativeDistanceMethod::computeTree(bool rooted)
{
	// Initialization:
	for(unsigned int i = 0; i < _matrix.size(); i++) {
		N * leaf = new N(i, _matrix.getName(i));
		leaf -> setInfos(1);
		_currentNodes[i] = leaf;
	}
	unsigned int idNextNode = _matrix.size();
	vector<double> newDist(_matrix.size());
	
	// Build tree:
	while(_currentNodes.size() > (rooted ? 2 : 3)) {
		vector<unsigned int> bestPair = getBestPair();
		vector<double> distances = computeBranchLengthsForPair(bestPair);
//cout << bestPair[0] << "\t" << bestPair[1] << endl;
		N * best1 = _currentNodes[bestPair[0]];
		N * best2 = _currentNodes[bestPair[1]];
//cout << "id\t" << best1 -> getId() << "\t" << best2 -> getId() << endl;
		N * parent = new N(idNextNode++);
		parent -> addSon(* best1);
		parent -> addSon(* best2);
		best1  -> setDistanceToFather(distances[0]);
		best2  -> setDistanceToFather(distances[1]);
		parent -> setInfos(best1 -> getInfos() + best2 -> getInfos());
		for(map<unsigned int, N *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++) {
			unsigned int id = i -> first;
			if(id != bestPair[0] && id != bestPair[1]) {
				newDist[id] = computeDistancesFromPair(bestPair, distances, id);
			} else {
				newDist[id] = 0;
			}
		}
		// Actualize _currentNodes:
		_currentNodes[bestPair[0]] = parent;
		_currentNodes.erase(bestPair[1]);
		for(map<unsigned int, N *>::iterator i = _currentNodes.begin(); i != _currentNodes.end(); i++) {
			unsigned int id = i -> first;
			_matrix(bestPair[0], id) = _matrix(id, bestPair[0]) = newDist[id];
		}
		
	}
	finalStep(idNextNode);
}


