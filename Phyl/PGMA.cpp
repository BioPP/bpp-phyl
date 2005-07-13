//
// File: PGMA.cpp
// Created by: Julien Dutheil
// Created on: Mon jul 11 11:41 2005
//

/*
Copyright ou � ou Copr. CNRS, (16 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant � fournir des classes
pour l'analyse de donn�es phylog�n�tiques.

Ce logiciel est r�gi par la licence CeCILL soumise au droit fran�ais et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffus�e par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilit� au code source et des droits de copie,
de modification et de redistribution accord�s par cette licence, il n'est
offert aux utilisateurs qu'une garantie limit�e.  Pour les m�mes raisons,
seule une responsabilit� restreinte p�se sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les conc�dants successifs.

A cet �gard  l'attention de l'utilisateur est attir�e sur les risques
associ�s au chargement,  � l'utilisation,  � la modification et/ou au
d�veloppement et � la reproduction du logiciel par l'utilisateur �tant 
donn� sa sp�cificit� de logiciel libre, qui peut le rendre complexe � 
manipuler et qui le r�serve donc � des d�veloppeurs et des professionnels
avertis poss�dant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invit�s � charger  et  tester  l'ad�quation  du
logiciel � leurs besoins dans des conditions permettant d'assurer la
s�curit� de leurs syst�mes et ou de leurs donn�es et, plus g�n�ralement, 
� l'utiliser et l'exploiter dans les m�mes conditions de s�curit�. 

Le fait que vous puissiez acc�der � cet en-t�te signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accept� les
termes.
*/

/*
Copyright or � or Copr. CNRS, (November 16, 2004)

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

#include "PGMA.h"
#include "NodeTemplate.h"
#include "Tree.h"
#include "TreeTemplate.h"
#include <cmath>
#include <iostream>
using namespace std;

#if defined(VIRTUAL_COV)
		TreeTemplate<Node> * 
#else
		Tree *
#endif
PGMA::getTree() const
{
	Node * root = TreeTools::cloneSubtree<Node>(* dynamic_cast<TreeTemplate<NodeTemplate<PGMAInfos> > *>(_tree) -> getRootNode());
	return new TreeTemplate<Node>(* root);
}
	
vector<unsigned int> PGMA::getBestPair()
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

