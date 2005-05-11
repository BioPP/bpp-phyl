//
// File: Tree.h
// Created by: Julien Dutheil <julien.dutheil@ens-lyon.fr>
// Created on: Thu Mar 13 12:03:18 2003
//

/*
Copyright ou � ou Copr. Julien Dutheil, (16 Novembre 2004) 

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
Copyright or � or Copr. Julien Dutheil, (November 16, 2004)

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

#include "Tree.h"
#include "TreeTools.h"

//From Utils:
#include <Utils/Exceptions.h>

//from the STL:
#include <algorithm>
#include <iostream>
using namespace std;

/******************************************************************************
 *                               The node class                               *
 ******************************************************************************/

/** Copy constructor: *********************************************************/
	
Node::Node(const Node & node)
{
	_id               = node._id;
	_name             = node.hasName() ? new string(* node._name) : NULL;
	_father           = node._father;
	_distanceToFather = node.hasDistanceToFather() ? new double(* node._distanceToFather) : NULL;
	_sons             = node._sons;
	for(map<string, Clonable *>::iterator i = node._properties.begin(); i != node._properties.end(); i++)
		_properties[i -> first] = i -> second -> clone();
}

/** Assignation operator: *****************************************************/

Node & Node::operator=(const Node & node)
{
	_id               = node._id;
	_name             = node.hasName() ? new string(* node._name) : NULL;
	_father           = node._father;
	_distanceToFather = node.hasDistanceToFather() ? new double(* node._distanceToFather) : NULL;
	_sons             = node._sons;
	for(map<string, Clonable *>::iterator i = node._properties.begin(); i != node._properties.end(); i++)
		_properties[i -> first] = i -> second -> clone();
	return * this;
}
			
/** Sons: *********************************************************************/
			
void Node::swap(unsigned int branch1, unsigned int branch2)
{
    Node* node1 = getSon(branch1);
    Node* node2 = getSon(branch2);
    removeSon(*node1);
    removeSon(*node2);
    addSon(branch1, *node2);
    addSon(branch2, *node1);
}

vector<const Node *> Node::getNeighbors() const
{
	vector<const Node *> neighbors(_sons.size() + 1);
	neighbors[0] = _father;
	for(unsigned int i = 0; i < _sons.size(); i++) neighbors[i+1] = _sons[i];
	return neighbors;
}
		
vector<Node *> Node::getNeighbors()
{
	vector<Node *> neighbors(_sons.size() + 1);
	neighbors[0] = _father;
	for(unsigned int i = 0; i < _sons.size(); i++) neighbors[i+1] = _sons[i];
	return neighbors;
}
	
/******************************************************************************/






/*******************************************************************************
 *                                The tree class                               *
 *******************************************************************************/

/*
void Tree::setNewOutgroup(Tree::Node & node) {
	//cout << "New outgroup " << node.getId() << endl;
	if(node == * _root) return; //can't reroot with root!
	Node * father = node.getFather();
	if(father == _root) return; //nothing to do, already rooted here!
	Node * grandFather = father -> getFather();
	if(grandFather != _root) setNewOutgroup(* father);//must reroot with father first!
	//cout << "dealing with node " << node.getId() << endl;
	
	father -> eraseSon(node);
	_root  -> eraseSon(* father);
	Node * uncle = (* _root)[0];
	_root  -> eraseSon(0);//remove uncle.
	uncle  -> setDistanceToFather(uncle -> getDistanceToFather() + father -> getDistanceToFather());
	father -> setDistanceToFather(node.getDistanceToFather() / 2);//Place new root at mid length between node and father.
	node   .  setDistanceToFather(node.getDistanceToFather() / 2);
	father -> addSon(* uncle);
	_root  -> addSon(node);
	_root  -> addSon(* father);//We add it again ;-)
}*/
