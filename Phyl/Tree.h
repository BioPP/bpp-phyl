//
// File: Tree.h
// Created by: Julien Dutheil
// Created on: Thu Mar 13 12:03:18 2003
//

/*
Copyright ou © ou Copr. Julien Dutheil, (16 Novembre 2004) 

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
Copyright or © or Copr. Julien Dutheil, (November 16, 2004)

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

#ifndef _TREE_H_
#define _TREE_H_

#include "TreeExceptions.h"
#include "TreeTools.h"


// From the STL:
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

using namespace std;

/**
 * @brief The phylogenetic tree class.
 * 
 * This class is part of the object implementation of phylogenetic trees. Tree are made
 * made of nodes, instances of the class Node.
 * 
 * Trees are oriented (rooted), i.e. each node has one <i>father node</i> and possibly
 * many <i>son nodes</i>. Leaves are nodes without descendant and root is defined has the without
 * father. Inner nodes will generally contain two descendants (the tree is then called
 * <i>bifurcating</i>), but mutlifurcating trees are also allowed with this kind of description.
 * In the rooted case, each inner node also defines a <i>subtree</i>.
 * This allows to work recursively on trees, which is very convenient in most cases.
 * To deal with non-rooted trees, we place an artificial root at a particular node:
 * hence the root node appears to be trifurcated. This is the way unrooted trees are
 * described in the parenthetic description, the so called Newick format.
 * 
 * To clone a tree from from another tree with a different template,
 * consider using the TreeTools::cloneSutree<N>() method:
 * <code>
 * Tree * t = new Tree<Node>(...)
 * NodeTemplate<int> * newRoot = TreeTools::cloneSubtree< NodeTemplate<int> >(* (t -> getRootNode()))
 * Tree< NodeTemplate<int> > * tt = new Tree< NodeTemplate<int> >(* newRoot);
 * </code>
 * 
 * @see Node
 * @see NodeTemplate
 * @see TreeTools
 */
template<class N=Node>
class Tree {

	/**
	 * Fields:
	 */
	protected:
		N * _root;
		string _name;

	public: // Constructors and destructor:
		
		Tree() { _root = NULL; }

		Tree(const Tree<N> & t)
		{
			//Perform a hard copy of the nodes:
			_root = TreeTools::cloneSubtree<N>(* t.getRootNode());
		}

		Tree(N & root) { _root = &root; }

		Tree<N> & operator=(const Tree<N> & t)
		{
			//Perform a hard copy of the nodes:
			_root = TreeTools::cloneSubtree<N>(* t.getRootNode());
    	return *this;
		}

		virtual ~Tree() { destroyNode(* _root); delete _root; }

			
	/**
	 * Methods:
 	 */
	
	public:
		
		virtual N * getRootNode() { return _root; }

		virtual const N * getRootNode() const { return _root; }

		virtual void setRootNode(N & root) { _root = & root; }
	
		virtual string getName() const { return _name; }
	
		virtual void setName(const string & name) { _name = name; }

		virtual unsigned int getNumberOfLeaves() const { return TreeTools::getNumberOfLeaves(* _root); }

		virtual vector<const N *> getLeaves() const { return TreeTools::getLeaves(* const_cast<const N *>(_root)); }

		virtual vector<      N *> getLeaves()       { return TreeTools::getLeaves(* _root); }


		virtual vector<const N *> getNodes() const { return TreeTools::getNodes (* const_cast<const N *>(_root)); }

		virtual vector<      N *> getNodes()       { return TreeTools::getNodes (* _root); }

		virtual vector<double> getBranchLengths() const { return TreeTools::getBranchLengths(* _root); }

		virtual vector<string> getLeavesNames() const { return TreeTools::getLeavesNames(* const_cast<const N *>( _root)); }

		//void setNewOutgroup(Tree::N & node);

		// Works on root:
		
		/**
		 * @brief <p></p>
		 * @param p_iNewRoot The node to be considered as the new root, i.e.
		 * defining the subtree to be considered as the new outgroup.
		 */
		void rootAt(N & p_iNewRoot)
	  {
			if (* _root == p_iNewRoot) return;
			vector<Node *> pathMatrix = TreeTools::getPathBetweenAnyTwoNodes(* _root, p_iNewRoot);
			//pathMatrix size is always bigger than 2.

			for (unsigned int i = 0; i < pathMatrix.size() - 1 ; i++) {
				pathMatrix[i] -> _father = pathMatrix[i + 1];
				//pathMatrix[i] -> setFather(*pathMatrix[i + 1]);
				pathMatrix[i] -> setDistanceToFather(pathMatrix[i + 1] -> getDistanceToFather());
				typename vector<Node *>::iterator vec_iter;
				vec_iter = remove(pathMatrix[i] -> _sons.begin(), pathMatrix[i] -> _sons.end(), pathMatrix[i + 1]);
				pathMatrix[i] -> _sons.erase(vec_iter, pathMatrix[i] -> _sons.end()); // pg 1170, primer.
	
				pathMatrix[i+1] -> _sons.push_back(pathMatrix[i + 1] -> getFather());
				pathMatrix[i+1] -> _father = NULL;
				//pathMatrix[i+1] -> deleteFather();
			}
			_root = & p_iNewRoot;
		}

		
		/**
		 * @brief Tell if the tree is rooted.
		 * 
		 * @return True if the tree is rooted.
		 */
		virtual bool isRooted() const { return _root -> getNumberOfSons() == 2; }
		
		/**
		 * @brief Unroot a rooted tree.
		 *
		 * @return True if the tree has been unrooted.
		 * @throw UnrootedTreeException If the tree is already rooted.
		 */
		virtual bool unroot() throw (UnrootedTreeException<N>)
		{
			if(!isRooted()) throw UnrootedTreeException<N>("Tree::unroot", this);

    		if(_root -> getNumberOfSons() == 2) {
				N* son1 = _root -> getSon(0);
				N* son2 = _root -> getSon(1);
				if(son1 -> isLeaf() && son2 -> isLeaf()) return false; // We can't unroot a single branch!
					// We manage to have a subtree in position 0:
					if(son1 -> isLeaf()) {
						_root -> swap(0, 1);
						son1 = _root -> getSon(0);
						son2 = _root -> getSon(1);
					}

					// Take care of branch lengths:
					if(son1 -> hasDistanceToFather()) {
						if(son2 -> hasDistanceToFather()) {
						// Both nodes have lengths, we sum them:
						son2 -> setDistanceToFather(son1 -> getDistanceToFather() + son2 -> getDistanceToFather());
					} else {
						// Only node 1 has length, we set it to node 2:
						son2 -> setDistanceToFather(son1 -> getDistanceToFather());
					}
					son1 -> deleteDistanceToFather();
				} // Else node 2 may or may not have a branch length, we do not care!

				// Remove the root:
				_root -> removeSons();
				son1 -> addSon(*son2);
				setRootNode(*son1);
				return true;
			} else return false; // Tree is already rooted.
		}


		
		/**
		 * @brief Number nodes.
		 */
		virtual void resetNodesId()
		{
			vector<N *> nodes = getNodes();
			for(unsigned int i = 0; i < nodes.size(); i++) nodes[i] -> setId(i);
		}
		
		// Works on (multi)furcations:
		
		/**
		 * @brief Tell if the tree is multifurcating.
		 * 
		 * @return True if the tree is multifurcating.
		 */
		virtual bool isMultifurcating() const
		{
			bool b = false;
			for(unsigned int i = 0; i < _root -> getNumberOfSons(); i++) {
				b = b || TreeTools::isMultifurcating(* _root -> getSon(i));
			}
			return b;
		}
		
		/**
		 * @brief Get all the branch lengths of a tree.
		 *
		 * @return A vector with all branch lengths.
		 * @throw NodeException If a branch length is lacking.
		 */
		virtual Vdouble getBranchLengths() throw (NodeException)
		{
			Vdouble brLen(1);
			for(unsigned int i = 0; i < _root -> getNumberOfSons(); i++) {
				Vdouble sonBrLen = TreeTools::getBranchLengths(* _root -> getSon(i));
				for(unsigned int j = 0; j < sonBrLen.size(); j++) brLen.push_back(sonBrLen[j]);
			}
			return brLen;
		}

		/**
		 * @brief Get the total length (sum of all branch lengths) of a tree.
		 *
		 * @return The total length of the subtree.
		 * @throw NodeException If a branch length is lacking.
		 */
		virtual double getTotalLength() throw (NodeException)
		{
			return TreeTools::getTotalLength(*_root);
		}

		/**
		 * @brief Set all the branch lengths of a tree.
		 *
		 * @param brLen The branch length to apply.
		 */
		virtual void setBranchLengths(double brLen)
		{
			TreeTools::setBranchLengths(*_root, brLen);
		}
		
		/**
		 * @brief Give a length to branches that don't have one in a tree.
		 *
		 * @param brLen The branch length to apply.
		 */
		virtual void setVoidBranchLengths(double brLen)
		{
			TreeTools::setVoidBranchLengths(*_root, brLen);
		}
	
		/**
		 * @brief Scale a given tree.
		 *
		 * Multiply all branch lengths by a given factor.
		 *
		 * @param factor The factor to multiply all branch lengths with.
		 * @throw NodeException If a branch length is lacking.
		 */
		virtual void scaleTree(double factor) throw (NodeException)
		{
			TreeTools::scaleTree(* _root, factor);
		}

	protected:
		
		virtual void destroyNode(const N & node)
		{
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
				destroyNode(* node[i]);
				delete node[i];
			}
		}
		
};

#endif	//_TREE_H_

