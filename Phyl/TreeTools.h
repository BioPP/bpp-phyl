//
// File: TreeTools.h
// Created by:  Julien Dutheil
// Created on: Wed Aug  6 13:45:28 2003
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

#ifndef _TREETOOLS_H_
#define _TREETOOLS_H_

#include "TreeExceptions.h"
#include "Node.h"
#include "DistanceMatrix.h"
template<class N> class TreeTemplate;

// From Utils:
#include <Utils/Exceptions.h>

// From NumCalc:
#include <NumCalc/VectorTools.h>

/**
 * @brief A set of tools to deal with trees.
 */
class TreeTools
{
	public:
		TreeTools();
		 ~TreeTools();
	
	public:
		 
		/**
		 * @brief Retrieve all leaves from a subtree.
		 *
		 * @param node The node that defines the subtree.
		 * @return A vector of pointers toward each leaf in the subtree.
		 */
		template<class N>
		static vector<N *> getLeaves(N & node)
		{
			vector<N *> leaves;
			getLeaves<N>(node, leaves);
			return leaves;
		}

		template<class N>
		static void getLeaves(N & node, vector<N *> & leaves)
		{
			if(node.isLeaf()) {
				leaves.push_back(& node);
			}
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
				getLeaves<N>(* node.getSon(i), leaves);
			}
		}

		/**
		 * @brief Retrieve all leaves ids from a subtree.
		 *
		 * @param node The node that defines the subtree.
		 * @return A vector of ids.
		 */
		static vector<int> getLeavesId(const Node & node)
		{
			vector<int> ids;
			getLeavesId(node, ids);
			return ids;
		}

		static void getLeavesId(const Node & node, vector<int> & ids)
		{
			if(node.isLeaf()) {
				ids.push_back(node.getId());
			}
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
				getLeavesId(* node.getSon(i), ids);
			}
		}

		/**
		 * @brief Get the id of a leaf given its name in a subtree.
		 *
		 * @param node The node defining the subtree to search.
		 * @param name The name of the node.
		 * @return The id of the node.
		 * @throw NodeNotFoundException If the node is not found.
		 */
		static int getLeafId(const Node & node, const string & name) throw (NodeNotFoundException)
		{
			int * id = NULL;
			searchLeaf(node, name, id);
			if(id == NULL) throw NodeNotFoundException("TreeTools::getLeafId().", name);
			else {
				int i = *id;
				delete id;
				return i;
			}
		}

		static void searchLeaf(const Node & node, const string & name, int * & id) throw (NodeNotFoundException)
		{
			if(node.isLeaf()) {
				if(node.getName() == name) {
					id = new int(node.getId());
					return;
				}
			}
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
				searchLeaf(* node.getSon(i), name, id);
			}
		}

		/**
		 * @brief Retrieve all son nodes from a subtree.
		 *
		 * @param node The node that defines the subtree.
		 * @return A vector of pointers toward each son node in the subtree.
		 */
		template<class N>
		static vector<N *> getNodes(N & node)
		{
			vector<N *> nodes;
			getNodes<N>(node, nodes);
			return nodes;
		}

		template<class N>
		static void getNodes(N & node, vector<N *> & nodes)
		{
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
				getNodes<N>(* node.getSon(i), nodes);
			}
			nodes.push_back(& node);
		}

		/**
		 * @brief Retrieve all son nodes ids from a subtree.
		 *
		 * @param node The node that defines the subtree.
		 * @return A vector of ids.
		 */
		static vector<int> getNodesId(const Node & node)
		{
			vector<int> ids;
			getNodesId(node, ids);
			return ids;
		}

		static void getNodesId(const Node & node, vector<int> & ids)
		{
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
				getNodesId(* node.getSon(i), ids);
			}
			ids.push_back(node.getId());
		}

		template<class N>
		static vector<N *> searchNodeWithId(N & node, int id)
		{
			vector<N *> nodes;
			searchNodeWithId<N>(node, id, nodes);
			return nodes;		
		}

		template<class N>
		static void searchNodeWithId(N & node, int id, vector<N *> & nodes)
		{
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
				searchNodeWithId<N>(* node.getSon(i), id, nodes);
			}
			if(node.getId() == id) nodes.push_back(& node);
		}

		/**
		 * @brief Tell if a particular node is the root of a tree
		 * i.e. if it has a father node.
		 *
		 * @param node The node to check.
		 * @return True if node has a father.
		 */
		static bool isRoot(const Node & node);

		/**
		 * @brief Get the number of leaves of a subtree defined by a particular node.
		 *
		 * @param node The node defining the subtree to check.
		 * @return The number of leaves.
		 */
		static unsigned int getNumberOfLeaves(const Node & node);

    /**
     * @brief Get the leaves names of a subtree defined by a particular node.
     *
		 * @param node The node defining the subtree to check.
		 * @return The lst of all leaves names.
		 */
		static vector<string> getLeavesNames(const Node & node);

		/**
		 * @brief Get the depth of the subtree defined by node 'node', i.e. the maximum
		 * number of sons 'generations'.
		 *
		 * ex:
		 * <code>
		 *    +----------A
		 *    |
		 * ---+ N1     +-------B
		 *    |        |
		 *    +--------+ N2
		 *             |
		 *             +------C
		 * </code>
		 * Depth of node 'N1' id 2, depth of node 'N2' is 1, depth of leaves is 0.
		 */
		static unsigned int getDepth(const Node & node);

		/**
		 * @name Conversion tools.
		 *
		 * Convert from Newick standard tree description.
		 * The description is for a node, and hence is to be rounded with
		 * parenthesis. ex: (A:0.001, (B:0.001, C:0.02):0.005):0.0005
		 *
		 * @{
		 */

	private: 
		struct Element
		{
			string content;
			double * length;
			double * bootstrap;
		};

		static Element getElement(string elt) throw (IOException);

	public:
		/**
		 * @brief Parse a string in the parenthesis format and convert it to
		 * a subtree.
		 *
		 * @param description the string to parse;
		 * @return A pointer toward a dynamically created subtree.
		 */
		static Node * parenthesisToNode(const string & description);
	
		/**
		 * @brief Parse a string in the parenthesis format and convert it to
		 * a tree.
		 *
		 * @param description the string to parse;
		 * @return A pointer toward a dynamically created tree.
		 */
		static TreeTemplate<Node> * parenthesisToTree(const string & description);
		
		/**
		 * @brief Get the parenthesis description of a subtree.
		 *
		 * @param node The node defining the subtree.
		 * @return A string in the parenthesis format.
		 */
		static string nodeToParenthesis(const Node & node);
		static string nodeToParenthesis(const Tree & tree, int nodeId);

		/**
		 * @brief Get the parenthesis description of a tree.
		 *
		 * @param tree The tree to convert.
		 * @return A string in the parenthesis format.
		 */
		static string treeToParenthesis(const TreeTemplate<Node> & tree);
		static string treeToParenthesis(const Tree & tree);
		
		/** @} */
		
		/**
		 * @brief Tell is a subtree is multifurcating.
		 *
		 * @param node The root node of the subtree.
		 * @return True is the subtree contains at least one multifurcating
		 * node (including the root node).
		 */
		static bool isMultifurcating(const Node & node);
		
		/**
		 * @name Act on branch lengths.
		 *
		 * @{
		 */
		
		/**
		 * @brief Get all the branch lengths of a subtree.
		 *
		 * @param node The root node of the subtree.
		 * @return A vector with all branch lengths.
		 * @throw NodeException If a branch length is lacking.
		 */
		static Vdouble getBranchLengths(const Node & node) throw (NodeException);
		
		/**
		 * @brief Get the total length (sum of all branch lengths) of a subtree.
		 *
		 * @param node The root node of the subtree.
		 * @return The total length of the subtree.
		 * @throw NodeException If a branch length is lacking.
		 */
		static double getTotalLength(const Node & node) throw (NodeException);

		/**
		 * @brief Set all the branch lengths of a subtree.
		 *
		 * @param node  The root node of the subtree.
		 * @param brLen The branch length to apply.
		 */
		static void setBranchLengths(Node & node, double brLen);
		
		/**
		 * @brief Give a length to branches that don't have one in a subtree.
		 *
		 * @param node  The root node of the subtree.
		 * @param brLen The branch length to apply.
		 */
		static void setVoidBranchLengths(Node & node, double brLen);
				
		/**
		 * @brief Scale a given tree.
		 *
		 * Multiply all branch lengths by a given factor.
		 *
		 * @param node   The root node of the subtree to scale.
		 * @param factor The factor to multiply all branch lengths with.
		 * @throw NodeException If a branch length is lacking.
		 */
		static void scaleTree(Node & node, double factor) throw (NodeException);
				
		/** @} */

		//The following methods are adapted from the tree class of the SEMPHY library:
		static vector<Node *> getPathBetweenAnyTwoNodes(Node & node1, Node & node2, bool includeAncestor = true);
		
		static vector<const Node *> getPathBetweenAnyTwoNodes(const Node & node1, const Node & node2, bool includeAncestor = true);
		
		static vector<int> getPathBetweenAnyTwoNodes(const Tree & tree, int nodeId1, int nodeId2, bool includeAncestor = true);
	
		template<class N>
		static N * cloneSubtree(const Node & node) 
		{
			//First we copy this node using default copy constuctor:
			N * clone = new N(node);
			//Remove all sons. This is possible since Tree is a friend of Node:
			//clone -> sons.resize(0);
			//Now we perform a hard copy:
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
				//clone -> addSon(* cloneSubtree<N>(* node[i]));
				clone -> setSon(i, * cloneSubtree<N>(* node[i]));
			}
			return clone;
		}
		
		/**
		 * @name Random trees
		 *
		 * @{
		 */

		/**
		 * @brief Draw a random tree from a list of taxa.
		 *
		 * @param leavesNames A list of taxa.
		 * @return A random tree with all corresponding taxa.
		 */
		static TreeTemplate<Node> * getRandomTree(vector<string> & leavesNames);

		/** @} */
		
		/**
		 * @name Some properties.
		 *
		 * @{
		 */
		 
		/**
		 * @brief Bootstrap tag.
		 */
		static string BOOTSTRAP;
		
		/** @} */

		static double getDistanceBetweenAnyTwoNodes(const Node & node1, const Node & node2);
		/**
		 * @brief Get the total distance between two nodes.
		 *
		 * Sum all branch lengths between two nodes.
		 *
		 * @param tree The tree to consider.
		 * @param nodeId1 First node id.
		 * @param nodeId2 Second node id.
		 * @return The sum of all branch lengths
		 */
		static double getDistanceBetweenAnyTwoNodes(const Tree & tree, int nodeId1, int nodeId2);
		
		/**
		 * @brief Compute a distance matrix from a tree.
		 *
		 * Compute all distances between each leaves and store them in a matrix.
		 * A new DistanceMatrix object is created, and a pointer toward it is returned.
		 * The destruction of this matrix is left up to the user.
		 *
		 * @see getDistanceBetweenAnyTwoNodes
		 *
		 * @param tree The tree to use.
		 * @return The distance matrix computed from tree.
		 */
		static DistanceMatrix * getDistanceMatrix(const Tree & tree); 	

};


#endif	//_TREETOOLS_H_

