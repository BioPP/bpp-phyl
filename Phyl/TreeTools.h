//
// File: TreeTools.h
// Created by:  Julien Dutheil
// Created on: Wed Aug  6 13:45:28 2003
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

#ifndef _TREETOOLS_H_
#define _TREETOOLS_H_

#include "TreeExceptions.h"
#include "Node.h"

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
			if(node.isLeaf()) {
				leaves.push_back(& node);
			}
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
				vector<N *> sonLeaves = getLeaves(* node[i]);
				for(unsigned int j = 0; j < sonLeaves.size(); j++) {
					leaves.push_back(sonLeaves[j]);
				}
			}
			return leaves;
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
			for(unsigned int i = 0; i < node.getNumberOfSons(); i++) {
				vector<N *> sonNodes = getNodes(* node[i]);
				for(unsigned int j = 0; j < sonNodes.size(); j++) {
					nodes.push_back(sonNodes[j]);
				}
			}
			nodes.push_back(& node);
			return nodes;
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
		static Tree<Node> * parenthesisToTree(const string & description);
		
		/**
		 * @brief Get the parenthesis description of a subtree.
		 *
		 * @param node The node defining the subtree.
		 * @return A string in the parenthesis format.
		 */
		static string nodeToParenthesis(const Node & node);

		/**
		 * @brief Get the parenthesis description of a tree.
		 *
		 * @param tree The tree to convert.
		 * @return A string in the parenthesis format.
		 */
		static string treeToParenthesis(const Tree<Node> & tree);
		
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
		static vector<Node *> getPathBetweenAnyTwoNodes(Node & node1, Node & node2);
	
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
		static Tree<Node> * getRandomTree(vector<string> & leavesNames);

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
};


#endif	//_TREETOOLS_H_
