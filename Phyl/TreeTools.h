//
// File: TreeTools.h
// Created by:  <@bogdanof>
// Created on: Wed Aug  6 13:45:28 2003
//

#ifndef _TREETOOLS_H_
#define _TREETOOLS_H_

//#include "Tree.h"
#include "TreeExceptions.h" // Include a declaration of classes Node and Tree<Node>.

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
		static vector<const Node *> getLeaves(const Node & node);

		/**
		 * @brief Retrieve all leaves from a subtree.
		 *
		 * @param node The node that defines the subtree.
		 * @return A vector of pointers toward each leaf in the subtree.
		 */
		static vector<Node *> getLeaves(Node & node);

		/**
		 * @brief Retrieve all son nodes from a subtree.
		 *
		 * @param node The node that defines the subtree.
		 * @return A vector of pointers toward each son node in the subtree.
		 */
		static vector<const Node *> getNodes(const Node & node);

		/**
		 * @brief Retrieve all son nodes from a subtree.
		 *
		 * @param node The node that defines the subtree.
		 * @return A vector of pointers toward each son node in the subtree.
		 */
		static vector<Node *> getNodes(Node & node);

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
		vector<Node *> getPathBetweenAnyTwoNodes(Node & node1, Node & node2);

		template<class N>
		static const N * getSon(const N & node, unsigned int i)
		{
			return dynamic_cast<const N *>(node.getSon(i));
		}
		
		template<class N>
		static N * getSon(N & node, unsigned int i)
		{
			return dynamic_cast<N *>(node.getSon(i));
		}
		
		template<class N>
		static const N * getFather(const N & node)
		{
			return dynamic_cast<const N *>(node.getFather());
		}
		
		template<class N>
		static N * getFather(N & node)
		{
			return dynamic_cast<N *>(node.getFather());
		}
		
		template<class N>
		static N * cloneSubtree(const N & node) 
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
