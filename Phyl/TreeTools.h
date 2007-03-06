//
// File: TreeTools.h
// Created by:  Julien Dutheil
// Created on: Wed Aug  6 13:45:28 2003
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

#ifndef _TREETOOLS_H_
#define _TREETOOLS_H_

#include "TreeExceptions.h"
#include "Node.h"
#include "DistanceMatrix.h"
#include "Tree.h"

// From Utils:
#include <Utils/Exceptions.h>

// From NumCalc:
#include <NumCalc/VectorTools.h>

/**
 * @brief Inner class for parsing strings in Newick format.
 */
class NodeTokenizer
{
  protected:
    vector<string> tokens;
    mutable unsigned int currentPosition;
  
  public:
    NodeTokenizer(const string & description) throw (IOException): tokens(), currentPosition(0)
    {
      //cout << "NODETOKENIZER: " << description << endl;
      unsigned int tokCount = 0;
      int parCount = 0;
      unsigned int i;
      for(i = 0; i < description.size(); i++)
      {
        if(description[i] == '(') parCount++; //Another open parenthesis
        if(description[i] == ')') parCount--; //Another close parenthesis
        if(parCount < 0) throw IOException("Invalid tree description: closing parenthesis with no opening one, in " + description);
        if(description[i] == ',' && parCount == 0)
        {
          //New token found:
          //cout << "NODETOENIZER: NEWTOKEN " << description.substr(tokCount, i - tokCount - 1) << endl;
          tokens.push_back(description.substr(tokCount, i - tokCount));
          tokCount = i + 1;
        }          
      }
      //Add last token:
      //cout << "NODETOKENIZER: NEWTOKEN " << description.substr(tokCount) << endl;
      tokens.push_back(description.substr(tokCount));
      
      currentPosition = 0;
    }
    
  public:
    string next() const
    {
      string s = tokens[currentPosition];
      currentPosition++;
      return s;
    }
    bool hasNext() const
    { 
      return currentPosition < tokens.size();
    }
};

/**
 * @brief Generic utilitary methods dealing with trees.
 *
 * These methods work with all Tree object.
 * However, depending on the tree implementation, they may not be the most efficient.
 *
 * @see TreeTemplateTools
 */
class TreeTools
{
  public:
    TreeTools() {}
    virtual ~TreeTools() {}
  
  public:

    /**
     * @name Retrieve topology information
     *
     * @{
     */
    
    /**
     * @brief Retrieve all leaves from a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @return A vector with the ids of all leaves in the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static vector<int> getLeavesId(const Tree & tree, int nodeId) throw (NodeNotFoundException);

    /**
     * @brief Retrieve all leaves from a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @param leaves A vector with the ids of all leaves in the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static void getLeavesId(const Tree & tree, int nodeId, vector<int> & leaves) throw (NodeNotFoundException);
 
    /**
     * @brief Get the id of a leaf given its name in a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @param name The name of the node.
     * @return The id of the node.
     * @throw NodeNotFoundException If the node is not found.
     */
    static int getLeafId(const Tree & tree, int nodeId, const string & name) throw (NodeNotFoundException);

    /**
     * @brief Get the id of a leaf given its name in a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @param name The name of the node.
     * @param id The id of the node.
     * @throw NodeNotFoundException If the node is not found.
     */
    static void searchLeaf(const Tree & tree, int nodeId, const string & name, int * & id) throw (NodeNotFoundException);

    /**
     * @brief Get a vector of ancestor nodes between to nodes.
     *
     * @param tree The tree to use.
     * @param nodeId1 Id of first node.
     * @param nodeId2 Id of second node.
     * @param includeAncestor Tell if the common ancestor must be included in the vector.
     * @return A vector of ancestor nodes ids.
     * @throw NodeNotFoundException If the node is not found.
     */
    static vector<int> getPathBetweenAnyTwoNodes(const Tree & tree, int nodeId1, int nodeId2, bool includeAncestor = true) throw (NodeNotFoundException);
  
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
     *
     * @param tree The tree.
     * @param nodeId The id of node defining the subtree.
     * @return The depth of the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static unsigned int getDepth(const Tree & tree, int nodeId) throw (NodeNotFoundException);

    /**
     * @brief Get the height of the subtree defined by node 'node', i.e. the maximum
     * distance between leaves and the root of the subtree.
     *
     * The distance do not include the branch length of the subtree root node.
     * The height of a leaf is hence 0.
     *
     * @param tree The tree.
     * @param nodeId The id of node defining the subtree.
     * @return The height of the subtree.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If a branch length is lacking.
     */ 
    static double getHeight(const Tree & tree, int nodeId) throw (NodeNotFoundException,NodeException);
    /** @} */

    /**
     * @name Act on branch lengths.
     *
     * @{
     */
    
     /**
     * @brief Get all the branch lengths of a subtree.
     *
     * @param tree The tree.
     * @param nodeId The node defining the subtree.
     * @return A vector with all branch lengths.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If a branch length is lacking.
     */
    static Vdouble getBranchLengths(const Tree & tree, int nodeId) throw (NodeNotFoundException,NodeException);
    
    /**
     * @brief Get the total length (sum of all branch lengths) of a subtree.
     *
     * @param tree The tree.
     * @param nodeId The node defining the subtree.
     * @return The total length of the subtree.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If a branch length is lacking.
     */
    static double getTotalLength(const Tree & tree, int nodeId) throw (NodeNotFoundException,NodeException);

    /**
     * @brief Set all the branch lengths of a subtree.
     *
     * @param tree The tree.
     * @param nodeId The node defining the subtree.
     * @param brLen The branch length to apply.
     * @throw NodeNotFoundException If the node is not found.
     */
    static void setBranchLengths(Tree & tree, int nodeId, double brLen) throw (NodeNotFoundException);
          
    /**
     * @brief Give a length to branches that don't have one in a subtree.
     *
     * @param tree The tree.
     * @param nodeId The node defining the subtree.
     * @param brLen The branch length to apply.
     * @throw NodeNotFoundException If the node is not found.
     */
    static void setVoidBranchLengths(Tree & tree, int nodeId, double brLen) throw (NodeNotFoundException);
        
    /**
     * @brief Scale a given tree.
     *
     * Multiply all branch lengths by a given factor.
     *
     * @param tree The tree.
     * @param nodeId The node defining the subtree.
     * @param factor The factor to multiply all branch lengths with.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If a branch length is lacking.
     */
    static void scaleTree(Tree & tree, int nodeId, double factor) throw (NodeNotFoundException,NodeException);

    /**
     * @brief Grafen's method to initialize branch lengths.
     *
     * Each height of the node (toatl distance from the leaves) is set equal to the number of
     * leaf nodes for the corresponding subtrees - 1 for inner nodes, 0 for leaves.
     * 
     * If the tree already has branch lengths, they will be ignored.
     * 
     * Reference:
     * Grafen A. The phylogenetic regression. Philos Trans R Soc Lond B Biol Sci. 1989; 326(1233):119-57
     * 
     * @param tree The tree.
     */ 
    static void initBranchLengthsGrafen(Tree & tree);
    
    /**
     * @brief Compute branch lengths using Grafen's method.
     *
     * The 'height' of each node is devided by the total height of the tree, and the ratio is raised at power 'rho'.
     * A value of rho=0 hence returns a star tree.
     *
     * Reference:
     * Grafen A. The phylogenetic regression. Philos Trans R Soc Lond B Biol Sci. 1989; 326(1233):119-57
     * 
     * @param tree The tree to use.
     * @param power The rho parameter.
     * @param init Tell if the height must be initialized by calling the initBranchLengthsGrafen() method.
     *             Otherwise use branch lengths.
     * @throw NodeException If init=false and one branch length is lacking.
     */
    static void computeBranchLengthsGrafen(Tree & tree, double power=1, bool init=true) throw (NodeException);
   
  private:
    static unsigned int initBranchLengthsGrafen(Tree & tree, int nodeId) throw (NodeNotFoundException);
    static void computeBranchLengthsGrafen(Tree & tree, int nodeId, double power, double total, double & height, double & heightRaised) throw (NodeNotFoundException,NodeException);

  public:
    /**
     * @brief Modify a tree's branch lengths to make a clock tree.
     *
     * The height of each node is set to the mean height of all son nodes.
     * This may however lead to negative branch lengths, since the mean heigth
     * may be inferior to one of the son heights, due to short branch lengths.
     * If the 'noneg' is set to yes, the mean height is checked against all son
     * heights. If it is inferior to one of the son heights, the maximum son
     * height is used instead. This results in a multifurcation.
     * 
     * This method is recursive and will be applied on all sons nodes.
     * 
     * @param tree The tree to use.
     * @param nodeId The node defining the subtree.
     * @param noneg Tell if the correction for non negative branch lengths must be used.
     * @return The modified height of the node.
     * @throw NodeNotFoundException If the node is not found.
     * @throw NodeException If one branch length is lacking.
     */
    static double convertToClockTree(Tree & tree, int nodeId, bool noneg=false) throw (NodeNotFoundException,NodeException);
    
  public:
    /**
     * @brief Get the total distance between two nodes.
     *
     * Sum all branch lengths between two nodes.
     *
     * @param tree The tree to consider.
     * @param nodeId1 First node id.
     * @param nodeId2 Second node id.
     * @return The sum of all branch lengths between the two nodes.
     * @throw NodeNotFoundException If the node is not found.
     */
    static double getDistanceBetweenAnyTwoNodes(const Tree & tree, int nodeId1, int nodeId2) throw (NodeNotFoundException);
    
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
    /** @} */

    /**
     * @name Conversion tools.
     *
     * Convert from Newick standard tree description.
     * The description is for a node, and hence is to be surrounded with
     * parenthesis. ex: (A:0.001, (B:0.001, C:0.02)90:0.005)50:0.0005
     *
     * @{
     */

    struct Element
    {
      string content;
      string length;
      string bootstrap;
    };

    static Element getElement(string elt) throw (IOException);

    
    /**
     * @brief Get the parenthesis description of a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @param writeId Tells if node ids must be printed.
     *                This will overwrite bootstrap values if there are ones.
     *                Leaves id will be added to the leave names, separated by a '_' character.
     * @return A string in the parenthesis format.
     * @throw NodeNotFoundException If the node is not found.
     */
    static string nodeToParenthesis(const Tree & tree, int nodeId, bool writeId = false) throw (NodeNotFoundException);

    /**
     * @brief Get the parenthesis description of a subtree.
     *
     * @param tree The tree
     * @param nodeId The node defining the subtree.
     * @param bootstrap Tell is bootstrap values must be writen.
     * If so, the content of the property with name TreeTools::BOOTSTRAP will be written as bootstrap value.
     * The property should be a Number<double> object.
     * Otherwise, the content of the property with name 'propertyName' will be written.
     * In this later case, the property should be a String object.
     * @param propertyName The name of the property to use. Only used if bootstrap = false.
     * @return A string in the parenthesis format.
     * @throw NodeNotFoundException If the node is not found.
     */
    static string nodeToParenthesis(const Tree & tree, int nodeId, bool bootstrap, const string & propertyName) throw (NodeNotFoundException);

    /**
     * @brief Get the parenthesis description of a tree.
     *
     * @param tree The tree to convert.
     * @param writeId Tells if node ids must be printed.
     *                This will overwrite bootstrap values if there are ones.
     *                Leaves id will be added to the leave names, separated by a '_' character.
     * @return A string in the parenthesis format.
     */
    static string treeToParenthesis(const Tree & tree, bool writeId = false);
    
    /**
     * @brief Get the parenthesis description of a tree.
     *
     * @param tree The tree to convert.
     * @param bootstrap Tell is bootstrap values must be writen.
     * If so, the content of the property with name TreeTools::BOOTSTRAP will be written as bootstrap value.
     * The property should be a Number<double> object.
     * Otherwise, the content of the property with name 'propertyName' will be written.
     * In this later case, the property should be a String object.
     * @param propertyName The name of the property to use. Only used if bootstrap = false.
     * @return A string in the parenthesis format.
     */
    static string treeToParenthesis(const Tree & tree, bool bootstrap, const string & propertyName);
    
    /** @} */

    /**
     * @name Deal with identifiers
     *
     * @{
     */

    /**
     * @brief Retrieve all son nodes from a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @return A vector of ids of each son node in the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static vector<int> getNodesId(const Tree & tree, int nodeId) throw (NodeNotFoundException);

    /**
     * @brief Retrieve all son nodes from a subtree.
     *
     * @param tree The tree
     * @param nodeId The id of node defining the subtree.
     * @param nodes A vector of ids of each son node in the subtree.
     * @throw NodeNotFoundException If the node is not found.
     */
    static void getNodesId(const Tree & tree, int nodeId, vector<int> & nodes) throw (NodeNotFoundException);

    /**
     * @brief Get the maximum identifier used in a (sub)tree.
     *
     * This is a recursive method.
     *
     * @param tree The tree to check.
     * @param id The identifier of the subtree from which the recursion will be performed.
     * Use id=tree.getRootNodeId() to search for the whole tree.
     * @return The identifier number with maximum value.
     */
    static int getMaxId(const Tree & tree, int id);

    /**
     * @brief Get the minimum positive non-used identifier in a (sub)tree.
     *
     * This method uses the recursive method getNodesId, and then sort the ids.
     *
     * @param tree The tree to check.
     * @param id The identifier of the subtree from which the recursion will be performed.
     * Use id=tree.getRootNodeId() to search for the whole tree.
     * @return A non-used identifier number.
     */
    static int getMPNUId(const Tree & tree, int id);

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

#endif  //_TREETOOLS_H_

