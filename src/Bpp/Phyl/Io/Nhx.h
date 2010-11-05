//
// File: Nhx.h
// Created by: Bastien Boussau
// Created on: Tu Oct 19 10:24:03 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _Nhx_H_
#define _Nhx_H_

#include "IoTree.h"
#include "../TreeTemplate.h"

namespace bpp
{

/**
 * @brief The so-called 'Nhx - New Hampshire eXtended' parenthetic format. 
 *
 * See http://www.phylosoft.org/Nhx/ for details.
 *
 * Branch lengths and node annotations are supported:
 *
 * ex:
 * <code>
 * ((A:0.1[&&NHX:S=human], B:0.15[&&NHX:S=chimp]):0.2[&&NHX:B=90:D=N:S=primates], C:0.27[&&NHX:S=mouse]);
 * </code>
 *
 * Where e.g. "S=human" means that the sequence A comes from species "human", "B=90" stands for a support value, 
 * and "D=N" means that there was no duplication at the node. Other tags are allowed, see http://www.phylosoft.org/NHX/.
 * By default, comments between "[" and "]" are removed, unless the opening bracket is followed by "&&NHX".
 *
 * Code example:
 * @code
 * #include <Phyl/Nhx.h>
 * #include <Phyl/Tree.h>
 * 
 * Nhx * NhxReader = new Nhx();
 * try {
 *   Tree * tree = nhxReader->read("MyTestTree.dnd"); // Tree in file MyTestTree.dnd
 *   cout << "Tree has " << tree->getNumberOfLeaves() << " leaves." << endl;
 * } catch (Exception e) {
 *  cout << "Error when reading tree." << endl;
 * }
 * delete tree;
 * delete nhxReader;
 * @endcode
 *
 * All node annotations are stored as node properties, with type bppString for all properties except for support values, where a Number is used.
 *
 */
class Nhx:
  public AbstractITree,
  public AbstractOTree,
  public AbstractIMultiTree,
  public AbstractOMultiTree
{
  public:
    
    /**
     * @brief Build a new Nhx reader/writer.
     *
     * Comments between hooks ('[' ']') are ignored.
     * 
     */
    Nhx() {}

    virtual ~Nhx() {}
  
  public:

  struct Element
  {
  public:
    std::string content;
    std::string length;
    std::string annotation;
  
  public:
    Element() : content(), length(), annotation() {}
  };
  
  
    /**
     * @name The IOTree interface
     *
     * @{
     */
    const std::string getFormatName() const;
    const std::string getFormatDescription() const;
    /* @} */

    /**
     * @name The ITree interface
     *
     * @{
     */    
    TreeTemplate<Node>* read(const std::string& path) const throw (Exception)
    {
      return dynamic_cast<TreeTemplate<Node>*>(AbstractITree::read(path));
    }
    
    TreeTemplate<Node>* read(std::istream& in) const throw (Exception);
    /** @} */

    /**
     * @name The OTree interface
     *
     * @{
     */
    void write(const Tree& tree, const std::string& path, bool overwrite = true) const throw (Exception)
    {
      AbstractOTree::write(tree, path, overwrite);
    }
    
    void write(const Tree& tree, std::ostream& out) const throw (Exception)
    {
      write_(tree, out);
    }
    /** @} */

    /**
     * @name The IMultiTree interface
     *
     * @{
     */
    void read(const std::string& path, std::vector<Tree*>& trees) const throw (Exception)
    {
      AbstractIMultiTree::read(path, trees);
    }
    void read(std::istream& in, std::vector<Tree*>& trees) const throw (Exception);
    /**@}*/

    /**
     * @name The OMultiTree interface
     *
     * @{
     */
    void write(const std::vector<Tree*>& trees, const std::string& path, bool overwrite = true) const throw (Exception)
    {
      AbstractOMultiTree::write(trees, path, overwrite);
    }
    void write(const std::vector<Tree*>& trees, std::ostream& out) const throw (Exception)
    {
      write_(trees, out);
    }
    /** @} */

    TreeTemplate<Node>* parenthesisToTree(const std::string& description) const throw (Exception);

    std::string treeToParenthesis(const TreeTemplate<Node>& tree) const;

  protected:
    void write_(const Tree& tree, std::ostream& out) const throw (Exception);
    
    template<class N>
    void write_(const TreeTemplate<N>& tree, std::ostream& out) const throw (Exception);
    
    void write_(const std::vector<Tree*>& trees, std::ostream& out) const throw (Exception);
    
    template<class N>
    void write_(const std::vector<TreeTemplate<N>*>& trees, std::ostream& out) const throw (Exception);

    Element getElement(const std::string& elt) const throw (IOException);

    Node* parenthesisToNode(const std::string& description) const;
  
    std::string propertiesToParenthesis(const Node& node) const;
  
    std::string nodeToParenthesis(const Node& node) const;
  
    bool setNodeProperties(Node& node, const std::string properties) const;
};

} //end of namespace bpp.

#endif  //_Nhx_H_

