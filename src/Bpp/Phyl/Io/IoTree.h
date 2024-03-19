// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_IOTREE_H
#define BPP_PHYL_IO_IOTREE_H


#include "../Tree/PhyloTree.h"
#include "../Tree/Tree.h"

// From the STL:
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>

#include <Bpp/Exceptions.h>
#include <Bpp/Io/IoFormat.h>

namespace bpp
{
/**
 * @brief General interface for tree I/O.
 */
class IOTree :
  public virtual IOFormat
{
protected:
  struct Element
  {
public:
    std::string content;
    std::string length;
    std::string annotation;
    bool isLeaf;

public:
    Element() : content(),
      length(),
      annotation(),
      isLeaf(false) {}
  };

public:
  IOTree() {}
  virtual ~IOTree() {}

public:
  virtual const std::string getDataType() const { return "Tree"; }
};

/**
 * @brief General interface for tree readers.
 */
class ITree :
  public virtual IOTree
{
  /*
   * @brief Basic element for branch description.
   *
   */

public:
  ITree() {}
  virtual ~ITree() {}

public:
  /**
   * @brief Read a tree from a file.
   *
   * @param path The file path.
   * @return A new tree object.
   * @throw Exception If an error occured.
   */
  virtual std::unique_ptr<Tree> readTree(const std::string& path) const = 0;

  /**
   * @brief Read a tree from a stream.
   *
   * @param in The input stream.
   * @return A new tree object.
   * @throw Exception If an error occured.
   */
  virtual std::unique_ptr<Tree> readTree(std::istream& in) const = 0;
};

/**
 * @brief General interface for tree readers.
 */
class IPhyloTree :
  public virtual IOTree
{
public:
  IPhyloTree() {}
  virtual ~IPhyloTree() {}

public:
  /**
   * @brief Read a tree from a file.
   *
   * @param path The file path.
   * @return A new tree object.
   * @throw Exception If an error occured.
   */
  virtual std::unique_ptr<PhyloTree> readPhyloTree(const std::string& path) const = 0;

  /**
   * @brief Read a tree from a stream.
   *
   * @param in The input stream.
   * @return A new tree object.
   * @throw Exception If an error occured.
   */
  virtual std::unique_ptr<PhyloTree> readPhyloTree(std::istream& in) const = 0;
};


/**
 * @brief General interface for tree writers.
 */
class OTree :
  public virtual IOTree
{
public:
  OTree() {}
  virtual ~OTree() {}

public:
  /**
   * @brief Write a tree to a file.
   *
   * @param tree A tree object.
   * @param path The file path.
   * @param overwrite Tell if existing file must be overwritten.
   * Otherwise append to the file.
   * @throw Exception If an error occured.
   */
  virtual void writeTree(const Tree& tree, const std::string& path, bool overwrite) const = 0;

  /**
   * @brief Write a tree to a stream.
   *
   * @param tree A tree object.
   * @param out The output stream.
   * @throw Exception If an error occured.
   */
  virtual void writeTree(const Tree& tree, std::ostream& out) const = 0;
};

/**
 * @brief General interface for tree writers.
 */
class OPhyloTree :
  public virtual IOTree
{
public:
  OPhyloTree() {}
  virtual ~OPhyloTree() {}

public:
  /**
   * @brief Write a tree to a file.
   *
   * @param tree A tree object.
   * @param path The file path.
   * @param overwrite Tell if existing file must be overwritten.
   * Otherwise append to the file.
   * @throw Exception If an error occured.
   */
  virtual void writePhyloTree(const PhyloTree& tree, const std::string& path, bool overwrite) const = 0;

  /**
   * @brief Write a tree to a stream.
   *
   * @param tree A tree object.
   * @param out The output stream.
   * @throw Exception If an error occured.
   */
  virtual void writePhyloTree(const PhyloTree& tree, std::ostream& out) const =  0;
};

/**
 * @brief Partial implementation of the ITree interface.
 */
class AbstractITree :
  public virtual ITree
{
public:
  AbstractITree() {}
  virtual ~AbstractITree() {}

public:
  std::unique_ptr<Tree> readTree(std::istream& in) const override
  {
    auto tree = readTreeTemplate(in);
    return std::unique_ptr<Tree>(tree.release());
  }

  std::unique_ptr<Tree> readTree(const std::string& path) const override
  {
    std::ifstream input(path.c_str(), std::ios::in);
    auto tree = readTree(input);
    input.close();
    return tree;
  }

  virtual std::unique_ptr<TreeTemplate<Node>> readTreeTemplate(std::istream& in) const = 0;

  virtual std::unique_ptr<TreeTemplate<Node>> readTreeTemplate(const std::string& path) const
  {
    std::ifstream input(path.c_str(), std::ios::in);
    auto tree = readTreeTemplate(input);
    input.close();
    return tree;
  }

  virtual Element getElement(const std::string& elt) const
  {
    return Element();
  }
};

class AbstractIPhyloTree :
  public virtual IPhyloTree
{
public:
  AbstractIPhyloTree() {}
  virtual ~AbstractIPhyloTree() {}

public:
  std::unique_ptr<PhyloTree> readPhyloTree(std::istream& in) const override = 0;

  std::unique_ptr<PhyloTree> readPhyloTree(const std::string& path) const override
  {
    std::ifstream input(path.c_str(), std::ios::in);
    if (!input)
      throw IOException ("AbstractIPhyloTree::readPhyloTree(path): failed to read from path " + path);
    auto tree = readPhyloTree(input);
    input.close();
    return tree;
  }

  Element getElement(const std::string& elt) const
  {
    return Element();
  }
};

/**
 * @brief Partial implementation of the OTree interface.
 */
class AbstractOTree :
  public virtual OTree
{
public:
  AbstractOTree() {}
  virtual ~AbstractOTree() {}

public:
  virtual void writeTree(const Tree& tree, std::ostream& out) const = 0;
  virtual void writeTree(const Tree& tree, const std::string& path, bool overwrite) const
  {
    try
    {
      // Open file in specified mode

      std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));
      if (!output)
        throw Exception("Problem opening file " + path + "in write Tree.");
      writeTree(tree, output);
      output.close();
    }
    catch (IOException& e)
    {
      std::stringstream ss;
      ss << e.what() << "\nProblem writing tree to file " << path << "\n Is the file path correct and do \
you have the proper authorizations? ";
      throw (IOException ( ss.str() ) );
    }
  }
};

/**
 * @brief Partial implementation of the OTree interface.
 */
class AbstractOPhyloTree :
  public virtual OPhyloTree
{
public:
  AbstractOPhyloTree() {}
  virtual ~AbstractOPhyloTree() {}

public:
  virtual void writePhyloTree(const PhyloTree& tree, std::ostream& out) const = 0;
  virtual void writePhyloTree(const PhyloTree& tree, const std::string& path, bool overwrite) const
  {
    try
    {
      // Open file in specified mode

      std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));
      if (!output)
        throw Exception("Problem opening file " + path + "in writePhyloTree.");
      writePhyloTree(tree, output);
      output.close();
    }
    catch (IOException& e)
    {
      std::stringstream ss;
      ss << e.what() << "\nProblem writing tree to file " << path << "\n Is the file path correct and do \
you have the proper authorizations? ";
      throw (IOException ( ss.str() ) );
    }
  }
};


/**
 * @brief General interface for multiple trees readers.
 */
class IMultiTree :
  public virtual IOTree
{
public:
  IMultiTree() {}
  virtual ~IMultiTree() {}

public:
  /**
   * @brief Read trees from a file.
   *
   * @param path The file path.
   * @param trees The output trees container.
   * @throw Exception If an error occured.
   */
  virtual void readTrees(const std::string& path, std::vector<std::unique_ptr<Tree>>& trees) const = 0;

  /**
   * @brief Read trees from a stream.
   *
   * @param in The input stream.
   * @param trees The output trees container.
   * @throw Exception If an error occured.
   */
  virtual void readTrees(std::istream& in, std::vector<std::unique_ptr<Tree>>& trees) const = 0;
};

class IMultiPhyloTree :
  public virtual IOTree
{
public:
  IMultiPhyloTree() {}
  virtual ~IMultiPhyloTree() {}

public:
  /**
   * @brief Read trees from a file.
   *
   * @param path The file path.
   * @param trees The output trees container.
   * @throw Exception If an error occured.
   */
  virtual void readPhyloTrees(const std::string& path, std::vector<std::unique_ptr<PhyloTree>>& trees) const = 0;

  /**
   * @brief Read trees from a stream.
   *
   * @param in The input stream.
   * @param trees The output trees container.
   * @throw Exception If an error occured.
   */
  virtual void readPhyloTrees(std::istream& in, std::vector<std::unique_ptr<PhyloTree>>& trees) const = 0;
};

/**
 * @brief General interface for tree writers.
 */
class OMultiTree :
  public virtual IOTree
{
public:
  OMultiTree() {}
  virtual ~OMultiTree() {}

public:
  /**
   * @brief Write trees to a file.
   *
   * @param trees A vector of tree objects.
   * @param path The file path.
   * @param overwrite Tell if existing file must be overwritten.
   * Otherwise append to the file.
   * @throw Exception If an error occured.
   */
  virtual void writeTrees(
    const std::vector<const Tree*>& trees,
    const std::string& path,
    bool overwrite) const = 0;

  /**
   * @brief Write trees to a stream.
   *
   * @param trees A vector of tree objects.
   * @param out The output stream.
   * @throw Exception If an error occured.
   */
  virtual void writeTrees(
    const std::vector<const Tree*>& trees,
    std::ostream& out) const = 0;
};

/**
 * @brief General interface for tree writers.
 */
class OMultiPhyloTree :
  public virtual IOTree
{
public:
  OMultiPhyloTree() {}
  virtual ~OMultiPhyloTree() {}

public:
  /**
   * @brief Write trees to a file.
   *
   * @param trees A vector of tree objects.
   * @param path The file path.
   * @param overwrite Tell if existing file must be overwritten.
   * Otherwise append to the file.
   * @throw Exception If an error occured.
   */
  virtual void writePhyloTrees(
    const std::vector<const PhyloTree*>& trees,
    const std::string& path,
    bool overwrite) const = 0;

  /**
   * @brief Write trees to a stream.
   *
   * @param trees A vector of tree objects.
   * @param out The output stream.
   * @throw Exception If an error occured.
   */
  virtual void writePhyloTrees(
    const std::vector<const PhyloTree*>& trees,
    std::ostream& out) const = 0;
};

/**
 * @brief Partial implementation of the IMultiTree interface.
 */
class AbstractIMultiTree :
  public virtual IMultiTree
{
public:
  AbstractIMultiTree() {}
  virtual ~AbstractIMultiTree() {}

public:
  virtual void readTrees(std::istream& in, std::vector<std::unique_ptr<Tree>>& trees) const override = 0;

  virtual void readTrees(const std::string& path, std::vector<std::unique_ptr<Tree>>& trees) const override
  {
    std::ifstream input(path.c_str(), std::ios::in);
    readTrees(input, trees);
    input.close();
  }
};

/**
 * @brief Partial implementation of the IMultiTree interface.
 */
class AbstractIMultiPhyloTree :
  public virtual IMultiPhyloTree
{
public:
  AbstractIMultiPhyloTree() {}
  virtual ~AbstractIMultiPhyloTree() {}

public:
  virtual void readPhyloTrees(std::istream& in, std::vector<std::unique_ptr<PhyloTree>>& trees) const override = 0;

  virtual void readPhyloTrees(const std::string& path, std::vector<std::unique_ptr<PhyloTree>>& trees) const override
  {
    std::ifstream input(path.c_str(), std::ios::in);
    readPhyloTrees(input, trees);
    input.close();
  }
};

/**
 * @brief Partial implementation of the OTree interface.
 */
class AbstractOMultiTree :
  public virtual OMultiTree
{
public:
  AbstractOMultiTree() {}
  virtual ~AbstractOMultiTree() {}

public:
  virtual void writeTrees(const std::vector<const Tree*>& trees, std::ostream& out) const = 0;
  virtual void writeTrees(const std::vector<const Tree*>& trees, const std::string& path, bool overwrite) const

  {
    // Open file in specified mode
    std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));
    writeTrees(trees, output);
    output.close();
  }
};

/**
 * @brief Partial implementation of the OTree interface.
 */
class AbstractOMultiPhyloTree :
  public virtual OMultiPhyloTree
{
public:
  AbstractOMultiPhyloTree() {}
  virtual ~AbstractOMultiPhyloTree() {}

public:
  virtual void writePhyloTrees(const std::vector<const PhyloTree*>& trees, std::ostream& out) const = 0;

  virtual void writePhyloTrees(const std::vector<const PhyloTree*>& trees, const std::string& path, bool overwrite) const

  {
    // Open file in specified mode
    std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));
    writePhyloTrees(trees, output);
    output.close();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_IOTREE_H
