// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_IODAGRAPH_H
#define BPP_PHYL_IO_IODAGRAPH_H


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
 * @brief General interface for DAG I/O.
 */
class IODAG :
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
  IODAG() {}
  virtual ~IODAG() {}

public:
  virtual const std::string getDataType() const { return "DAG"; }
};

/**
 * @brief General interface for DAG readers.
 */
class IPhyloDAG :
  public virtual IODAG
{
public:
  IPhyloDAG() {}
  virtual ~IPhyloDAG() {}

public:
  /**
   * @brief Read a DAG from a file.
   *
   * @param path The file path.
   * @return a new PhyloDAG object.
   */

  virtual std::unique_ptr<PhyloDAG> readPhyloDAG(const std::string& path) const = 0;

  /**
   * @brief Read a DAG from a stream.
   *
   * @param in The input stream.
   * @return A new PhyloDAG object.
   */

  virtual std::unique_ptr<PhyloDAG> readPhyloDAG(std::istream& in) const = 0;
};


/**
 * @brief General interface for DAG writers.
 */

class OPhyloDAG :
  public virtual IODAG
{
public:
  OPhyloDAG() {}
  virtual ~OPhyloDAG() {}

public:
  /**
   * @brief Write a DAG to a file.
   *
   * @param dag a PhyloDAG object.
   * @param path The file path.
   * @param overwrite Tell if existing file must be overwritten.
   * Otherwise append to the file.
   * @throw Exception If an error occurred.
   */
  virtual void writePhyloDAG(const PhyloDAG& dag, const std::string& path, bool overwrite) const = 0;

  /**
   * @brief Write a PhyloDAG to a stream.
   *
   * @param dag A PhyloDAG object.
   * @param out The output stream.
   * @throw Exception If an error occurred.
   */
  virtual void writePhyloDAG(const PhyloDAG& dag, std::ostream& out) const =  0;
};

/**
 * @brief Partial implementation of the IDAG interface.
 */

class AbstractIPhyloDAG :
  public virtual IPhyloDAG
{
public:
  AbstractIPhyloDAG() {}
  virtual ~AbstractIPhyloDAG() {}

public:
  std::unique_ptr<PhyloDAG> readPhyloDAG(std::istream& in) const override = 0;

  std::unique_ptr<PhyloDAG> readPhyloDAG(const std::string& path) const override
  {
    std::ifstream input(path.c_str(), std::ios::in);
    if (!input)
      throw IOException ("AbstractIPhyloDAG::readPhyloDAG(path): failed to read from path " + path);
    auto dag = readPhyloDAG(input);
    input.close();
    return dag;
  }

  virtual Element getElement(const std::string& elt) const
  {
    return Element();
  }
};

/**
 * @brief Partial implementation of the ODAG interface.
 */
class AbstractOPhyloDAG :
  public virtual OPhyloDAG
{
public:
  AbstractOPhyloDAG() {}
  virtual ~AbstractOPhyloDAG() {}

public:
  virtual void writePhyloDAG(const PhyloDAG& dag, std::ostream& out) const = 0;
  virtual void writePhyloDAG(const PhyloDAG& dag, const std::string& path, bool overwrite) const
  {
    try
    {
      // Open file in specified mode

      std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));
      if (!output)
        throw Exception("Problem opening file " + path + "in writePhyloDAG.");
      writePhyloDAG(dag, output);
      output.close();
    }
    catch (IOException& e)
    {
      std::stringstream ss;
      ss << e.what() << "\nProblem writing dag to file " << path << "\n Is the file path correct and do \
you have the proper authorizations? ";
      throw (IOException ( ss.str() ) );
    }
  }
};


/**
 * @brief General interface for multiple trees readers.
 */
class IMultiPhyloDAG :
  public virtual IODAG
{
public:
  IMultiPhyloDAG() {}
  virtual ~IMultiPhyloDAG() {}

public:
  /**
   * @brief Read dags from a file.
   *
   * @param path The file path.
   * @param dags The output dags vector.
   * @throw Exception If an error occurred.
   */
  virtual void readPhyloDAGs(const std::string& path, std::vector<std::unique_ptr<PhyloDAG>>& dags) const = 0;

  /**
   * @brief Read dags from a stream.
   *
   * @param in The input stream.
   * @param dags The output dags container.
   * @throw Exception If an error occurred.
   */
  virtual void readPhyloDAGs(std::istream& in, std::vector<std::unique_ptr<PhyloDAG>>& dags) const = 0;
};

/**
 * @brief General interface for tree writers.
 */
class OMultiPhyloDAG :
  public virtual IODAG
{
public:
  OMultiPhyloDAG() {}
  virtual ~OMultiPhyloDAG() {}

public:
  /**
   * @brief Write dags to a file.
   *
   * @param dags A vector of dag objects.
   * @param path The file path.
   * @param overwrite Tell if existing file must be overwritten.
   * Otherwise append to the file.
   */
  virtual void writePhyloDAGs(
      const std::vector<const PhyloDAG*>& dags,
      const std::string& path,
      bool overwrite) const = 0;

  /**
   * @brief Write dags to a stream.
   *
   * @param dags A vector of dag objects.
   * @param out The output stream.
   * @throw Exception If an error occurred.
   */
  virtual void writePhyloDAGs(
      const std::vector<const PhyloDAG*>& dags,
      std::ostream& out) const = 0;
};

/**
 * @brief Partial implementation of the IMultiDAG interface.
 */
class AbstractIMultiPhyloDAG :
  public virtual IMultiPhyloDAG
{
public:
  AbstractIMultiPhyloDAG() {}
  virtual ~AbstractIMultiPhyloDAG() {}

public:
  virtual void readPhyloDAGs(std::istream& in, std::vector<std::unique_ptr<PhyloDAG>>& dags) const override = 0;

  virtual void readPhyloDAGs(const std::string& path, std::vector<std::unique_ptr<PhyloDAG>>& dags) const override
  {
    std::ifstream input(path.c_str(), std::ios::in);
    readPhyloDAGs(input, dags);
    input.close();
  }
};

/**
 * @brief Partial implementation of the ODAG interface.
 */
class AbstractOMultiPhyloDAG :
  public virtual OMultiPhyloDAG
{
public:
  AbstractOMultiPhyloDAG() {}
  virtual ~AbstractOMultiPhyloDAG() {}

public:
  virtual void writePhyloDAGs(const std::vector<const PhyloDAG*>& dags, std::ostream& out) const = 0;

  virtual void writePhyloDAGs(const std::vector<const PhyloDAG*>& dags, const std::string& path, bool overwrite) const

  {
    // Open file in specified mode
    std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));
    writePhyloDAGs(dags, output);
    output.close();
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_IODAGRAPH_H
