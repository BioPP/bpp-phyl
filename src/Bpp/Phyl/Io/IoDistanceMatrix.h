// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_IODISTANCEMATRIX_H
#define BPP_PHYL_IO_IODISTANCEMATRIX_H

#include <Bpp/Io/IoFormat.h>
#include <Bpp/Seq/DistanceMatrix.h>

// From the STL:
#include <iostream>
#include <fstream>

namespace bpp
{
/**
 * @brief General interface for distance matrix I/O.
 */
class IODistanceMatrix :
  public virtual IOFormat
{
public:
  IODistanceMatrix() {}
  virtual ~IODistanceMatrix() {}

public:
  virtual const std::string getDataType() const { return "Distance matrix"; }
};

/**
 * @brief General interface for distance matrix readers.
 */
class IDistanceMatrix :
  public virtual IODistanceMatrix
{
public:
  IDistanceMatrix() {}
  virtual ~IDistanceMatrix() {}

public:
  /**
   * @brief Read a distance matrix from a file.
   *
   * @param path The file path.
   * @return A new distance matrix object.
   * @throw Exception If an error occured.
   */
  virtual std::unique_ptr<DistanceMatrix> readDistanceMatrix(const std::string& path) const = 0;
  /**
   * @brief Read a distance matrix from a stream.
   *
   * @param in The input stream.
   * @return A new distance matrix object.
   * @throw Exception If an error occured.
   */
  virtual std::unique_ptr<DistanceMatrix> readDistanceMatrix(std::istream& in) const = 0;
};

/**
 * @brief General interface for distance matrix writers.
 */
class ODistanceMatrix :
  public virtual IODistanceMatrix
{
public:
  ODistanceMatrix() {}
  virtual ~ODistanceMatrix() {}

public:
  /**
   * @brief Write a distance matrix to a file.
   *
   * @param dist A distance matrix object.
   * @param path The file path.
   * @param overwrite Tell if existing file must be overwritten.
   * Otherwise append to the file.
   * @throw Exception If an error occured.
   */
  virtual void writeDistanceMatrix(const DistanceMatrix& dist, const std::string& path, bool overwrite) const = 0;
  /**
   * @brief Write a distance matrix to a stream.
   *
   * @param dist A distance matrix object.
   * @param out The output stream.
   * @throw Exception If an error occured.
   */
  virtual void writeDistanceMatrix(const DistanceMatrix& dist, std::ostream& out) const = 0;
};

/**
 * @brief Partial implementation of the IDistanceMatrix interface.
 */
class AbstractIDistanceMatrix :
  public virtual IDistanceMatrix
{
public:
  AbstractIDistanceMatrix() {}
  virtual ~AbstractIDistanceMatrix() {}

public:
  virtual std::unique_ptr<DistanceMatrix> readDistanceMatrix(const std::string& path) const
  {
    std::ifstream input(path.c_str(), std::ios::in);
    auto mat = readDistanceMatrix(input);
    input.close();
    return mat;
  }
  virtual std::unique_ptr<DistanceMatrix> readDistanceMatrix(std::istream& in) const = 0;
};

/**
 * @brief Partial implementation of the ODistanceMatrix interface.
 */
class AbstractODistanceMatrix :
  public virtual ODistanceMatrix
{
public:
  AbstractODistanceMatrix() {}
  virtual ~AbstractODistanceMatrix() {}

public:
  virtual void writeDistanceMatrix(const DistanceMatrix& dist, const std::string& path, bool overwrite) const
  {
    // Open file in specified mode
    std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out | std::ios::app));
    writeDistanceMatrix(dist, output);
    output.close();
  }
  virtual void writeDistanceMatrix(const DistanceMatrix& dist, std::ostream& out) const = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_IODISTANCEMATRIX_H
