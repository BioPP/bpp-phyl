// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_IO_PHYLIPDISTANCEMATRIXFORMAT_H
#define BPP_PHYL_IO_PHYLIPDISTANCEMATRIXFORMAT_H


#include "IoDistanceMatrix.h"

namespace bpp
{
/**
 * @brief Distance matrix I/O in Phylip format.
 *
 * Entry names must be 10 characters long. If 'extended' is set to true, then
 * entry names can be of any size, and should be separated from the data by at least two spaces.
 * Names should therefor not contian more than one consecutive space.
 */
class PhylipDistanceMatrixFormat :
  public AbstractIDistanceMatrix,
  public AbstractODistanceMatrix
{
private:
  bool extended_;

public:
  PhylipDistanceMatrixFormat(bool extended = false) : extended_(extended) {}
  virtual ~PhylipDistanceMatrixFormat() {}

public:
  const std::string getFormatName() const { return "Phylip"; }

  const std::string getFormatDescription() const {  return "Multiline space-delimited columns."; }

  std::unique_ptr<DistanceMatrix> readDistanceMatrix(const std::string& path) const
  {
    return AbstractIDistanceMatrix::readDistanceMatrix(path);
  }

  std::unique_ptr<DistanceMatrix> readDistanceMatrix(std::istream& in) const;

  void writeDistanceMatrix(const DistanceMatrix& dist, const std::string& path, bool overwrite = true) const
  {
    AbstractODistanceMatrix::writeDistanceMatrix(dist, path, overwrite);
  }

  void writeDistanceMatrix(const DistanceMatrix& dist, std::ostream& out) const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_IO_PHYLIPDISTANCEMATRIXFORMAT_H
