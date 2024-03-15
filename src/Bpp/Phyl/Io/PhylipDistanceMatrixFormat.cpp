// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>

#include "PhylipDistanceMatrixFormat.h"

// From SeqLib:
#include <Bpp/Seq/DistanceMatrix.h>

using namespace bpp;

// From the STL:
#include <iomanip>

using namespace std;

unique_ptr<DistanceMatrix> PhylipDistanceMatrixFormat::readDistanceMatrix(istream& in) const
{
  string s = FileTools::getNextLine(in);
  // the size of the matrix:
  unsigned int n = TextTools::fromString<unsigned int>(s);
  auto dist = make_unique<DistanceMatrix>(n);
  unsigned int rowNumber = 0;
  unsigned int colNumber = 0;
  s = FileTools::getNextLine(in);
  while (in)
  {
    if (colNumber == 0)
    { // New row
      if (extended_)
      {
        size_t pos = s.find("  ");
        if (pos == string::npos)
          throw Exception("PhylipDistanceMatrixFormat::read. Bad format, probably not 'extended' Phylip.");
        dist->setName(rowNumber, s.substr(0, pos));
        s = s.substr(pos + 2);
      }
      else
      {
        dist->setName(rowNumber, s.substr(0, 10));
        s = s.substr(11);
      }
    }
    StringTokenizer st(s, "\t ");
    for ( ; colNumber < n && st.hasMoreToken(); colNumber++)
    {
      double d = TextTools::fromString<double>(st.nextToken());
      (*dist)(rowNumber, colNumber) = d;
    }
    if (colNumber == n)
    {
      colNumber = 0;
      rowNumber++;
    }
    s = FileTools::getNextLine(in);
  }
  return dist;
}

void PhylipDistanceMatrixFormat::writeDistanceMatrix(const DistanceMatrix& dist, ostream& out) const
{
  size_t n = dist.size();
  out << "   " << n << endl;
  size_t offset = 10;
  if (extended_)
  {
    offset = 0;
    for (size_t i = 0; i < n; ++i)
    {
      size_t s = dist.getName(i).size();
      if (s > offset)
        offset = s;
    }
  }
  for (unsigned int i = 0; i < n; i++)
  {
    out << TextTools::resizeRight(dist.getName(i), offset, ' ');
    if (extended_)
    {
      out << "  ";
    }
    else
    {
      out << " ";
    }
    for (unsigned int j = 0; j < n; j++)
    {
      if (j > 0)
        out << " ";
      out << setprecision(8) << dist(i, j);
    }
    out << endl;
  }
}
