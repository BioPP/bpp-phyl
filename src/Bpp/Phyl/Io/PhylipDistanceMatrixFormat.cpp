//
// File: PhylipDistanceMatrixFormat.cpp
// Created by: Julien Dutheil
// Created on: Wed Jun 08 15:57 2005
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

#include "PhylipDistanceMatrixFormat.h"

#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>

// From SeqLib:
#include <Bpp/Seq/DistanceMatrix.h>

using namespace bpp;

// From the STL:
#include <iomanip>

using namespace std;

DistanceMatrix * PhylipDistanceMatrixFormat::read(istream& in) const throw (Exception)
{
	string s = FileTools::getNextLine(in);
	// the size of the matrix:
	unsigned int n = TextTools::fromString<unsigned int>(s);
	DistanceMatrix * dist = new DistanceMatrix(n);
	unsigned int rowNumber = 0;
	unsigned int colNumber = 0;
	s = FileTools::getNextLine(in);
	while (in)
  {
		if (colNumber == 0)
    { // New row
      if (extended_) {
        size_t pos = s.find("  ");
        if (pos == string::npos)
          throw Exception("PhylipDistanceMatrixFormat::read. Bad format, probably not 'extended' Phylip.");
        dist->setName(rowNumber, s.substr(0, pos));
        s = s.substr(pos + 2);
      } else {
        dist->setName(rowNumber, s.substr(0, 10));
        s = s.substr(11);
      }
    }
		StringTokenizer st(s, "\t ");
		for (; colNumber < n && st.hasMoreToken(); colNumber++)
    {
			double d = TextTools::fromString<double>(st.nextToken());
			(* dist)(rowNumber, colNumber) = d;
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

void PhylipDistanceMatrixFormat::write(const DistanceMatrix& dist, ostream& out) const throw (Exception)
{
	size_t n = dist.size();
	out << "   " << n << endl;
  size_t offset = 10;
  if (extended_) {
    offset = 0;
    for (size_t i = 0; i < n; ++i) {
      size_t s = dist.getName(i).size();
      if (s > offset) offset = s;
    }
  }
	for (unsigned int i = 0; i < n; i++)
  {
    out << TextTools::resizeRight(dist.getName(i), offset, ' ');
    if (extended_) {
      out << "  ";
    } else {
      out << " ";
    }
    for (unsigned int j = 0; j < n; j++) {
      if (j > 0) out << " ";
		  out << setprecision(8) << dist(i, j); 
    }
		out << endl;
	}
}

