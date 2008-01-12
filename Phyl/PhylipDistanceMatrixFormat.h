//
// File: PhylipDistanceMatrixFormat.h
// Created by: Julien Dutheil
// Created on: Wed Jun 08 15:57 2005
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

#ifndef _PHYLIPDISTANCEMATRIXFORMAT_H_
#define _PHYLIPDISTANCEMATRIXFORMAT_H_

#include "IODistanceMatrix.h"

namespace bpp
{

/**
 * @brief Distance matrix I/O in Phylip format.
 */
class PhylipDistanceMatrixFormat:
  public virtual AbstractIDistanceMatrix,
  public virtual AbstractODistanceMatrix
{
	public:
		PhylipDistanceMatrixFormat() {}
		virtual ~PhylipDistanceMatrixFormat() {}

	public:
		const string getFormatName() const { return "Phylip"; }

		const string getFormatDescription() const {	return "Multiline space-delimited columns."; }
		DistanceMatrix * read(const string & path) const throw (Exception)
		{
			return AbstractIDistanceMatrix::read(path);
		}
		DistanceMatrix * read(istream & in) const throw (Exception);
		
		void write(const DistanceMatrix & dist, const string & path, bool overwrite = true) const throw (Exception)
		{
			AbstractODistanceMatrix::write(dist, path, overwrite);
		}
		void write(const DistanceMatrix & dist, ostream & out) const throw (Exception);

};

} //end of namespace bpp.

#endif //_PHYLIPDISTANCEMATRIXFORMAT_H_

