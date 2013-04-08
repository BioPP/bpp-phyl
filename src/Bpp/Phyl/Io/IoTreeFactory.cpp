//
// File IOTreeFactory.cpp
// Created by: Julien Dutheil
// Created on: Tue 18/04/06 10:24
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

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

#include "IoTreeFactory.h"
#include "Newick.h"
#include "NexusIoTree.h"
#include "Nhx.h"

using namespace bpp;

const std::string IOTreeFactory::NEWICK_FORMAT = "Newick"; 
const std::string IOTreeFactory::NEXUS_FORMAT = "Nexus"; 
const std::string IOTreeFactory::NHX_FORMAT = "Nhx"; 

ITree* IOTreeFactory::createReader(const std::string& format) throw (Exception)
{
       if (format == NEWICK_FORMAT) return new Newick();
  else if (format == NEXUS_FORMAT) return new NexusIOTree();
  else if (format == NHX_FORMAT) return new Nhx();
  else throw Exception("Format " + format + " is not supported for input.");
}
  
OTree* IOTreeFactory::createWriter(const std::string& format) throw (Exception)
{
       if (format == NEWICK_FORMAT) return new Newick();
  else if (format == NEXUS_FORMAT) return new NexusIOTree();
  else if (format == NHX_FORMAT) return new Nhx();
  else throw Exception("Format " + format + " is not supported for output.");
}

