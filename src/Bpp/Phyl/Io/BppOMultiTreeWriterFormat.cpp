//
// File: BppOMultiTreeWriterFormat.cpp
// Created by: Julien Dutheil
// Created on: Sun Jun 19th, 13:16
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 2016)

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

#include "BppOMultiTreeWriterFormat.h"
#include "Newick.h"
#include "NexusIOTree.h"
#include "Nhx.h"

#include <Bpp/Text/KeyvalTools.h>

#include <string>
#include <memory>

using namespace bpp;
using namespace std;

OMultiTree* BppOMultiTreeWriterFormat::read(const std::string& description) throw (Exception)
{
  unparsedArguments_.clear();
  string format = "";
  KeyvalTools::parseProcedure(description, format, unparsedArguments_);
  unique_ptr<OMultiTree> oTrees;
  if (format == "Newick")
  {
    bool allowComments = ApplicationTools::getBooleanParameter("allow_comments", unparsedArguments_, false, "", true, warningLevel_);
    oTrees.reset(new Newick(allowComments));
  }
  else if (format == "Nhx")
  {
    oTrees.reset(new Nhx());
  }
  else if (format == "Nexus")
  {
    oTrees.reset(new NexusIOTree());
  }
  else
  {
    throw Exception("Trees format '" + format + "' unknown.");
  }

  return oTrees.release();
}

