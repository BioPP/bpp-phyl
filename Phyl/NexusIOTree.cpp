//
// File: NexusIOTree.cpp
// Created by: Julien Dutheil
// Created on: Wed May 27 19:06 2009
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

#include "NexusIOTree.h"
#include "Tree.h"
#include "TreeTemplate.h"
#include "TreeTemplateTools.h"

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/FileTools.h>
#include <Utils/StringTokenizer.h>
#include <Utils/NestedStringTokenizer.h>

//From SeqLib:
#include <Seq/NexusTools.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>

using namespace std;

/******************************************************************************/

const string NexusIOTree::getFormatName() const { return "Nexus"; }

/******************************************************************************/

const string NexusIOTree::getFormatDescription() const
{
	return string("Nexus format (trees only). ");
}

/******************************************************************************/

#if defined(NO_VIRTUAL_COV)
		Tree *
#else
		TreeTemplate<Node> * 
#endif
NexusIOTree::read(istream & in) const throw (Exception)
{
	// Checking the existence of specified file
	if (! in) { throw IOException ("NexusIOTree::read(). Failed to read from stream"); }
	
  //Look for the TREES block:
  string line = "";
  while(line != "BEGIN TREES;")
  {
    if(in.eof())
      throw Exception("NexusIOTree::read(). No trees block was found.");
    line = TextTools::removeSurroundingWhiteSpaces(FileTools::getNextLine(in));
  }
  
  string cmdName = "", cmdArgs = "";
  NexusTools::getNextCommand(in, cmdName, cmdArgs, false);

  //Look for the TRANSLATE command:
  map<string, string> translation;
  bool hasTranslation = false;
  if(cmdName == "TRANSLATE")
  {
    //Parse translation:
    StringTokenizer st(cmdArgs, ",");
    while (st.hasMoreToken())
    {
      string tok = TextTools::removeSurroundingWhiteSpaces(st.nextToken());
      NestedStringTokenizer nst(tok, "'", "'", " \t");
      if (nst.numberOfRemainingTokens() != 2)
        throw Exception("NexusIOTree::read(). Unvalid translation description.");
      string name = nst.nextToken();
      string tln  = nst.nextToken();
      translation[name] = tln;
    }
    hasTranslation = true;
    NexusTools::getNextCommand(in, cmdName, cmdArgs, false);
  }

  //Now parse the tree:
  if(cmdName != "TREE")
    throw Exception("NexusIOTree::read(). Unvalid command found: " + cmdName);
  string::size_type i = cmdArgs.find("=");
  if(i == string::npos)
    throw Exception("NexusIOTree::read(). unvalid format, should be tree-name=tree-description");
  string description = cmdArgs.substr(i + 1);
	TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree(description + ";", true);

  //Now translate leaf names if there is a translation:
  if(hasTranslation)
  {
    vector<Node*> leaves = tree->getLeaves();
    for(unsigned int i = 0; i < leaves.size(); i++)
    {
      string name = leaves[i]->getName();
      if(translation.find(name) == translation.end())
      {
        throw Exception("NexusIOTree::read(). No translation was given for this leaf: " + name);
      }
      leaves[i]->setName(translation[name]);
    }
  }

  return tree;
}

/******************************************************************************/

