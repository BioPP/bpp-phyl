//
// File: Newick.cpp
// Created by: Julien Dutheil
// Created on: Thu Oct 23 15:35:03 2003
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

#include "Newick.h"
#include "../Tree/Tree.h"
#include "../Tree/TreeTemplate.h"
#include "../Tree/TreeTemplateTools.h"

#include <Bpp/Text/TextTools.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>

using namespace std;

/******************************************************************************/

const string Newick::getFormatName() const { return "Newick"; }

/******************************************************************************/

const string Newick::getFormatDescription() const
{
  return string("New hampshire parenthesis format. ") +
    "See http://evolution.genetics.washington.edu/phylip/newicktree.html for more info.";
}

/******************************************************************************/

#if defined(NO_VIRTUAL_COV)
    Tree*
#else
    TreeTemplate<Node> * 
#endif
Newick::read(istream& in) const throw (Exception)
{
  // Checking the existence of specified file
  if (! in) { throw IOException ("Newick::read: failed to read from stream"); }
  
  //We concatenate all line in file till we reach the ending semi colon:
  string temp, description;// Initialization
  // Main loop : for all file lines
  while (getline(in, temp, '\n'))
  {
    string::size_type index = temp.find(";");
    if(index != string::npos)
    {
      description += temp.substr(0, index + 1);
      break;
    }
    else description += temp;
  }

  if (allowComments_) description = TextTools::removeSubstrings(description, '[', ']');
  if (TextTools::isEmpty(description))
    throw IOException("Newick::read: no tree was found!");
  return TreeTemplateTools::parenthesisToTree(description, useBootstrap_, bootstrapPropertyName_, false, verbose_);
}

/******************************************************************************/

void Newick::write_(const Tree& tree, ostream& out) const throw (Exception)
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (! out) { throw IOException ("Newick::writeTree: failed to write to stream"); }
  if(useBootstrap_)
  {
    out << TreeTools::treeToParenthesis(tree, writeId_);
  }
  else
  {
    out << TreeTools::treeToParenthesis(tree, false, bootstrapPropertyName_);
  }
}

/******************************************************************************/

template<class N>
void Newick::write_(const TreeTemplate<N>& tree, ostream& out) const throw (Exception)
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (! out) { throw IOException ("Newick::writeTree: failed to write to stream"); }
  if(useBootstrap_)
  {
    out << TreeTemplateTools::treeToParenthesis(tree, writeId_);
  }
  else
  {
    out << TreeTemplateTools::treeToParenthesis(tree, false, bootstrapPropertyName_);
  }
}

/******************************************************************************/

void Newick::read(istream& in, vector<Tree*>& trees) const throw (Exception)
{
  // Checking the existence of specified file
  if (! in) { throw IOException ("Newick::read: failed to read from stream"); }
  
  // Main loop : for all file lines
  string temp, description;// Initialization
  string::size_type index;
  //We concatenate all line in file till we reach the ending semi colon:
  while (getline(in, temp, '\n'))
  {
    index = temp.find(";");
    if (index != string::npos)
    {
      description += temp.substr(0, index + 1);
      if (allowComments_) description = TextTools::removeSubstrings(description, '[', ']');
      trees.push_back(TreeTemplateTools::parenthesisToTree(description, useBootstrap_, bootstrapPropertyName_, false, verbose_));
      description = temp.substr(index + 1);
    }
    else description += temp;
  }
  //In case the file is empty, the method will not add any neww tree to the vector.
}

/******************************************************************************/

void Newick::write_(const vector<const Tree*>& trees, ostream& out) const throw (Exception)
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (! out) { throw IOException ("Newick::write: failed to write to stream"); }
  for(unsigned int i = 0; i < trees.size(); i++)
  {
    if(useBootstrap_)
    {
      out << TreeTools::treeToParenthesis(*trees[i], writeId_);
    }
    else
    {
      out << TreeTools::treeToParenthesis(*trees[i], false, bootstrapPropertyName_);
    }
  }
}

/******************************************************************************/

template<class N>
void Newick::write_(const vector<TreeTemplate<N>*>& trees, ostream& out) const throw (Exception)
{
  // Checking the existence of specified file, and possibility to open it in write mode
  if (! out) { throw IOException ("Newick::write: failed to write to stream"); }
  for(unsigned int i = 0; i < trees.size(); i++)
  {
    if(useBootstrap_)
    {
      out << TreeTemplateTools::treeToParenthesis(*trees[i], writeId_);
    }
    else
    {
      out << TreeTemplateTools::treeToParenthesis(*trees[i], false, bootstrapPropertyName_);
    }
  }
}

/******************************************************************************/

