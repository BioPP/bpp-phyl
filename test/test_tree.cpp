//
// File: test_tree.cpp
// Created by: Julien Dutheil
// Created on: Sun Nov 14 10:20 2010
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

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

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <string>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  //Get some leaf names:
  vector<string> leaves(100);
  for (size_t i = 0; i < leaves.size(); ++i)
    leaves[i] = "leaf" + TextTools::toString(i);
  
  for (unsigned int j = 0; j < 1000; ++j) {
    //Generate a random tree, without branch lengths:
    TreeTemplate<Node>* tree = TreeTemplateTools::getRandomTree(leaves, true);
    TreeTemplate<Node>* tree2 = new TreeTemplate<Node>(*tree);
    if (!tree->hasSameTopologyAs(*tree2))
      return 1; //Error!!!
    tree2->getRootNode()->swap(0,1);
    cout << "First test passed." << endl;
    if (!tree->hasSameTopologyAs(*tree2))
      return 1; //Error!!!
    cout << "Second test passed." << endl;
  
    //Convert tree to string and read it again:
    string newick = TreeTemplateTools::treeToParenthesis(*tree);
    TreeTemplate<Node>* tree3 = TreeTemplateTools::parenthesisToTree(newick);
    if (!tree->hasSameTopologyAs(*tree3))
      return 1; //Error!!!
    cout << "Third test passed." << endl;
    
    //-------------
    delete tree;
    delete tree2;
    delete tree3;
  }

  return 0;
}
