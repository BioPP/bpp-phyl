//
// File: test_nhx.cpp
// Created by: Julien Dutheil
// Created on: Fri Dec 31 15:43 2010
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
#include <Bpp/Phyl/Io/Nhx.h>
#include <string>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  //Get some leaf names:
  vector<string> leaves(5);
  for (size_t i = 0; i < leaves.size(); ++i)
    leaves[i] = "leaf" + TextTools::toString(i);
  
  //Generate a random tree, without branch lengths:
  TreeTemplate<Node>* tree = TreeTemplateTools::getRandomTree(leaves, true);

  //Now assign random properties:
  vector<Node*> nodes = tree->getNodes();
  for (size_t i = 0; i < nodes.size(); ++i) {
    nodes[i]->setDistanceToFather(RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0));
    nodes[i]->setNodeProperty("GN", BppString("Gene" + TextTools::toString(i)));
    nodes[i]->setNodeProperty("AC", BppString("XXXXXX"));
    nodes[i]->setBranchProperty("B", Number<double>(floor(RandomTools::giveRandomNumberBetweenZeroAndEntry(100.) + 0.5)));
    nodes[i]->setBranchProperty("W", Number<int>(static_cast<int>(floor(RandomTools::giveRandomNumberBetweenZeroAndEntry(5.) + 0.5))));
  }
  
  //Convert tree to string and read it again:
  Nhx nhxParser(true);
  ofstream out("randomTree.nhx", ios::out);
  nhxParser.write(*tree, out);
  out.close();
  TreeTemplate<Node>* tree2 = nhxParser.read("randomTree.nhx");
  ofstream out2("randomTree2.nhx", ios::out);
  nhxParser.write(*tree2, out2);
  out2.close();

  delete tree;
  delete tree2;
  return 0;
}
