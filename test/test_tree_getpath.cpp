//
// File: test_tree_getpath.cpp
// Created by: Julien Dutheil
// Created on: Mon Oct 20 09:56 2014
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

#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Tree/TreeTemplateTools.h>
#include <string>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree ("(A:0.1, (B :0.2, C:0.4):0.1);" );
  vector <int> Leaves_select = tree->getLeavesId();

  Node* node1 = tree->getNode(Leaves_select[0]);
  Node* node2 = tree->getNode(Leaves_select[1]);
  Node* node3 = tree->getNode(Leaves_select[2]);

  vector<Node*> vecNode = TreeTemplateTools::getPathBetweenAnyTwoNodes(*node1, *node2, true);
  cout << "Id node1 " << node1->getId() << endl;
  cout << "Id node2 " << node2->getId() << endl;
  cout << "Ids of path:" << endl;
  for (size_t i = 0; i < vecNode.size(); i++){
    cout << vecNode[i]->getId() << endl;
  }
  if (vecNode[0]->getId() != 0) return -1;
  if (vecNode[1]->getId() != 4) return -1;
  if (vecNode[2]->getId() != 3) return -1;
  if (vecNode[3]->getId() != 1) return -1;

  vecNode = TreeTemplateTools::getPathBetweenAnyTwoNodes(*node1, *node3, true);
  cout << "Id node1 " << node1->getId() << endl;
  cout << "Id node3 " << node3->getId() << endl;
  cout << "Ids of path: " << endl;
  for (size_t i = 0; i < vecNode.size(); i++){
    cout << vecNode[i]->getId() << endl;
  }
  if (vecNode[0]->getId() != 0) return -1;
  if (vecNode[1]->getId() != 4) return -1;
  if (vecNode[2]->getId() != 3) return -1;
  if (vecNode[3]->getId() != 2) return -1;
 
  return 0;
}
