// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Tree/TreeTemplateTools.h>
#include <string>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  auto tree = TreeTemplateTools::parenthesisToTree ("(A:0.1, (B :0.2, C:0.4):0.1);" );
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
