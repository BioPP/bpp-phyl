// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Tree/TreeTemplateTools.h>
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
  auto tree = TreeTemplateTools::getRandomTree(leaves, true);

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
  nhxParser.writeTree(*tree, out);
  out.close();
  auto tree2 = nhxParser.readTree("randomTree.nhx");
  ofstream out2("randomTree2.nhx", ios::out);
  nhxParser.writeTree(*tree2, out2);
  out2.close();

  return 0;
}
