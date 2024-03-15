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
  //Get some leaf names:
  vector<string> leaves(100);

  for (size_t i = 0; i < leaves.size(); ++i)
    leaves[i] = "leaf" + TextTools::toString(i);
  
  //Testing rerooting:
  //Get a random topology:
  unique_ptr< TreeTemplate<Node> > tr(TreeTemplateTools::getRandomTree(leaves, true));
  //Set random branch lengths:

  vector<Node*> nodes = tr->getNodes();
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (nodes[i]->hasFather())
      nodes[i]->setDistanceToFather(RandomTools::giveRandomNumberBetweenZeroAndEntry(1.));
  }
  
  double totalLen = TreeTools::getTotalLength(*tr.get(), tr->getRootId(), false);
  ApplicationTools::displayResult("Total length", totalLen);

  for (unsigned int i = 0; i < 100; ++i) {
    size_t pos = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(100);
    tr->rootAt(nodes[pos]);
    double l = TreeTools::getTotalLength(*tr.get(), tr->getRootId(), false);
    if ((l - totalLen) / totalLen > 0.00000001) {
      cerr << "Error, rerooting gave incorrect branch lengths :(: " << l << " vs. " << totalLen << endl;
      return 1;
    }
  }

  return 0;
}
