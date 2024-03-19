// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Phyl/Tree/TreeTemplate.h>
// #include <Bpp/Phyl/Tree/TreeIterator.h>
#include <Bpp/Phyl/Tree/TreeTemplateTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <string>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;

void giveNamesToInternalNodes(Tree* tree)
{
  TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(tree);
  vector<Node*> nodes = ttree->getNodes();
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    if (!nodes[i]->hasName())
      nodes[i]->setName("N" + TextTools::toString(nodes[i]->getId()));
  }
}

int main()
{
  // parse a string:
/*  TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("(((S1:1,S2:1):3,S3:2):1,(S4:1,S5:1):2);");
   cout << TreeTemplateTools::treeToParenthesis(*tree) << endl;
   giveNamesToInternalNodes(tree); // give internal names to nodes in post-order
   TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(tree);

   // iterate in preorder
   string expectedOrder1[9] = {"N8", "N4", "N2", "S1", "S2", "S3", "N7", "S4", "S5"};
   PreOrderTreeIterator* treeIt1 = new PreOrderTreeIterator(*ttree);
   int counter = 0;
   for (const Node* node = treeIt1->begin(); node != treeIt1->end(); node = treeIt1->next()) {
      if (node->getName().compare(expectedOrder1[counter]) != 0)
      {
        cerr << "Preorder traversion failed at step " << counter << ": returned " << node->getName() << " instead of " << expectedOrder1[counter] << endl;
        return 1;
      }
      counter += 1;
   }
   delete(treeIt1);

   // iterate in inorder
   string expectedOrder2[9] = {"S1", "N2", "S2", "N4", "S3", "N8", "S4", "N7", "S5"};
   InOrderTreeIterator* treeIt2 = new InOrderTreeIterator(*ttree);
   counter = 0;
   for (const Node* node = treeIt2->begin(); node != treeIt2->end(); node = treeIt2->next()) {
      if (node->getName().compare(expectedOrder2[counter]) != 0)
      {
        cerr << "Inorder traversion failed at step " << counter << ": returned " << node->getName() << " instead of " << expectedOrder2[counter] << endl;
        return 1;
      }
      counter += 1;
   }
   delete(treeIt2);

   // iterate in postorder
   string expectedOrder3[9] = {"S1", "S2", "N2", "S3", "N4", "S4", "S5", "N7", "N8"};
   PostOrderTreeIterator* treeIt3 = new PostOrderTreeIterator(*ttree);
   counter = 0;
   for (const Node* node = treeIt3->begin(); node != treeIt3->end(); node = treeIt3->next()) {
      if (node->getName().compare(expectedOrder3[counter]) != 0)
      {
        cerr << "Postorder traversion failed at step " << counter << ": returned " << node->getName() << " instead of " << expectedOrder3[counter] << endl;
        return 1;
      }
      counter += 1;
   }
   delete(treeIt3);

   delete(tree);
 */return 0;
}
