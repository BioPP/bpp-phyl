//
// File: test_tree_iterator.cpp
// Created by: Keren Halabi
// Created on: Sun Apr 29 16:48 2019
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
#include <Bpp/Phyl/TreeIterator.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
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
    for (size_t i=0; i<nodes.size(); ++i) {
        if (!nodes[i]->hasName())
            nodes[i]->setName("N" + TextTools::toString(nodes[i]->getId()));
    }  
}

int main() {

  // parse a string:
  TreeTemplate<Node>* tree = TreeTemplateTools::parenthesisToTree("(((S1:1,S2:1):3,S3:2):1,(S4:1,S5:1):2);");
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
  return 0;
}
