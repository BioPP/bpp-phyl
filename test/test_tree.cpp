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
#include <Bpp/Phyl/Io/Newick.h>
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
    //cout << "First test passed." << endl;
    if (!tree->hasSameTopologyAs(*tree2))
      return 1; //Error!!!
    //cout << "Second test passed." << endl;
  
    //Convert tree to string and read it again:
    string newick = TreeTemplateTools::treeToParenthesis(*tree);
    TreeTemplate<Node>* tree3 = TreeTemplateTools::parenthesisToTree(newick);
    if (!tree->hasSameTopologyAs(*tree3))
      return 1; //Error!!!
    //cout << "Third test passed." << endl;
    
    //-------------
    delete tree;
    delete tree2;
    delete tree3;
  }

  //Try to parse a string:
  TreeTemplate<Node>* tree4 = TreeTemplateTools::parenthesisToTree("((A:1,B:2):3,C:4);");
  cout << TreeTemplateTools::treeToParenthesis(*tree4) << endl;
  delete tree4;

  TreeTemplate<Node>* tree5 = TreeTemplateTools::parenthesisToTree("((A:1,B:2):3,C:4):5;");
  cout << TreeTemplateTools::treeToParenthesis(*tree5) << endl;
  delete tree5;

  Newick tReader;
  istringstream iss6("((A,B),C);");
  TreeTemplate<Node>* tree6 = tReader.read(iss6);
  cout << TreeTemplateTools::treeToParenthesis(*tree6) << endl;
  delete tree6;
  
  istringstream iss7("((A:1,B:2):3,C:4):5;");
  TreeTemplate<Node>* tree7 = tReader.read(iss7);
  cout << TreeTemplateTools::treeToParenthesis(*tree7) << endl;
  delete tree7;

  istringstream iss8("((A,B)aa,C)2;");
  tReader.enableExtendedBootstrapProperty("ESS");
  TreeTemplate<Node>* tree8 = tReader.read(iss8);
  cout << TreeTemplateTools::treeToParenthesis(*tree8) << endl;
  vector<int> ids = tree8->getNodesId();
  for (size_t i = 0; i < ids.size(); ++i) {
    cout << "Node " << ids[i] << ":" << endl;
    vector<string> nodePpt = tree8->getNode(ids[i])->getNodePropertyNames();
    for (size_t j = 0; j < nodePpt.size(); ++j)
      if (tree8->getNode(ids[i])->hasNodeProperty(nodePpt[j]))
        cout << "N: " << nodePpt[j] << "=" << dynamic_cast<BppString*>(tree8->getNode(ids[i])->getNodeProperty(nodePpt[j]))->toSTL() << endl;
    vector<string> branchPpt = tree8->getNode(ids[i])->getBranchPropertyNames();
    for (size_t j = 0; j < branchPpt.size(); ++j)
      if (tree8->getNode(ids[i])->hasBranchProperty(branchPpt[j]))
        cout << "B: " << branchPpt[j] << "=" << dynamic_cast<BppString*>(tree8->getNode(ids[i])->getBranchProperty(branchPpt[j]))->toSTL() << endl;
  }
  delete tree8;

  //Test file parsing:
  TreeTemplate<Node>* tree9 = TreeTemplateTools::getRandomTree(leaves, true);
  Newick tWriter;
  tWriter.write(*tree9, "tmp_tree.dnd");
  Tree* test = tReader.read("tmp_tree.dnd");
  if (!TreeTools::haveSameTopology(*tree9, *test))
    return 1;
  cout << "Newick I/O ok." << endl;

  //Multiple trees:
  vector<Tree *> trees;
  for (unsigned int i = 0; i < 100; ++i) {
    trees.push_back(TreeTemplateTools::getRandomTree(leaves, true));
  }
  tWriter.write(trees, "tmp_trees.dnd");

  vector<Tree *> trees2;
  tReader.read("tmp_trees.dnd", trees2);

  for (unsigned int i = 0; i < 100; ++i) {
    if (!TreeTools::haveSameTopology(*trees[i], *trees2[i]))
    {
      cerr << "Tree " << i << " failed to write and/or read!" << endl;
      return 1;
    }
  }
  cout << "Newick multiple I/O ok." << endl;

  for (unsigned int i = 0; i < 100; ++i) {
    delete trees[i];
    delete trees2[i];
  }

  //Try newick read on non-file:
  cout << "Testing parsing a directory..." << endl;
  try {
    Tree* tmp = tReader.read("test/");
    cerr << "Arg, reading on directory should fail!" << endl;
    if (tmp != NULL) {
      cerr << "Output of read on directory is not NULL!" << endl;
    }
    return 1;
  } catch(Exception& ex) {
    cout << "Ok, reading on directory throws exception!" << endl;
  }

  cout << "Testing parsing a directory for multiple trees..." << endl;
  try {
    vector<Tree*> treesTmp;
    tReader.read("test/", treesTmp);
    if (treesTmp.size() != 0) {
      cerr << "Output of multiple read on directory is not 0!" << endl;
      return 1;
    } else {
      cout << "Ok, reading on directory returns a vector of size 0!" << endl;
    }
  } catch(Exception& ex) {
    cout << "Error, no exception should be thrown here!" << endl;
  }

  //Now try some weird cases, to see if we handle them properly:
  //single node tree:
  cout << "Testing a tree with a node of degree 2:" << endl;
  TreeTemplate<Node>* weird1 = TreeTemplateTools::parenthesisToTree("((A:1):2.0,B);");
  if (weird1->getNodes().size() != 4) {
    cout << "Error, tree has " << weird1->getNodes().size() << " node(s) instead of 4!" << endl;
    VectorTools::print(weird1->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird1) << endl;
  delete weird1;

  cout << "Testing a tree with a node of degree 2, without branch length:" << endl;
  TreeTemplate<Node>* weird2 = TreeTemplateTools::parenthesisToTree("((A),B);");
  if (weird2->getNodes().size() != 4) {
    cout << "Error, tree has " << weird2->getNodes().size() << " node(s) instead of 4!" << endl;
    VectorTools::print(weird2->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird2) << endl;
  delete weird2;

  cout << "Testing a tree with several single nodes:" << endl;
  TreeTemplate<Node>* weird3 = TreeTemplateTools::parenthesisToTree("((((((A)):1)):3),B);");
  if (weird3->getNodes().size() != 8) {
    cout << "Error, tree has " << weird3->getNodes().size() << " node(s) instead of 8!" << endl;
    VectorTools::print(weird3->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird3) << endl;
  delete weird3;

  cout << "Testing a tree with a single leaf:" << endl;
  TreeTemplate<Node>* weird4 = TreeTemplateTools::parenthesisToTree("(A:1.0);");
  if (weird4->getNodes().size() != 2) {
    cout << "Error, tree has " << weird4->getNodes().size() << " node(s) instead of 2!" << endl;
    VectorTools::print(weird4->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird4) << endl;
  delete weird4;

  cout << "Testing a tree with a single node:" << endl;
  TreeTemplate<Node>* weird5 = TreeTemplateTools::parenthesisToTree("((A:1.0));");
  if (weird5->getNodes().size() != 3) {
    cout << "Error, tree has " << weird5->getNodes().size() << " node(s) instead of 3!" << endl;
    VectorTools::print(weird5->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird5) << endl;
  delete weird5;

  cout << "Testing a tree with a single node and branch lengths:" << endl;
  TreeTemplate<Node>* weird6 = TreeTemplateTools::parenthesisToTree("((A:1.0):2.0);");
  if (weird6->getNodes().size() != 3) {
    cout << "Error, tree has " << weird6->getNodes().size() << " node(s) instead of 3!" << endl;
    VectorTools::print(weird6->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird6) << endl;
  delete weird6;




  return 0;
}
