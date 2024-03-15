// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Tree/TreeTemplateTools.h>
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
    auto tree = TreeTemplateTools::getRandomTree(leaves, true);
    auto tree2 = new TreeTemplate<Node>(*tree);
    if (!tree->hasSameTopologyAs(*tree2))
      return 1; //Error!!!
    tree2->getRootNode()->swap(0,1);
    //cout << "First test passed." << endl;
    if (!tree->hasSameTopologyAs(*tree2))
      return 1; //Error!!!
    //cout << "Second test passed." << endl;
  
    //Convert tree to string and read it again:
    string newick = TreeTemplateTools::treeToParenthesis(*tree);
    auto tree3 = TreeTemplateTools::parenthesisToTree(newick, true, TreeTools::BOOTSTRAP, false, false);
    if (!tree->hasSameTopologyAs(*tree3))
      return 1; //Error!!!
    //cout << "Third test passed." << endl;
    
    //-------------
  }

  //Try to parse a string:
  auto tree4 = TreeTemplateTools::parenthesisToTree("((A:1,B:2):3,C:4);");
  cout << TreeTemplateTools::treeToParenthesis(*tree4) << endl;

  auto tree5 = TreeTemplateTools::parenthesisToTree("((A:1,B:2):3,C:4):5;");
  cout << TreeTemplateTools::treeToParenthesis(*tree5) << endl;

  Newick tReader;
  istringstream iss6("((A,B),C);");
  auto tree6 = tReader.readTree(iss6);
  cout << TreeTemplateTools::treeToParenthesis(*tree6) << endl;
  
  istringstream iss7("((A:1,B:2):3,C:4):5;");
  auto tree7 = tReader.readTree(iss7);
  cout << TreeTemplateTools::treeToParenthesis(*tree7) << endl;

  istringstream iss8("((A:1,B:2)80:3,C:4)2:5;");
  auto tree8 = tReader.readTreeTemplate(iss8);
  cout << TreeTemplateTools::treeToParenthesis(*tree8) << endl;
  vector<int> ids = tree8->getNodesId();
  for (size_t i = 0; i < ids.size(); ++i) {
    cout << "Node " << ids[i] << ":" << endl;
    if (tree8->getNode(ids[i])->hasBranchProperty(TreeTools::BOOTSTRAP))
      cout << "N: BOOTSTRAP=" << dynamic_cast<Number<double>*>(tree8->getNode(ids[i])->getBranchProperty(TreeTools::BOOTSTRAP))->getValue() << endl;
    vector<string> branchPpt = tree8->getNode(ids[i])->getBranchPropertyNames();
  }

  istringstream iss9("((A,B)aa,C)2;");
  tReader.enableExtendedBootstrapProperty("ESS");
  auto tree9 = tReader.readTreeTemplate(iss9);
  cout << TreeTemplateTools::treeToParenthesis(*tree9) << endl;
  ids = tree9->getNodesId();
  for (size_t i = 0; i < ids.size(); ++i) {
    cout << "Node " << ids[i] << ":" << endl;
    vector<string> nodePpt = tree9->getNode(ids[i])->getNodePropertyNames();
    for (size_t j = 0; j < nodePpt.size(); ++j)
      if (tree9->getNode(ids[i])->hasNodeProperty(nodePpt[j]))
        cout << "N: " << nodePpt[j] << "=" << dynamic_cast<BppString*>(tree9->getNode(ids[i])->getNodeProperty(nodePpt[j]))->toSTL() << endl;
    vector<string> branchPpt = tree9->getNode(ids[i])->getBranchPropertyNames();
    for (size_t j = 0; j < branchPpt.size(); ++j)
      if (tree9->getNode(ids[i])->hasBranchProperty(branchPpt[j]))
        cout << "B: " << branchPpt[j] << "=" << dynamic_cast<BppString*>(tree9->getNode(ids[i])->getBranchProperty(branchPpt[j]))->toSTL() << endl;
  }

  //Test file parsing:
  auto tree10 = TreeTemplateTools::getRandomTree(leaves, true);
  Newick tWriter;
  tWriter.writeTree(*tree10, "tmp_tree.dnd", true);
  auto test = tReader.readTree("tmp_tree.dnd");
  if (!TreeTools::haveSameTopology(*tree10, *test))
    return 1;
  cout << "Newick I/O ok." << endl;

  //Multiple trees:
  vector<const Tree *> trees;
  for (unsigned int i = 0; i < 100; ++i) {
    trees.push_back(TreeTemplateTools::getRandomTree(leaves, true).release());
  }
  tWriter.writeTrees(trees, "tmp_trees.dnd", true);

  vector<unique_ptr<Tree>> trees2;
  tReader.readTrees("tmp_trees.dnd", trees2);

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
  }

  //Try newick read on non-file:
  cout << "Testing parsing a directory..." << endl;
  try {
    auto tmp = tReader.readTree("test/");
    cerr << "Arg, reading on directory should fail!" << endl;
    if (tmp) {
      cerr << "Output of read on directory is not NULL!" << endl;
    }
    return 1;
  } catch (Exception& ex) {
    cout << "Ok, reading on directory throws exception!" << endl;
  }

  cout << "Testing parsing a directory for multiple trees..." << endl;
  try {
    vector<unique_ptr<Tree>> treesTmp;
    tReader.readTrees("test/", treesTmp);
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
  auto weird1 = TreeTemplateTools::parenthesisToTree("((A:1):2.0,B);");
  if (weird1->getNodes().size() != 4) {
    cout << "Error, tree has " << weird1->getNodes().size() << " node(s) instead of 4!" << endl;
    VectorTools::print(weird1->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird1) << endl;

  cout << "Testing a tree with a node of degree 2, without branch length:" << endl;
  auto weird2 = TreeTemplateTools::parenthesisToTree("((A),B);");
  if (weird2->getNodes().size() != 4) {
    cout << "Error, tree has " << weird2->getNodes().size() << " node(s) instead of 4!" << endl;
    VectorTools::print(weird2->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird2) << endl;

  cout << "Testing a tree with several single nodes:" << endl;
  auto weird3 = TreeTemplateTools::parenthesisToTree("((((((A)):1)):3),B);");
  if (weird3->getNodes().size() != 8) {
    cout << "Error, tree has " << weird3->getNodes().size() << " node(s) instead of 8!" << endl;
    VectorTools::print(weird3->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird3) << endl;

  cout << "Testing a tree with a single leaf:" << endl;
  auto weird4 = TreeTemplateTools::parenthesisToTree("(A:1.0);");
  if (weird4->getNodes().size() != 2) {
    cout << "Error, tree has " << weird4->getNodes().size() << " node(s) instead of 2!" << endl;
    VectorTools::print(weird4->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird4) << endl;

  cout << "Testing a tree with a single node:" << endl;
  auto weird5 = TreeTemplateTools::parenthesisToTree("((A:1.0));");
  if (weird5->getNodes().size() != 3) {
    cout << "Error, tree has " << weird5->getNodes().size() << " node(s) instead of 3!" << endl;
    VectorTools::print(weird5->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird5) << endl;

  cout << "Testing a tree with a single node and branch lengths:" << endl;
  auto weird6 = TreeTemplateTools::parenthesisToTree("((A:1.0):2.0);");
  if (weird6->getNodes().size() != 3) {
    cout << "Error, tree has " << weird6->getNodes().size() << " node(s) instead of 3!" << endl;
    VectorTools::print(weird6->getLeavesNames());
    return 1;
  }
  cout << TreeTemplateTools::treeToParenthesis(*weird6) << endl;
 
  return 0;
}
