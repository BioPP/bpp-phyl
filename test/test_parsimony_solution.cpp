// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>
#include <iostream>

using namespace bpp;
using namespace std;

void giveNamesToInternalNodes(TreeTemplate<Node>& ttree)
{
  vector<Node*> nodes = ttree.getNodes();
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    if (!nodes[i]->hasName())
      nodes[i]->setName("N" + TextTools::toString(nodes[i]->getId()));
  }
}

int main()
{
  try
  {
    // process tree
    shared_ptr<TreeTemplate<Node>> ttree = TreeTemplateTools::parenthesisToTree("(((((((S1:1,S2:1):1,S3:2):1,S4:3):1,S5:4):1,S6:5):1,S7:6):1,S8:7);");
    cout << TreeTemplateTools::treeToParenthesis(*ttree) << endl;

    // process character data
    shared_ptr<const Alphabet> alphabet = make_shared<BinaryAlphabet>();
    auto sites = make_shared<VectorSiteContainer>(alphabet);
    auto seq1 = make_unique<Sequence>("S1", "0", alphabet);
    sites->addSequence("S1", seq1);
    auto seq2 = make_unique<Sequence>("S2", "1", alphabet);
    sites->addSequence("S2", seq2);
    auto seq3 = make_unique<Sequence>("S3", "0", alphabet);
    sites->addSequence("S3", seq3);
    auto seq4 = make_unique<Sequence>("S4", "1", alphabet);
    sites->addSequence("S4", seq4);
    auto seq5 = make_unique<Sequence>("S5", "1", alphabet);
    sites->addSequence("S5", seq5);
    auto seq6 = make_unique<Sequence>("S6", "1", alphabet);
    sites->addSequence("S6", seq6);
    auto seq7 = make_unique<Sequence>("S7", "0", alphabet);
    sites->addSequence("S7", seq7);
    auto seq8 = make_unique<Sequence>("S8", "0", alphabet);
    sites->addSequence("S8", seq8);

    // compute the maxmum parsimony score and solution according to ACCTRAN approach
    auto mpData = make_shared<DRTreeParsimonyScore>(ttree,  sites);

    // make sure the score is 3
    if (mpData->getScore() != 3)
    {
      cerr << "Error! The compacted maximum parsimony score is incorrect" << endl;
      return 1;
    }

    // make sure the solution is: (((((((S1{0},S2{1})N2{0},S3{0})N4{0},S4{1})N6{1},S5{1})N8{1},S6{1})N10{1},S7{0})N12{0},S8{0})N14{0}
    mpData->computeSolution();
    auto solution = make_shared<TreeTemplate<Node>>(mpData->tree());
    giveNamesToInternalNodes(*solution); // give internal names to nodes in post-order
    map<string, int> nodeToState;
    nodeToState["S1"] = 0;
    nodeToState["S2"] = 1;
    nodeToState["N2"] = 0;
    nodeToState["S3"] = 0;
    nodeToState["N4"] = 0;
    nodeToState["S4"] = 1;
    nodeToState["N6"] = 1;
    nodeToState["S5"] = 1;
    nodeToState["N8"] = 1;
    nodeToState["S6"] = 1;
    nodeToState["N10"] = 1;
    nodeToState["S7"] = 0;
    nodeToState["N12"] = 0;
    nodeToState["S8"] = 0;
    nodeToState["N14"] = 0;
    vector<Node*> nodes = solution->getNodes();
    string nodeName;
    int nodeState;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
      nodeName = nodes[i]->getName();
      nodeState = (int) mpData->getNodeState(nodes[i]);
      if (nodeState != nodeToState[nodeName])
      {
        cerr << "Error! assignment of state in node " << nodeName << " is " << nodeState << " instead of " << nodeToState[nodeName] << endl;
        return 1;
      }
    }
    cout << "all good :)" << endl;
  }
  catch (Exception& ex)
  {
    cerr << ex.what() << endl;
    return 1;
  }

  return 0;
}
