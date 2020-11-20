//
// File: test_parsimony.cpp
// Created by: Julien Dutheil
// Created on: Wed 03/10 16:47 2012
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

#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>
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
  try {
 
        // process tree
        TreeTemplate<Node>* ttree = TreeTemplateTools::parenthesisToTree("(((((((S1:1,S2:1):1,S3:2):1,S4:3):1,S5:4):1,S6:5):1,S7:6):1,S8:7);");
        cout << TreeTemplateTools::treeToParenthesis(*ttree) << endl;
        Tree* tree = dynamic_cast<Tree*>(ttree); 
        
        // process character data
        const BinaryAlphabet* alphabet = new BinaryAlphabet();
        VectorSiteContainer sites(alphabet);
        sites.addSequence(BasicSequence("S1", "0", alphabet));
        sites.addSequence(BasicSequence("S2", "1", alphabet));
        sites.addSequence(BasicSequence("S3", "0", alphabet));
        sites.addSequence(BasicSequence("S4", "1", alphabet));
        sites.addSequence(BasicSequence("S5", "1", alphabet));
        sites.addSequence(BasicSequence("S6", "1", alphabet));
        sites.addSequence(BasicSequence("S7", "0", alphabet));
        sites.addSequence(BasicSequence("S8", "0", alphabet));

        // compute the maxmum parsimony score and solution according to ACCTRAN approach
        DRTreeParsimonyScore* mpData  = new DRTreeParsimonyScore(*tree,  dynamic_cast<const SiteContainer&>(sites)); 

        // make sure the score is 3
        if (mpData->getScore() != 3)
        {
            cerr << "Error! The compated maximum parsminoy score is incorrect" << endl;
            return 1;
        }

        // make sure the solution is: (((((((S1{0},S2{1})N2{0},S3{0})N4{0},S4{1})N6{1},S5{1})N8{1},S6{1})N10{1},S7{0})N12{0},S8{0})N14{0}
        
//        mpData->computeScores(); // Should be called at construction 
        Tree* solution = mpData->getTree().clone();
        giveNamesToInternalNodes(solution); // give internal names to nodes in post-order
        map<string,int> nodeToState;
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
        vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(solution))->getNodes();
        string nodeName;
        int nodeState;
        for (size_t i=0; i<nodes.size(); ++i)
        {
            nodeName = nodes[i]->getName();
            throw Exception("test_parsimony_solution: missing method getNodeState");
            // nodeState = (int) mpData->getNodeState(nodes[i]);
            if (nodeState != nodeToState[nodeName])
            {
                cerr << "Error! assignment of state in node " << nodeName << " is " << nodeState << " instead of " << nodeToState[nodeName] << endl;
                return 1;
            }
        }

        delete mpData;
        delete solution;
        delete ttree;
        
  } catch (Exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }  

  return 0;
}
