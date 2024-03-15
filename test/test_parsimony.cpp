// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  try {
    Newick treeReader;
    shared_ptr<TreeTemplate<Node>> tree = treeReader.readTreeTemplate("example1.mp.dnd");

    Phylip alnReader(false, false);
    shared_ptr<SiteContainerInterface> sites = alnReader.readAlignment("example1.ph", AlphabetTools::DNA_ALPHABET);

    DRTreeParsimonyScore pars(tree, sites, true, true);
  
    cout << "Parsimony score: " << pars.getScore() << endl;

    if (pars.getScore() != 9) return 1;
    
  } catch (Exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }  

  return 0;
}
