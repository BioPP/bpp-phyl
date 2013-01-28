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
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>
#include <iostream>

using namespace bpp;
using namespace std;

int main() {
  try {
    Newick treeReader;
    auto_ptr<Tree> tree(treeReader.read("example1.mp.dnd"));
    Phylip alnReader(false, false);
    auto_ptr<SiteContainer> sites(alnReader.readAlignment("example1.ph", &AlphabetTools::DNA_ALPHABET));

    DRTreeParsimonyScore pars(*tree, *sites, true, true);
  
    cout << "Parsimony score: " << pars.getScore() << endl;

    if (pars.getScore() != 9) return 1;
    
  } catch (Exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }  

  return 0;
}
