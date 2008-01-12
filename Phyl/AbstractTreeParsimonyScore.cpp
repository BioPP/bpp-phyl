//
// File: AbstractTreeParsimonyScore.cpp
// Created by: Julien Dutheil
// Created on: Thu Jul 29 18:11 2005
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

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

#include "AbstractTreeParsimonyScore.h"
#include "PatternTools.h"
#include "TreeTemplateTools.h"

// From Utils:
#include <Utils/ApplicationTools.h>

using namespace bpp;

AbstractTreeParsimonyScore::AbstractTreeParsimonyScore(
	const Tree & tree,
	const SiteContainer & data,
	bool verbose)
	throw (Exception):
    _tree(new TreeTemplate<Node>(tree)), _data(NULL), _alphabet(data.getAlphabet()), _nbStates(_alphabet->getSize())
{
	if(_tree->isRooted())
  {
		if(verbose) ApplicationTools::displayWarning("Tree has been unrooted.");
		_tree->unroot();
	}
  TreeTemplateTools::deleteBranchLengths(*_tree->getRootNode());
	
	//Sequences will be in the same order than in the tree:
	_data = PatternTools::getSequenceSubset(data, * _tree->getRootNode());
	if(_data->getNumberOfSequences() == 1) throw Exception("Error, only 1 sequence!");
	if(_data->getNumberOfSequences() == 0) throw Exception("Error, no sequence!");
	if(_data->getAlphabet()->getSize() > 20) throw Exception("Error, only alphabet with size <= 20 are supported. See the source file of AbstractTreeParsimonyScore.");
}

vector<unsigned int> AbstractTreeParsimonyScore::getScoreForEachSite() const
{
	vector<unsigned int> scores(_data->getNumberOfSites());
	for(unsigned int i = 0; i < scores.size(); i++)
  {
		scores[i] = getScoreForSite(i);
	}
	return scores;
}

