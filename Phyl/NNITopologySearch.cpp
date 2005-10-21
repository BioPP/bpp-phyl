//
// File: NNITopologySearch.h
// Created by: Julien Dutheil
// Created on: Wed Oct 12 10:52 2005
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

#include "NNITopologySearch.h"

//From Utils:
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>

//From NumCalc:
#include <NumCalc/VectorTools.h>
using namespace VectorFunctions;

const string NNITopologySearch::FAST   = "Fast";
const string NNITopologySearch::BETTER = "Better";
const string NNITopologySearch::PHYML  = "PhyML";

void NNITopologySearch::search() throw (Exception)
{
	     if(_algorithm == FAST)   searchFast();
	else if(_algorithm == BETTER) searchBetter();
	else if(_algorithm == PHYML)  searchPhyML();
  else throw Exception("Unknown NNI algorithm: " + _algorithm + ".\n");
}

void NNITopologySearch::searchFast() throw (Exception)
{
	TreeTemplate<Node> * tree = dynamic_cast<TreeTemplate<Node> *>(_searchableTree->getTree());
	vector<Node *> nodes = tree->getNodes();

	bool test = true;
	do { 
		
		vector<Node *> nodesSub = nodes;
		for(unsigned int i = nodesSub.size(); i>0; i--) {// !!! must not reach i==0 because of unsigned int
			if(!(nodesSub[i-1]->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove root node.	
			else if(!(nodesSub[i-1]->getFather()->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove son of root node.	
		}
		
		// Test all NNIs:
		test = false;
		for(unsigned int i = 0; !test && i < nodesSub.size(); i++) {
			Node * node = nodesSub[i];
			double diff = _searchableTree->testNNI(node->getFather(), node);
			if(_verbose >= 2) {
				ApplicationTools::displayMessage(TextTools::toString(node->getId())
						                    + "->" + TextTools::toString(node->getFather()->getId())
																+ ": " + TextTools::toString(diff));
			}
			
			if(diff < 0.) { //Good NNI found...
				if(_verbose == 1) {
					ApplicationTools::displayMessage(TextTools::toString(node->getId())
							                    + "->" + TextTools::toString(node->getFather()->getId())
																	+ ": " + TextTools::toString(diff));
				}
				_searchableTree->doNNI(node->getFather(), node);
				// Notify:
				_searchableTree->topologyChangePerformed(TopologyChangeEvent());
				test = true;
			}
		}
		if(_verbose >= 2) ApplicationTools::displayTaskDone();
		
	} while(test);
}

void NNITopologySearch::searchBetter() throw (Exception)
{
	TreeTemplate<Node> * tree = dynamic_cast<TreeTemplate<Node> *>(_searchableTree->getTree());
	vector<Node *> nodes = tree->getNodes();

	bool test = true;
	do { 
		if(_verbose >= 2) ApplicationTools::displayTask("Test all possible NNIs...");
		
		vector<Node *> nodesSub = nodes;
		for(unsigned int i = nodesSub.size(); i>0; i--) {// !!! must not reach i==0 because of unsigned int
			if(!(nodesSub[i-1]->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove root node.	
			else if(!(nodesSub[i-1]->getFather()->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove son of root node.	
		}
		
		// Test all NNIs:
		vector<Node *> improving;
		vector<double> improvement;
		if(_verbose >= 2) ApplicationTools::message << endl;
		for(unsigned int i = 0; i < nodesSub.size(); i++) {
			Node * node = nodesSub[i];
			double diff = _searchableTree->testNNI(node->getFather(), node);
			if(_verbose >= 2) {
				ApplicationTools::displayMessage(TextTools::toString(node->getId())
						                    + "->" + TextTools::toString(node->getFather()->getId())
																+ ": " + TextTools::toString(diff));
			}
			
			if(diff < 0.) {
				improving.push_back(node);
				improvement.push_back(diff);
			}
		}
		if(_verbose >= 2) ApplicationTools::displayTaskDone();
		test = improving.size() > 0;
		if(test) {
			unsigned int nodeMin = posmin(improvement);
			Node * node = improving[nodeMin];
			if(_verbose >= 1) {
				ApplicationTools::displayMessage(TextTools::toString(node->getId()) + ": " + TextTools::toString(improvement[nodeMin]));
			}
			_searchableTree->doNNI(node->getFather(), node);
			
			// Notify:
			_searchableTree->topologyChangePerformed(TopologyChangeEvent());
		}
	} while(test);
	
}

void NNITopologySearch::searchPhyML() throw (Exception)
{
	TreeTemplate<Node> * tree = dynamic_cast<TreeTemplate<Node> *>(_searchableTree->getTree());
	vector<Node *> nodes = tree->getNodes();

	bool test = true;
	do { 
		if(_verbose >= 2) ApplicationTools::displayTask("Test all possible NNIs...");
		
		vector<Node *> nodesSub = nodes;
		for(unsigned int i = nodesSub.size(); i>0; i--) {// !!! must not reach i==0 because of unsigned int
			if(!(nodesSub[i-1]->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove root node.	
			else if(!(nodesSub[i-1]->getFather()->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove son of root node.	
		}
		
		// Test all NNIs:
		vector<Node *> improving;
		vector<double> improvement;
		if(_verbose >= 2) ApplicationTools::message << endl;
		for(unsigned int i = 0; i < nodesSub.size(); i++) {
			Node * node = nodesSub[i];
			double diff = _searchableTree->testNNI(node->getFather(), node);
			if(_verbose >= 2) {
				ApplicationTools::displayMessage(TextTools::toString(node->getId())
						                    + "->" + TextTools::toString(node->getFather()->getId())
																+ ": " + TextTools::toString(diff));
			}
			
			if(diff < 0.) {
				// Must test for brother NNIs...
				bool ok = true;
				for(unsigned int j = 0; j < improving.size(); j++) {
					if(improving[j]->getFather() == node->getFather()) {
						//These are brother NNIs. We only keep the best:
						if(diff < improvement[j]) { //Replace node
							improving[j] = node;
							improvement[j] = diff;
						} //Otherwise forget about this NNI.
						ok = false;
					}
				}
				if(ok) {//No brother NNI found. We add this NNI to the list.
					improving.push_back(node);
					improvement.push_back(diff);
				}
			}
		}
		if(_verbose >= 2) ApplicationTools::displayTaskDone();
		test = improving.size() > 0;
		if(_verbose >= 1 && test)
			ApplicationTools::displayMessage(TextTools::toString<unsigned int>(improving.size())+" rearrangements:");
		for(unsigned int i = 0; i < improving.size(); i++) {
			Node * node = improving[i];
			if(_verbose >= 1) {
				ApplicationTools::displayMessage("\t" + TextTools::toString(node->getId()) + ": " + TextTools::toString(improvement[i]));
			}
			_searchableTree->doNNI(node->getFather(), node);
		}
		
		// Notify:
		_searchableTree->topologyChangePerformed(TopologyChangeEvent());
	} while(test);
}

