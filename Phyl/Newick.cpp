//
// File: Newick.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thu Oct 23 15:35:03 2003
//

#include "Newick.h"
#include "TreeTools.h"

// From Utils:
#include <Utils/TextTools.h>

// From the STL:
#include <iostream>
#include <fstream>

using namespace std;

/******************************************************************************/

Newick::Newick(bool allowComments): _allowComments(allowComments) {}

Newick::~Newick() {}
	
/******************************************************************************/

const string Newick::getFormatName() const { return "Newick"; }

/******************************************************************************/

const string Newick::getFormatDescription() const {
	return string("New hampshire parenthesis format. ") +
		"See http://evolution.genetics.washington.edu/phylip/newicktree.html for more info.";
}

/******************************************************************************/

Tree * Newick::read(const string & path) const throw (Exception) {
	// Checking the existence of specified file
	ifstream file(path.c_str(), ios::in);
	if (! file) { throw IOException ("Newick::read : failed to open file"); }
	
	//We concatenate all line in file till we reach the ending semi colon:
	string temp, description;// Initialization
	// Main loop : for all file lines
	while (! file.eof()) {
		getline(file, temp, '\n');  // Copy current line in temporary string
		int index = temp.find(";");
		if(index >= 0 && index < (int)temp.size()) {
			description += temp.substr(0, index + 1);
			break;
		} else description += temp;
	}
	file.close();
	if(_allowComments) description = TextTools::removeSubstrings(description, '[', ']');
	return TreeTools::parenthesisToTree(description);
}

/******************************************************************************/

void Newick::write(const Tree & tree, const string & path, bool overwrite) const throw (Exception) {
	// Open file in specified mode
	ofstream file(path.c_str(), overwrite ? (ios::out) : (ios::out|ios::app));

	// Checking the existence of specified file, and possibility to open it in write mode
	if (! file) { throw IOException ("Newick::write : failed to open file"); }
	file << TreeTools::treeToParenthesis(tree);
	file.close();
}

/******************************************************************************/
