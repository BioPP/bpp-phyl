//
// File: Newick.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thu Oct 23 15:35:03 2003
//

#ifndef _NEWICK_H_
#define _NEWICK_H_

#include "IOTree.h"

class Newick: public ITree, public OTree
{
	protected:
		bool _allowComments;
	
	public:
		Newick(bool allowComments = false);
		virtual ~Newick();
	
	public:
		const string getFormatName() const;
		const string getFormatDescription() const;		
		Tree * read(const string & path) const throw (Exception);
		void write(const Tree & tree, const string & path, bool overwrite = true) const throw (Exception);
};


#endif	//_NEWICK_H_
