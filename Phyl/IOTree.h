//
// File: IOTree.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thu Oct 23 15:19:06 2003
//

#ifndef _IOTREE_H_
#define _IOTREE_H_

#include "Tree.h"

// From the STL:
#include <string>

using namespace std;

//From Utils:
#include <Utils/Exceptions.h>

class IOTree
{
	public:
		IOTree();
		virtual ~IOTree();

	public:
		virtual const string getFormatName() const = 0;
		virtual const string getFormatDescription() const = 0;		
};

class ITree: public IOTree
{
	public:
		ITree();
		virtual ~ITree();

	public:
		virtual Tree * read(const string & path) const throw (Exception) = 0;
};

class OTree: public IOTree
{
	public:
		OTree();
		virtual ~OTree();

	public:
		virtual void write(const Tree & tree, const string & path, bool overwrite) const throw (Exception) = 0;
};

#endif	//_IOTREE_H_
