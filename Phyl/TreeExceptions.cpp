//
// File: TreeExceptions.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov  3 17:04:46 2003
//

#include "TreeExceptions.h"
#include "Tree.h"

/******************************************************************************
 *                            Tree exceptions:                                *
 ******************************************************************************/

NodeException::NodeException(const char *   text, const Node * node) :
	Exception("NodeException: " + string(text) + (node != NULL ? "(" + node -> getName() + ")" : "")),
	node(node) {};
		
NodeException::NodeException(const string & text, const Node * node) :
	Exception("NodeException: " + text + (node != NULL ? "(" + node -> getName() + ")" : "")),
	node(node) {};
		
NodeException::~NodeException() throw() {};
	
const Node * NodeException::getNode() const { return node; }

/******************************************************************************/

UnrootedTreeException::UnrootedTreeException(const char *   text, const Tree * tree) :
	Exception("UnrootedTreeException: " + string(text) + (tree != NULL ? "(" + tree -> getName() + ")" : "")),
	tree(tree) {};
		
UnrootedTreeException::UnrootedTreeException(const string & text, const Tree * tree) :
	Exception("UnrootedTreeException: " + text + (tree != NULL ? "(" + tree -> getName() + ")" : "")),
	tree(tree) {};
		
UnrootedTreeException::~UnrootedTreeException() throw() {};
	
const Tree * UnrootedTreeException::getTree() const { return tree; }

/******************************************************************************/
