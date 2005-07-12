//
// File: TreeExceptions.h
// Created by: Julien Dutheil
// Created on: Mon Nov  3 17:04:46 2003
//

#include "TreeExceptions.h"
#include "Node.h"
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

NodeNotFoundException::NodeNotFoundException(const char *   text, const string & id) :
	Exception("NodeNotFoundException: " + string(text) + "(" + id + ")"), _id(id) {};
		
NodeNotFoundException::NodeNotFoundException(const string & text, const string & id) :
	Exception("NodeNotFoundException: " + text + "(" + id + ")"), _id(id) {};
		
NodeNotFoundException::~NodeNotFoundException() throw() {};
	
string NodeNotFoundException::getId() const { return _id; }

/******************************************************************************/

UnrootedTreeException::UnrootedTreeException(const char *   text, const Tree * tree) :
			Exception("UnrootedTreeException: " + string(text) + (tree != NULL ? "(" + tree -> getName() + ")" : "")),
			tree(tree) {};

UnrootedTreeException::UnrootedTreeException(const string & text, const Tree * tree) :
			Exception("UnrootedTreeException: " + text + (tree != NULL ? "(" + tree -> getName() + ")" : "")),
			tree(tree) {};

/******************************************************************************/

