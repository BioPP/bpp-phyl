//
// File: TreeExceptions.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Mon Nov  3 17:04:46 2003
//

#ifndef _TREEEXCEPTIONS_H_
#define _TREEEXCEPTIONS_H_

// From Utils:
#include <Utils/Exceptions.h>

class Node;
class Tree;

/******************************************************************************
 *                            Tree exceptions:                                *
 ******************************************************************************/

/**
 * @brief <p>Exception thrown when something is wrong with a particular node.</p>
 */
class NodeException : public Exception {

	protected:
		const Node * node;
			
	public:	// Class constructors and destructor:

		/**
		 * @brief Build a new NodeException.
		 * @param text A message to be passed to the exception hierarchy.
		 * @param node A const pointer toward the node that threw the exception.
		 */
		NodeException(const char * text, const Node * node = NULL);

		/**
		 * @brief Build a new NodeException.
		 * @param text A message to be passed to the exception hierarchy.
		 * @param node A const pointer toward the node that threw the exception.
		 */
		NodeException(const string & text, const Node * node = NULL);
	
		virtual ~NodeException() throw ();
	
	public:
		/**
		 * @brief <p>Get the node that threw the exception.</p>
		 * @return A pointer toward the faulty node.
		 */
		virtual const Node * getNode() const;
};

/******************************************************************************/

class UnrootedTreeException : public Exception {

	protected:
		const Tree * tree;
			
	public: // Class constructors and destructor:

		/**
		 * @brief Build a new UnrootedTreeException.
		 * 
		 * @param text A message to be passed to the exception hierarchy.
		 * @param tree A const pointer toward the tree that threw the exception.
		 */
		UnrootedTreeException(const char *   text, const Tree * tree = NULL);
	
		/**
		 * @brief Build a new UnrootedTreeException.
		 * 
		 * @param text A message to be passed to the exception hierarchy.
		 * @param tree A const pointer toward the tree that threw the exception.
		 */
		UnrootedTreeException(const string & text, const Tree * tree = NULL);
	
		virtual ~UnrootedTreeException() throw ();
	
	public:
		/**
		 * @brief <p>Get the tree that threw the exception.</p>
		 * @return The faulty tree/
		 */
		virtual const Tree * getTree() const;
};

/******************************************************************************/

#endif	//_TREEEXCEPTIONS_H_
