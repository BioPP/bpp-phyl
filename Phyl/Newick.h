//
// File: Newick.h
// Created by: jdutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thu Oct 23 15:35:03 2003
//

#ifndef _NEWICK_H_
#define _NEWICK_H_

#include "IOTree.h"

/**
 * @brief The so-called 'newick' parenthetic format.
 *
 * Branch lengths and bootstraps are supported:
 *
 * ex:
 * <code>
 * ((A:0.1, B:0.15)90:0.2, C:0.27);
 * </code>
 *
 * Code example:
 * @code
 * #include <Phyl/Newick.h>
 * #include <Phyl/Tree.h>
 * 
 * Newick * newickReader = new Newick(false); //No comment allowed!
 * try {
 * 	Tree * tree = newickReader -> read("MyTestTree.dnd"); // Tree in file MyTestTree.dnd
 * 	cout << "Tree has " << tree -> getNumberOfLeaves() << " leaves." << endl;
 * } catch (Exception e) {
 *	cout << "Error when reading tree." << endl;
 * }
 * delete tree;
 * delete newickReader;
 * @endcode
 */
class Newick: public ITree, public OTree
{
	protected:
		bool _allowComments;
	
	public:
		/**
		 * @brief Build a new Newick reader/writer.
		 *
		 * Some newick format allow comments between hooks ('[' ']').
		 * 
		 * @param allowComments Tell if comments between [] are allowed in file.
		 */
		Newick(bool allowComments = false);
		virtual ~Newick();
	
	public:

		/**
		 * @name The IOTree interface
		 *
		 * @{
		 */
		const string getFormatName() const;
		const string getFormatDescription() const;
		/* @} */

		/**
		 * @name The ITree interface
		 *
		 * @{
		 */
		Tree * read(const string & path) const throw (Exception);
		/** @} */

		/**
		 * @name The OTree interface
		 *
		 * @{
		 */
		void write(const Tree & tree, const string & path, bool overwrite = true) const throw (Exception);
		/** @} */
};


#endif	//_NEWICK_H_
