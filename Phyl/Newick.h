//
// File: Newick.h
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Thu Oct 23 15:35:03 2003
//

/*
Copyright ou � ou Copr. CNRS, (16 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant � fournir des classes
pour l'analyse de donn�es phylog�n�tiques.

Ce logiciel est r�gi par la licence CeCILL soumise au droit fran�ais et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffus�e par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilit� au code source et des droits de copie,
de modification et de redistribution accord�s par cette licence, il n'est
offert aux utilisateurs qu'une garantie limit�e.  Pour les m�mes raisons,
seule une responsabilit� restreinte p�se sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les conc�dants successifs.

A cet �gard  l'attention de l'utilisateur est attir�e sur les risques
associ�s au chargement,  � l'utilisation,  � la modification et/ou au
d�veloppement et � la reproduction du logiciel par l'utilisateur �tant 
donn� sa sp�cificit� de logiciel libre, qui peut le rendre complexe � 
manipuler et qui le r�serve donc � des d�veloppeurs et des professionnels
avertis poss�dant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invit�s � charger  et  tester  l'ad�quation  du
logiciel � leurs besoins dans des conditions permettant d'assurer la
s�curit� de leurs syst�mes et ou de leurs donn�es et, plus g�n�ralement, 
� l'utiliser et l'exploiter dans les m�mes conditions de s�curit�. 

Le fait que vous puissiez acc�der � cet en-t�te signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accept� les
termes.
*/

/*
Copyright or � or Copr. CNRS, (November 16, 2004)

Julien.Dutheil@univ-montp2.fr

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
 * 	Tree<Node> * tree = newickReader -> read("MyTestTree.dnd"); // Tree in file MyTestTree.dnd
 * 	cout << "Tree has " << tree -> getNumberOfLeaves() << " leaves." << endl;
 * } catch (Exception e) {
 *	cout << "Error when reading tree." << endl;
 * }
 * delete tree;
 * delete newickReader;
 * @endcode
 */
class Newick: public virtual AbstractITree, public virtual AbstractOTree
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
		Tree<Node> * read(const string & path) const throw (Exception)
		{
			return AbstractITree::read(path);
		}
		Tree<Node> * read(istream & in) const throw (Exception);
		/** @} */

		/**
		 * @name The OTree interface
		 *
		 * @{
		 */
		void write(const Tree<Node> & tree, const string & path, bool overwrite = true) const throw (Exception)
		{
			AbstractOTree::write(tree, path, overwrite);
		}
		void write(const Tree<Node> & tree, ostream & out) const throw (Exception);
		/** @} */
};


#endif	//_NEWICK_H_
