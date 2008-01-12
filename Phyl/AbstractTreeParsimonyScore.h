//
// File: AbstractTreeParsimonyScore.h
// Created by: Julien Dutheil
// Created on: Thu Jul 28 17:25 2005
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

#ifndef _ABSTRACTTREEPARSIMONYSCORE_H_
#define _ABSTRACTTREEPARSIMONYSCORE_H_

#include "TreeParsimonyScore.h"
#include "Node.h"

// From SeqLib:
#include <Seq/SiteContainer.h>

namespace bpp
{

/**
 * @brief Partial implementation of the TreeParsimonyScore interface.
 */
class AbstractTreeParsimonyScore :
	public virtual TreeParsimonyScore
{
	protected:
		TreeTemplate<Node> * _tree;
		const SiteContainer * _data;
		const Alphabet * _alphabet;
		unsigned int _nbStates;
		
	public:
		AbstractTreeParsimonyScore(
			const Tree & tree,
			const SiteContainer & data,
			bool verbose)
			throw (Exception);

    AbstractTreeParsimonyScore(const AbstractTreeParsimonyScore & tp):
      _tree(NULL), _data(NULL), _alphabet(tp._alphabet), _nbStates(tp._nbStates)
    {
      _tree     = tp._tree->clone();
      _data     = dynamic_cast<SiteContainer *>(tp._data->clone());
    }
    
    AbstractTreeParsimonyScore & operator=(const AbstractTreeParsimonyScore & tp)
    {
      _tree     = tp._tree->clone();
      _data     = dynamic_cast<SiteContainer *>(tp._data->clone());
      _alphabet = tp._alphabet;
      _nbStates = tp._nbStates;
      return *this;
    }

		virtual ~AbstractTreeParsimonyScore()
    {
      delete _tree;
      delete _data;
    }

	public:
		virtual const Tree * getTree() const { return _tree; }
		virtual vector<unsigned int> getScoreForEachSite() const;
};

} //end of namespace bpp.

#endif //_ABSTRACTTREEPARSIMONYSCORE_H_

