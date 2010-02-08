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

#ifndef ABSTRACTTREEPARSIMONYSCORE_H__
#define ABSTRACTTREEPARSIMONYSCORE_H__

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
	private:
		TreeTemplate<Node>* tree_;
		const SiteContainer* data_;
		const Alphabet* alphabet_;
		unsigned int nbStates_;
		
	public:
		AbstractTreeParsimonyScore(
			const Tree& tree,
			const SiteContainer & data,
			bool verbose)
			throw (Exception);

    AbstractTreeParsimonyScore(const AbstractTreeParsimonyScore& tp):
      tree_(0), data_(0), alphabet_(tp.alphabet_), nbStates_(tp.nbStates_)
    {
      tree_     = tp.tree_->clone();
      data_     = dynamic_cast<SiteContainer *>(tp.data_->clone());
    }
    
    AbstractTreeParsimonyScore& operator=(const AbstractTreeParsimonyScore & tp)
    {
      tree_     = dynamic_cast<TreeTemplate<Node>*>(tp.tree_->clone());
      data_     = dynamic_cast<SiteContainer *>(tp.data_->clone());
      alphabet_ = tp.alphabet_;
      nbStates_ = tp.nbStates_;
      return *this;
    }

		virtual ~AbstractTreeParsimonyScore()
    {
      delete tree_;
      if (data_) delete data_;
    }

	public:
		virtual const Tree& getTree() const { return *tree_; }
		virtual std::vector<unsigned int> getScoreForEachSite() const;

  protected:
    const TreeTemplate<Node>* getTreeP_() const { return tree_; }
    TreeTemplate<Node>* getTreeP_() { return tree_; }
};

} //end of namespace bpp.

#endif //ABSTRACTTREEPARSIMONYSCORE_H__

