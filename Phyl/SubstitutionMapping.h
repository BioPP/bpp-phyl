//
// File: SubstitutionMapping.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 09:51 2005
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004, 2005, 2006)

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

#ifndef _SUBSTITUTIONMAPPING_H_
#define _SUBSTITUTIONMAPPING_H_

#include "Tree.h"
#include "TreeTemplate.h"

//From Utils:
#include <Utils/Clonable.h>

namespace bpp
{

/**
 * @brief General interface for storing mapping data.
 *
 * There are several kinds of mapping:
 * - Exact mapping, storing the positions of each substitution onto each branch,
 * - Probabilistic mapping, storing the number of substitutions onto each branch.
 *
 * Since only probabilistic substitution mapping is implemented for now, the basal 
 * interfac only contains one method.
 * More methods are expected to be added later.
 */
class SubstitutionMapping:
  public virtual Clonable
{

  public:
    SubstitutionMapping() {}
    virtual ~SubstitutionMapping() {}

#ifndef NO_VIRTUAL_COV
    SubstitutionMapping * clone() const = 0;
#endif

  public:
    
    /**
     * @return Get the phylogenetic tree associated to this mapping.
     */
    virtual const Tree * getTree() const = 0;
    /**
     * @return The number of sites mapped.
     */
    virtual unsigned int getNumberOfSites() const = 0;
    /**
     * @return The number of branches mapped.
     */
    virtual unsigned int getNumberOfBranches() const = 0;
    
    /**
     * @param index The site index.
     * @return The site position corresponding to the index.
     */
    virtual int getSitePosition(unsigned int index) const = 0;
};






/**
 * @brief Partial implementation of the substitution mapping interface.
 */
class AbstractSubstitutionMapping:
  public SubstitutionMapping
{
  protected:
    const TreeTemplate<Node> * _tree;
    vector<int> _sitesPostions;

  public:
    AbstractSubstitutionMapping() {}
    virtual ~AbstractSubstitutionMapping() {}

  public:

		virtual const	TreeTemplate<Node> * getTree() const { return _tree; }
 
    int getSitePosition(unsigned int index) const { return _sitesPostions[index]; }
};

} //end of namespace bpp.

#endif //_SUBSTITUTIONMAPPING_H_

