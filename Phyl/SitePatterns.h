//
// File: SitePatterns.h
// Created by: Julien Dutheil
// Created on: Tue Nov 29 15:37 2005
//  from file PatternTools.h
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
 
#ifndef _SITEPATTERNS_H_
#define _SITEPATTERNS_H_

#include "Tree.h"

// From Utils:
#include <Utils/Clonable.h>

// From NumCalc:
#include <NumCalc/VectorTools.h>

//From SeqLib:
#include <Seq/SiteContainer.h>
#include <Seq/Site.h>

// From the STL:
#include <map>

using namespace std;

namespace bpp
{

/**
 * @brief Data structure for site patterns.
 * 
 * 'names' are the sequence names
 * 'sites' points toward a unique site
 * 'weights' is the number of sites identical to this sites
 * 'indices' are the positions in the original container
 */
class SitePatterns:
  public virtual Clonable
{
  private:
    /**
     * @brief Class used for site pattern sorting.
     */
    class SortableSite
    {
	    public:
		    SortableSite() {}
		    virtual ~SortableSite() {}
		
	    public:
		    string siteS; 
		    const Site * siteP;
		    unsigned int originalPosition;
    };

    /**
     * @brief Class used for site pattern sorting.
     */
    struct SSComparator : binary_function<SortableSite, SortableSite, bool>
    {
	    bool operator()(SortableSite * ss1, SortableSite * ss2) const { return ss1->siteS < ss2->siteS; }
    };

  protected: 
  	vector<string> _names;
	  vector<const Site *> _sites;
	  vector<unsigned int> _weights;
	  vector<unsigned int> _indices;
    const SiteContainer * _sequences;
    const Alphabet * _alpha;
    bool _own;

  public:
   /**
     * @brief Build a new SitePattern object.
     *
     * Look for patterns (unique sites) within a site container.
     *
     * @param sequences The container to look in.
     * @param own       Tel is the class own the sequence container.
     * If yes, the sequences wll be deleted together with this instance.
     */
    SitePatterns(const SiteContainer * sequences, bool own = false);

    virtual ~SitePatterns()
    {
      if(_own) delete _sequences;
    }

    SitePatterns(const SitePatterns & patterns)
    {
  	  _names     = patterns._names;
	    _sites     = patterns._sites;
	    _weights   = patterns._weights;
	    _indices   = patterns._indices;
      if(!patterns._own) _sequences = patterns._sequences;
      else               _sequences = dynamic_cast<SiteContainer *>(patterns._sequences->clone());
      _alpha     = patterns._alpha;
      _own       = patterns._own;
    }
    SitePatterns& operator=(const SitePatterns & patterns)
    {
  	  _names     = patterns._names;
	    _sites     = patterns._sites;
	    _weights   = patterns._weights;
	    _indices   = patterns._indices;
      if(!patterns._own) _sequences = patterns._sequences;
      else               _sequences = dynamic_cast<SiteContainer *>(patterns._sequences->clone());
      _alpha     = patterns._alpha;
      _own       = patterns._own;
      return *this;
    }

#ifdef NO_VIRTUAL_COV
    Clonable*
#else
    SitePatterns *
#endif
    clone() const { return new SitePatterns(*this); }

  public:
    /**
     * @return The number of times each unique site was found.
     */
		const vector<unsigned int> getWeights() const { return _weights; }
    /**
     * @return The position of each unique site.
     */
		const vector<unsigned int> getIndices() const { return _indices; }

    /**
     * @return A new container with each unique site.
     */
		SiteContainer * getSites() const;
    
};

} //end of namespace bpp.

#endif // _SITEPATTERNS_H_

