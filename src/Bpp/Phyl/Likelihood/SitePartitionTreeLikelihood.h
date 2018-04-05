//
// File: SitePartitionTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Mon Dec 7 11:02 2009
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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


#ifndef _SITEPARTITIONTREELIKELIHOOD_H_
#define _SITEPARTITIONTREELIKELIHOOD_H_

#include "TreeLikelihood.h"

namespace bpp
{

/**
 * @brief Specialization of the TreeLikelihood interface for partition models, homogeneous case.
 *
 * These models allow the distinct sites of an alignment to have a different model.
 * The substitution model is however assumed to be the same along the tree.
 * Such models are hence homogeneous in time.
 */
class SitePartitionHomogeneousTreeLikelihood :
  public virtual TreeLikelihood
{
  public:

    SitePartitionHomogeneousTreeLikelihood* clone() const = 0;

  public:
    const TransitionModel* getModel(int nodeId, size_t siteIndex) const
    {
      return getModelForSite(siteIndex);
    }

    TransitionModel* getModel(int nodeId, size_t siteIndex)
    {
      return getModelForSite(siteIndex);
    }

    /**
     * @brief Get the substitution model associated to a given node.
     *
     * @param siteIndex The position in the alignment.
     * @return A pointer toward the corresponding model.
     */
    virtual const TransitionModel* getModelForSite(size_t siteIndex) const = 0;

    /**
     * @brief Get the substitution model associated to a given node.
     *
     * @param siteIndex The position in the alignment.
     * @return A pointer toward the corresponding model.
     */
    virtual TransitionModel* getModelForSite(size_t siteIndex) = 0;

};

} //end of namespace bpp.

#endif  //_SITEPARTITIONTREELIKELIHOOD_H_

