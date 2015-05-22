//
// File: MultiPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: jeudi 14 mai 2015, à 20h 58
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

#ifndef _MULTI_PHYLOLIKELIHOOD_H_
#define _MULTI_PHYLOLIKELIHOOD_H_

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>

#include "PhyloLikelihood.h"

namespace bpp
{
  namespace newlik
  {

    /**
     * @brief The MultiPhyloLikelihood interface, for phylogenetic likelihood.
     *
     * This interface defines the common methods needed to compute a likelihood
     * from several "simple"  phylolikelihoods.
     */
    
    class MultiPhyloLikelihood:
      public virtual PhyloLikelihood
    {
    public:
      MultiPhyloLikelihood() {}
      virtual ~MultiPhyloLikelihood() {}

      virtual MultiPhyloLikelihood* clone() const = 0;
      

    public:

      /**
       *
       * @name The single data Phylolikelihood storage.
       *
       * @{
       */

      virtual void addPhylolikelihood(size_t, PhyloLikelihood*) = 0;

      virtual std::vector<size_t> getNumbersOfPhyloLikelihoods() const =0;

      virtual const PhyloLikelihood* getPhylolikelihood(size_t) const = 0;

      virtual PhyloLikelihood* getPhylolikelihood(size_t) = 0;
      
      /**
       *
       * @}
       *
       */

    };

  } //end of namespace newlik.
} //end of namespace bpp.

#endif  //_MULTI_PHYLOLIKELIHOOD_H_

