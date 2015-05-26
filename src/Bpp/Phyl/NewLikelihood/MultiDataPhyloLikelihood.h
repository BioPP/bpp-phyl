//
// File: MultiDataPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: jeudi 11 juillet 2013, à 14h 05
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

#ifndef _MULTI_DATA_PHYLOLIKELIHOOD_H_
#define _MULTI_DATA_PHYLOLIKELIHOOD_H_

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>

#include "PhyloLikelihood.h"
#include "SingleDataPhyloLikelihood.h"

namespace bpp
{
    /**
     * @brief The MultiDataPhyloLikelihood interface, for phylogenetic likelihood.
     *
     * This interface defines the common methods needed to compute a likelihood
     * from a sequence alignement, usually involving one or more phylogenetic trees.
     */
    
    class MultiDataPhyloLikelihood:
      public virtual PhyloLikelihood
    {
    public:
      MultiDataPhyloLikelihood() {}
      virtual ~MultiDataPhyloLikelihood() {}

      virtual MultiDataPhyloLikelihood* clone() const = 0;
      

    public:

      /**
       *
       * @name The data functions
       *
       * @{
       */
      
      /**
       * @brief Set the dataset for which the likelihood must be evaluated.
       *
       * @param i the number of the SinglePhylolikelihood that use the data
       * @param sites The data set to use.
       */

      virtual void setData(const SiteContainer& sites, size_t i) = 0;
      
    
      /**
       * @brief Get the dataset for which the likelihood must be evaluated.
       *
       *
       * @return A pointer toward the site container where the sequences are stored.
       */
      
      virtual const SiteContainer* getData(size_t i) const = 0;

      /**
       * @}
       */

      /**
       *
       * @name The single data Phylolikelihood storage.
       *
       * @{
       */

      virtual void addSingleDataPhylolikelihood(size_t, SingleDataPhyloLikelihood*) = 0;

      virtual std::vector<size_t> getNumbersOfSingleDataPhyloLikelihoods() const =0;

      virtual const SingleDataPhyloLikelihood* getSingleDataPhylolikelihood(size_t) const = 0;

      virtual SingleDataPhyloLikelihood* getSingleDataPhylolikelihood(size_t) = 0;
      
      /**
       *
       * @}
       *
       */

    };

} //end of namespace bpp.

#endif  //_MULTI_DATA_PHYLOLIKELIHOOD_H_

