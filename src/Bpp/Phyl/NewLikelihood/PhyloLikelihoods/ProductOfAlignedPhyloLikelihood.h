//
// File: ProductOfAlignedPhyloLikelihood.h
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

#ifndef _PRODUCT_OFALIGNEDPHYLOLIKELIHOOD_H_
#define _PRODUCT_OFALIGNEDPHYLOLIKELIHOOD_H_

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignedValuesContainer.h>

#include "SetOfAlignedPhyloLikelihood.h"

namespace bpp
{

    /**
     * @brief The ProductOfAlignedPhyloLikelihood class, for phylogenetic
     * likelihood on several independent data.
     *
     */
    
    class ProductOfAlignedPhyloLikelihood:
    public SetOfAlignedPhyloLikelihood
    {
    public:
      ProductOfAlignedPhyloLikelihood(PhyloLikelihoodContainer* pC);

      ProductOfAlignedPhyloLikelihood(PhyloLikelihoodContainer* pC, const std::vector<size_t>& nPhylo);
      
      ~ProductOfAlignedPhyloLikelihood() {};
      
      ProductOfAlignedPhyloLikelihood* clone() const
      {
        return new ProductOfAlignedPhyloLikelihood(*this);
      }

      ProductOfAlignedPhyloLikelihood(const ProductOfAlignedPhyloLikelihood& sd);
        
      ProductOfAlignedPhyloLikelihood& operator=(const ProductOfAlignedPhyloLikelihood& sd);

    protected:
      
      void computeDLogLikelihood_(const std::string& variable) const
      {
        SetOfAbstractPhyloLikelihood::computeDLogLikelihood_(variable);
      }

      void computeD2LogLikelihood_(const std::string& variable) const
      {
        SetOfAbstractPhyloLikelihood::computeD2LogLikelihood_(variable);
      }
      
    public:

      /**
       *
       * @name Inherited from PhyloLikelihood
       *
       * @{
       */
      
      /**
       * @name The site likelihood functions.
       *
       */

      double getLogLikelihood() const;

      double getDLogLikelihood(const std::string& variable) const;

      double getD2LogLikelihood(const std::string& variable) const;
      

       /**
       * @brief Get the logarithm of the likelihood for any site.
       *
       * @return The logarithm of the likelihood of the data at this site.
       */

      virtual double getLikelihoodForASite(size_t site) const
      {
        updateLikelihood();
        computeLikelihood();
        
        double x=1;

        const std::vector<size_t>& nPhylo=getNumbersOfPhyloLikelihoods();
        
        for (size_t i=0; i<nPhylo.size(); i++)
          x *= getPhyloLikelihood(nPhylo[i])->getLikelihoodForASite(site);
        
        return x;
      }

      virtual double getLogLikelihoodForASite(size_t site) const
      {
        updateLikelihood();
        computeLikelihood();
        
        double x=0;

        const std::vector<size_t>& nPhylo=getNumbersOfPhyloLikelihoods();
        
        for (size_t i=0; i<nPhylo.size(); i++)
          x += getPhyloLikelihood(nPhylo[i])->getLogLikelihoodForASite(site);
        
        return x;
      }

      /**
       * @brief Get the derivates of the LogLikelihood at a Site
       *
       */

      double getDLogLikelihoodForASite(const std::string& variable, size_t site) const;
      
      double getD2LogLikelihoodForASite(const std::string& variable, size_t site) const;
      
      /** @} */

    };

} //end of namespace bpp.

#endif  //_PRODUCT_OFALIGNEDPHYLOLIKELIHOOD_H_


