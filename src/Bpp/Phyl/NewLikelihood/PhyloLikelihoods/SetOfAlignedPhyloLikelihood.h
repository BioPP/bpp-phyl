//
// File: SetOfAlignedPhyloLikelihood.h
// Created by: Laurent Guéguen
// Created on: mercredi 7 octobre 2015, à 13h 30
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

#ifndef _SET_OF_ALIGNED_PHYLOLIKELIHOOD_H_
#define _SET_OF_ALIGNED_PHYLOLIKELIHOOD_H_

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>

#include "AlignedPhyloLikelihood.h"
#include "SetOfAbstractPhyloLikelihood.h"

namespace bpp
{
    /**
     * @brief The SetOfAlignedPhyloLikelihood abstract class.
     *
     * This class defines the common methods needed to compute a
     * likelihood from aligned phylogenies.
     */
    
    class SetOfAlignedPhyloLikelihood:
    public SetOfAbstractPhyloLikelihood,
    virtual public AbstractAlignedPhyloLikelihood
    {
    public:
      SetOfAlignedPhyloLikelihood(PhyloLikelihoodContainer* pC, const std::string& prefix = "");
      
      SetOfAlignedPhyloLikelihood(PhyloLikelihoodContainer* pC, const std::vector<size_t>& nPhylo, const std::string& prefix = "");
      
      SetOfAlignedPhyloLikelihood(const SetOfAlignedPhyloLikelihood& soap) :
      AbstractPhyloLikelihood(soap),
      AbstractAlignedPhyloLikelihood(soap),
      SetOfAbstractPhyloLikelihood(soap)
      {
      }

      SetOfAlignedPhyloLikelihood& operator=(const SetOfAlignedPhyloLikelihood& soap)
      {
        SetOfAbstractPhyloLikelihood::operator=(soap);
        AbstractAlignedPhyloLikelihood::operator=(soap);

        return *this;
      }

      virtual SetOfAlignedPhyloLikelihood* clone() const = 0;

      virtual ~SetOfAlignedPhyloLikelihood()
      {
      }

    public:

      /**
       *
       * @brief adds a PhyloLikelihood already stored in the m ap, iff
       * it is an AlignedPhyloLikelihood of the same size 
       *
       * @return if the PhyloLikelihood has been added.
       */
      
      virtual bool addPhyloLikelihood(size_t nPhyl);

      const AbstractAlignedPhyloLikelihood* getAbstractPhyloLikelihood(size_t nPhyl) const
      {
        return dynamic_cast<const AbstractAlignedPhyloLikelihood*>((*pPhyloCont_)[nPhyl]);
      }
      
      
      AbstractAlignedPhyloLikelihood* getAbstractPhyloLikelihood(size_t nPhyl)
      {
        return dynamic_cast<AbstractAlignedPhyloLikelihood*>((*pPhyloCont_)[nPhyl]);
      }

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
       * @brief To be defined in inheriting classes.
       *
       * @brief Get the likelihood for a site.
       *
       * @param site The site index to analyse.
       * @return The likelihood for site <i>site</i>.
       */

      virtual double getLikelihoodForASite(size_t site) const = 0;

      virtual double getLogLikelihoodForASite(size_t site) const = 0;

      /**
       * @brief Get the likelihood for a site for an aligned
       * phyloLikelihood 
       *
       * @param site The site index to analyse.
       * @param p the phyloLikelihood index.
       * @return The likelihood for site <i>site</i>.
       */

  
      double getLikelihoodForASiteForAPhyloLikelihood(size_t site, size_t nPhyl) const
      {
        return getAbstractPhyloLikelihood(nPhyl)->getLikelihoodForASite(site);
      }

      /**
       * @brief Get the first derivative of the likelihood for a site
       * for a process.
       *
       * @param site The site index to analyse.
       * @param p the phyloLikelihood index.
       * @return The likelihood for site <i>site</i>.
       */

      double getDLogLikelihoodForASiteForAPhyloLikelihood(size_t site, size_t nPhyl) const
      {
        return getAbstractPhyloLikelihood(nPhyl)->getDLogLikelihoodForASite(site);
      }
      
      /**
       * @brief Get the second derivative of the likelihood for a site for
       * a process.
       *
       * @param site The site index to analyse.
       * @param p the phyloLikelihood index.
       * @return The likelihood for site <i>site</i>.
       */

      double getD2LogLikelihoodForASiteForAPhyloLikelihood(size_t site, size_t nPhyl) const
      {
        return getAbstractPhyloLikelihood(nPhyl)->getD2LogLikelihoodForASite(site);
      }

    };
  
  
} //end of namespace bpp.

#endif  //_SET_OF_ALIGNED_PHYLOLIKELIHOOD_H_

