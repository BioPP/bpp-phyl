//
// File: MixturePhyloLikelihood.h
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

#ifndef _MIXTUREPHYLOLIKELIHOOD_H_
#define _MIXTUREPHYLOLIKELIHOOD_H_


#include "MultiPhyloLikelihood.h"

#include <Bpp/Numeric/Prob/Simplex.h>

//From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>


namespace bpp
{
  namespace newlik
  {

    /**
     * @brief Likelihood Collection made of a mixture of
     * TreeLikelihoods. The resulting likelihood is the mean value of
     * the TreeLikelihoods, ponderated with parametrized probabilities
     * (through a Simplex).
     *
     */
    
    class MixtureLikelihoodCollection:
      public LikelihoodCollection
    {
    private:

      Simplex simplex_;
      
    public:
      MixtureLikelihoodCollection(SubstitutionProcessCollection* processColl,
                                  bool verbose = true,
                                  bool patterns = true);

      MixtureLikelihoodCollection(const SiteContainer& data,
                                  SubstitutionProcessCollection* processColl,
                                  bool verbose = true,
                                  bool patterns = true);
      
      MixtureLikelihoodCollection(const MixtureLikelihoodCollection& mlc) : LikelihoodCollection(mlc), simplex_(mlc.simplex_) {}

      MixtureLikelihoodCollection& operator=(const MixtureLikelihoodCollection& mlc)
      {
        LikelihoodCollection::operator=(mlc);
        simplex_=mlc.simplex_;
        return *this;
      }

      virtual ~MixtureLikelihoodCollection() {}

      MixtureLikelihoodCollection* clone() const { return new MixtureLikelihoodCollection(*this);}

    public:

      void fireParameterChanged(const ParameterList & parameters);

      ParameterList getDerivableParameters() const { return getBranchLengthsParameters(); }

      ParameterList getNonDerivableParameters() const;

      /**
       * @name The likelihood functions.
       *
       * @{
       */

      /**
       * @brief Get the likelihood for a site.
       *
       * @param site The site index to analyse.
       * @return The likelihood for site <i>site</i>.
       */
      double getLikelihoodForASite(size_t site) const;
 
      double getDLikelihoodForASite(size_t site) const;

      double getD2LikelihoodForASite(size_t site) const;

      
      /*
       * @}
       */

    };

  } //end of namespace newlik.
} //end of namespace bpp.

#endif  //_MIXTUREPHYLOLIKELIHOOD_H_

