//
// File: AlignedPhyloLikelihood_DF.h
// Authors: Laurent Guéguen (2019)
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef ALIGNED_PHYLOLIKELIHOOD_DF_H
#define ALIGNED_PHYLOLIKELIHOOD_DF_H

#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/AlignedPhyloLikelihood.h>

#include "PhyloLikelihood_DF.h"

namespace bpp {

  namespace dataflow
  {
    
    /* Wraps a dataflow graph as a function: resultNode = f(variableNodes).
     *
     */
  
    class AlignedPhyloLikelihood_DF :
      virtual public PhyloLikelihood_DF,
      virtual public AlignedPhyloLikelihood
    {
    public:
      AlignedPhyloLikelihood_DF (dataflow::Context & context,
                                 std::shared_ptr<LikelihoodCalculationSingleProcess> likCal,
                                 const ParameterList & variableNodes)
        : PhyloLikelihood_DF(context, likCal, variableNodes)
      {}

      /*
       * @brief: the parameters are those of the LikelihoodCalculationSingleProcess
       */
      
      AlignedPhyloLikelihood_DF (dataflow::Context & context, std::shared_ptr<LikelihoodCalculationSingleProcess> likCal)
        : PhyloLikelihood_DF(context, likCal)
      {}

      // Legacy boilerplate
      AlignedPhyloLikelihood_DF * clone () const override { return new AlignedPhyloLikelihood_DF (*this); }

      size_t getNumberOfSites() const {
        return getLikelihoodCalculation()->getNumberOfSites();
      }

      size_t getNumberOfDistinctSites() const {
        return getLikelihoodCalculation()->getNumberOfDistinctSites();
      }
      
      /**
       * @brief Get the likelihood for a site.
       *
       * @param site The site index to analyse.
       * @return The likelihood for site <i>site</i>.
       */

      double getLikelihoodForASite(size_t site) const
      {
        return getLikelihoodCalculation()->getSiteLikelihoods()->getTargetValue()[site];
      }
      
      /**
       * @brief Get the log likelihood for a site, and its derivatives.
       *
       * @param site The site index to analyse.
       * @return The (D)log likelihood for site <i>site</i>.
       */

      double getLogLikelihoodForASite(size_t site) const
      {
        return std::log(getLikelihoodForASite(site));
      }
    
      double getDLogLikelihoodForASite(const std::string& variable, size_t site) const
      {
        return getFirstOrderDerivativeVector(variable)->getTargetValue()[site];
      }
      
      double getD2LogLikelihoodForASite(const std::string& variable, size_t site) const
      {
        return getSecondOrderDerivativeVector(variable)->getTargetValue()[site];
      }
    
      /**
       * @brief Get the likelihood for each site.
       *
       * @return A vector with all likelihoods for each site.
       */

      Vdouble getLikelihoodPerSite() const
      {
        auto vLik=getLikelihoodCalculation()->getSiteLikelihoods()->getTargetValue();
        
        Vdouble v(vLik.size());
        
        Eigen::VectorXd::Map(&v[0], v.size()) = vLik;
        return v;
      }
      
/** @} */

    };

  }
  
} // namespace bpp

#endif // DATA_FLOW_FUNCTION_H
