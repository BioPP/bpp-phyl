//
// File: LikelihoodCalculation.h
// Authors:
//   FranÃ§ois Gindraud, Laurent GuÃ©guen (2018)
// Created: jeudi 28 fÃ©vrier 2019, Ã  07h 22
//

/*
  Copyright or Â© or Copr. Bio++ Development Team, (November 16, 2004)
  
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

#ifndef BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_LIKELIHOODCALCULATION_H
#define BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_LIKELIHOODCALCULATION_H

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Eigen/Core>

#include "Bpp/Phyl/NewLikelihood/DataFlow/DataFlowNumeric.h"
#include "Definitions.h"

namespace bpp
{
/** Base classes for Likelihood calculations. Also used to store shared_ptr.
 *
 */

class LikelihoodCalculation :
  public AbstractParameterAliasable
{
private:
  Context& context_;

protected:
  /******************************************/
  /** Likelihoods  **/

  mutable ValueRef<DataLik> likelihood_;

public:
  LikelihoodCalculation(Context& context) :
    AbstractParameterAliasable(""),
    context_(context)
  {}

  LikelihoodCalculation(Context& context,
                        ParameterList& paramList) :
    AbstractParameterAliasable(""),
    context_(context)
  {
    shareParameters_(paramList);
  }


  /*
   * @brief Copy the likelihood calculation IN THE SAME CONTEXT.
   *
   */

  LikelihoodCalculation(const LikelihoodCalculation& lik) :
    AbstractParameterAliasable(lik),
    context_(lik.context_),
    likelihood_()
  {}

  LikelihoodCalculation* clone() const
  {
    throw bpp::Exception("LikelihoodCalculation clone should not happen.");
  }

  ValueRef<DataLik> getLikelihoodNode()
  {
    if (!likelihood_)
      makeLikelihoods();

    return likelihood_;
  }

  virtual bool isInitialized() const
  {
    return true;
  }

  const Context& getContext() const
  {
    return context_;
  }

  /*
   * @brief Return likelihood_ without any computation, used to build likelihood_.
   *
   */
  void setLikelihoodNode(ValueRef<DataLik> ll)
  {
    likelihood_ = ll;
  }

  /*
   * @brief fix Factor such that valRef value becomes normal.
   *
   */

protected:
  /* @brief Build the likelihood DF  */

  virtual void makeLikelihoods()
  {
    throw Exception("LikelihoodCalculation:: makeLikelihoods should not be called. Probably likelihood_ is null.");
  }

  Context& getContext_()
  {
    return context_;
  }


  ValueRef<DataLik> getLikelihoodNode_()
  {
    return likelihood_;
  }
};

using SiteLikelihoods = Value<RowLik>;
using SiteLikelihoodsRef = ValueRef<RowLik>;

class AlignedLikelihoodCalculation :
  public LikelihoodCalculation
{
public:
  AlignedLikelihoodCalculation(Context& context) :
    LikelihoodCalculation(context) {}

  AlignedLikelihoodCalculation(const AlignedLikelihoodCalculation& lik) :
    LikelihoodCalculation(lik)
  {}


  AlignedLikelihoodCalculation* clone() const
  {
    throw bpp::Exception("AlignedLikelihoodCalculation clone should not happen.");
  }

protected:
  /******************************************/
  /** Site Likelihoods  **/

  /* For Data used for output (ie non shrunk if ever)
   *
   */

  mutable SiteLikelihoodsRef siteLikelihoods_;

  /* For Data used for computation (ie shrunked data if ever)
   *
   */

  mutable SiteLikelihoodsRef patternedSiteLikelihoods_;

public:
  /*
   * @brief Return the ref to the Sitelikelihoods_ vector on data
   * (shrunked or not).
   *
   * @brief shrunk: bool true if vector on shrunked data (default: false)
   *
   */
  SiteLikelihoodsRef getSiteLikelihoods(bool shrunk = false)
  {
    if (!(siteLikelihoods_ || patternedSiteLikelihoods_))
      makeLikelihoods();

    if (shrunk && patternedSiteLikelihoods_)
      return patternedSiteLikelihoods_;
    else
      return siteLikelihoods_;
  }

  /*
   * @brief Return the ref to the Sitelikelihoods_ vector on data
   * (shrunked or not).
   *
   * @brief ll: Site Likelihoods
   * @brief shrunk: given Likelihoods are on shrunked data (default false)
   *
   */
  void setSiteLikelihoods(SiteLikelihoodsRef ll,
                          bool shrunk = false)
  {
    if (shrunk)
      patternedSiteLikelihoods_ = ll;
    else
      siteLikelihoods_ = ll;
  }

  /*
   * @brief Return Likelihood on a site
   * @param pos : site position
   * @param shrunk : on shrunked data (default false)
   */
  DataLik getLikelihoodForASite(size_t pos, bool shrunk = false)
  {
    if (!(siteLikelihoods_ || patternedSiteLikelihoods_))
      makeLikelihoods();

    if (shrunk && patternedSiteLikelihoods_)
      return patternedSiteLikelihoods_->getTargetValue()(Eigen::Index(pos));
    else
      return siteLikelihoods_->getTargetValue()(Eigen::Index(pos));
  }

  DataLik getLogLikelihoodForASite(size_t pos, bool shrunk = false)
  {
    using namespace std;
    return log(getLikelihoodForASite(pos, shrunk));
  }

  /**
   * @brief Get the likelihood for each site.
   *
   *@return A vector with all site likelihoods.
   *
   */
  VDataLik getLikelihoodPerSite()
  {
    auto vLik = getSiteLikelihoods(false)->getTargetValue();
    VDataLik v;
    copyEigenToBpp(vLik, v);
    return v;
  }

protected:
  virtual void makeLikelihoods()
  {
    throw Exception("AlignedLikelihoodCalculation:: makeLikelihoods should not be called. Probably both siteLikelihoods_ and patternedSiteLikelihoods_ are null.");
  }
};
} // namespace bpp
#endif // BPP_PHYL_NEWLIKELIHOOD_DATAFLOW_LIKELIHOODCALCULATION_H
