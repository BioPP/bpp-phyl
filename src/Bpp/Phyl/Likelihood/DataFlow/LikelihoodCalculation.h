// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_DATAFLOW_LIKELIHOODCALCULATION_H
#define BPP_PHYL_LIKELIHOOD_DATAFLOW_LIKELIHOODCALCULATION_H

#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Eigen/Core>

#include "Bpp/Phyl/Likelihood/DataFlow/DataFlowNumeric.h"
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
    context_(context),
    likelihood_()
  {}

  LikelihoodCalculation(Context& context,
      ParameterList& paramList) :
    AbstractParameterAliasable(""),
    context_(context),
    likelihood_()
  {
    shareParameters_(paramList);
  }


  /**
   * @brief Copy the likelihood calculation IN THE SAME CONTEXT.
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
    return likelihood_ != 0;
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

  virtual void cleanAllLikelihoods()
  {
    if (likelihood_)
      getContext_().erase(likelihood_);
    likelihood_.reset();
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
    LikelihoodCalculation(context),
    siteLikelihoods_(),
    patternedSiteLikelihoods_()
  {}

  AlignedLikelihoodCalculation(const AlignedLikelihoodCalculation& lik) :
    LikelihoodCalculation(lik),
    siteLikelihoods_(),
    patternedSiteLikelihoods_()
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
      return patternedSiteLikelihoods_->targetValue()(Eigen::Index(pos));
    else
      return siteLikelihoods_->targetValue()(Eigen::Index(pos));
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
    auto vLik = getSiteLikelihoods(false)->targetValue();
    VDataLik v;
    copyEigenToBpp(vLik, v);
    return v;
  }

  /*
   * @brief Clean all the existing likelihoods (ie remove them from
   * Context), and all unique dependencies.
   *
   *
   */
  void cleanAllLikelihoods()
  {
    LikelihoodCalculation::cleanAllLikelihoods();

    if (siteLikelihoods_)
    {
      getContext_().erase(siteLikelihoods_);
      siteLikelihoods_.reset();
    }
    if (patternedSiteLikelihoods_)
    {
      getContext_().erase(patternedSiteLikelihoods_);
      patternedSiteLikelihoods_.reset();
    }
  }

protected:
  virtual void makeLikelihoods()
  {
    throw Exception("AlignedLikelihoodCalculation:: makeLikelihoods should not be called. Probably both siteLikelihoods_ and patternedSiteLikelihoods_ are null.");
  }
};
} // namespace bpp
#endif // BPP_PHYL_LIKELIHOOD_DATAFLOW_LIKELIHOODCALCULATION_H
