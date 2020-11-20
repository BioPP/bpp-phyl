//
// File: MultinomialFromTransitionModel.h
// Created by: Laurent Gueguen
// Created on: dimanche 26 janvier 2020, à 07h 52
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

#ifndef _MULTINOMIAL_FROM_TRANSITION_MODEL_H_
#define _MULTINOMIAL_FROM_TRANSITION_MODEL_H_

#include "AbstractWrappedModel.h"
#include <unsupported/Eigen/SpecialFunctions>

// namespace std
// {
// //  template<class U = Eigen::VectorXd>
//   template<>


// }

namespace bpp
{
/**
 * @brief From a model, compute the likelihood of counts given an
 * ancestral state.
 *
 * It has the same parameters as the SubModel.
 *
 * If counts are not integers, rounded values are used.
 *
 */

  class MultinomialFromTransitionModel :
    virtual public AbstractParameterAliasable,
    virtual public AbstractWrappedModel
  {
    typedef bool lessEigenType(Eigen::VectorXd const&, Eigen::VectorXd const&);

  private:
   /*
    * @brief The related model.
    *
    */
    
    std::unique_ptr<TransitionModel> subModel_;

    /*
     * The number of states
     *
     */
    
    size_t size_;

    /*
     * @brief Reference time to avoid recomuputation of transition
     * matrix when time has not changed. If <0, it means that
     * transition matrix should be recomputed (for ex if parameters
     * have changed).
     *
     */

    mutable double tref_;

    /*
     * @brief Transition Matrices owned by the submodel.
     *
     */

    /*
     * @brief These ones are for bookkeeping:
     */
    
    mutable const Matrix<double>* Pij_t, *dPij_dt, *d2Pij_dt2;

    /*
     * @brief Used return vectors
     *
     */

    mutable Eigen::VectorXd Pi_, dPi_, d2Pi_;

    /*
     * @brief map to store constant values of multinomial
     *
     */

    mutable std::map<Eigen::VectorXd, double, bool(*)(const Eigen::VectorXd&, const Eigen::VectorXd&)> mapFact_;
    
    static bool lessEigen(Eigen::VectorXd const& a, Eigen::VectorXd const& b)
    {
      assert(a.size()==b.size());
      for(auto i=0;i<a.size();++i)
      {
        if(a[i]<b[i]) return true;
        if(a[i]>b[i]) return false;
      }
      return false;
    };
  
  protected:
    BranchModel& getModel() 
    {
      return *subModel_.get();
    }

    TransitionModel& getTransitionModel()
    {
      return *subModel_.get();
    }

  public:
    MultinomialFromTransitionModel(const TransitionModel& originalModel) :
      AbstractParameterAliasable("MultinomialFrom."+originalModel.getNamespace()),
      AbstractWrappedModel(),
      subModel_(std::unique_ptr<TransitionModel>(originalModel.clone())),
      size_(originalModel.getNumberOfStates()),
      tref_(NumConstants::MINF()), Pij_t(0), dPij_dt(0), d2Pij_dt2(0), Pi_(size_), dPi_(size_), d2Pi_(size_), mapFact_(&lessEigen)
    {
      subModel_->setNamespace(getNamespace());
      addParameters_(subModel_->getParameters());
    }
    
    MultinomialFromTransitionModel(const MultinomialFromTransitionModel& fmsm) :
      AbstractParameterAliasable(fmsm),
      AbstractWrappedModel(fmsm),
      subModel_(fmsm.subModel_->clone()),
      size_(fmsm.size_),
      tref_(NumConstants::MINF()), Pij_t(0), dPij_dt(0), d2Pij_dt2(0), Pi_(size_), dPi_(size_), d2Pi_(size_), mapFact_(&lessEigen)
    {
    }

    MultinomialFromTransitionModel& operator=(const MultinomialFromTransitionModel& fmsm)
    {
      AbstractParameterAliasable::operator=(fmsm);

      subModel_ = std::unique_ptr<TransitionModel>(fmsm.subModel_->clone());
      size_ = fmsm.size_;
      Pi_.resize(size_);
      dPi_.resize(size_);
      d2Pi_.resize(size_);
      
      tref_=NumConstants::MINF();
      mapFact_.clear();
      
      return *this;
    }
    
    ~MultinomialFromTransitionModel() {}

    MultinomialFromTransitionModel* clone() const { return new MultinomialFromTransitionModel(*this); }

  public:
    void fireParameterChanged(const ParameterList& parameters)
    {
      AbstractParameterAliasable::fireParameterChanged(parameters);
      if (getModel().matchParametersValues(parameters))
        tref_=NumConstants::MINF();
    }

    const BranchModel& getModel() const
    {
      return *subModel_.get();
    }

    const TransitionModel& getTransitionModel() const
    {
      return *subModel_.get();
    }

    const Eigen::VectorXd& Lik_t    (const Eigen::VectorXd& from, double t) const;
    const Eigen::VectorXd& dLik_dt  (const Eigen::VectorXd& from, double t) const;
    const Eigen::VectorXd& d2Lik_dt2(const Eigen::VectorXd& from, double t) const;
        
    void setFreqFromData(const SequencedValuesContainer& data, double pseudoCount)
    {
      getTransitionModel().setFreqFromData(data, pseudoCount);
    }

    virtual void setFreq(std::map<int, double>& m)
    {  
      getTransitionModel().setFreq(m);
    }

    double getRate() const { return getTransitionModel().getRate(); }

    void setRate(double rate) { return getTransitionModel().setRate(rate); }

    double getInitValue(size_t i, int state) const { return getTransitionModel().getInitValue(i, state); }

    std::string getName() const
    {
      return "MultinomialFrom";
    }

    void addRateParameter()
    {
      getModel().addRateParameter();
      addParameter_(new Parameter(getNamespace() + "rate", getModel().getRate(), Parameter::R_PLUS_STAR));
    }

  private:

    /*
     * @brief Fills the res vector of the likelihoods of the counts
     * given ancestral states & transition matrix.
     *
     */

    
    void compute_Multinomial_(const Eigen::VectorXd& counts) const;
    void compute_dMultinomial_dt_(const Eigen::VectorXd& counts) const;
    void compute_d2Multinomial_dt2_(const Eigen::VectorXd& counts) const;

    static uint factorial(uint n){
      return (n==0 || n==1)?1:factorial(n-1) * n;
    }

    double getFact_(const Eigen::VectorXd& counts) const
    {
      auto it=mapFact_.find(counts);
      if (it==mapFact_.end())
      {
        auto lsto(std::lgamma(counts.sum()+1));
        auto lr((counts.array()+1.).lgamma().sum());
        auto fact=std::exp(lsto-lr);
        
        mapFact_[counts]=fact;
        return fact;
      }
      else
        return it->second;
    }

  };
} // end of namespace bpp.

#endif  // _MULTINOMIAL_FROM_TRANSITION_MODEL_H_
