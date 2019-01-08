//
// File: DistanceEstimation.h
// Created by: Julien Dutheil
//             Vincent Ranwez
// Created on: Wed jun 08 10:39 2005
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

#ifndef _DISTANCEESTIMATION_H_
#define _DISTANCEESTIMATION_H_

#include "../Model/SubstitutionModel.h"
#include "../Likelihood/AbstractTreeLikelihood.h"
#include "../Likelihood/DRHomogeneousTreeLikelihood.h"
#include "../PseudoNewtonOptimizer.h"

#include <Bpp/Clonable.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Function/SimpleMultiDimensions.h>
#include <Bpp/Numeric/Function/MetaOptimizer.h>

// From bpp-seq:
#include <Bpp/Seq/Container/AlignedValuesContainer.h>

// From the STL:
#include <memory>

namespace bpp
{

/**
 * @brief This class is a simplified version of DRHomogeneousTreeLikelihood for 2-Trees.
 */
  class TwoTreeLikelihood:
    public AbstractDiscreteRatesAcrossSitesTreeLikelihood  
  {
  private:
    std::shared_ptr<AlignedValuesContainer> shrunkData_;
    std::vector<std::string> seqnames_;
    TransitionModel* model_;
    ParameterList brLenParameters_;
    
    mutable VVVdouble pxy_;

    mutable VVVdouble dpxy_;

    mutable VVVdouble d2pxy_;

    /**
     * @brief As previous, but for the global container.
     *
     * The size of this vector is equal to the number of sites in the container,
     * each element corresponds to a site in the container and points to the
     * corresponding column in the likelihood array of the root node.
     * If the container contains no repeated site, there will be a strict
     * equivalence between each site and the likelihood array of the root node.
     * However, if this is not the case, some pointers may point toward the same
     * element in the likelihood array.
     */
    std::vector<size_t> rootPatternLinks_;

    /**
     * @brief The frequency of each site.
     */
    std::vector<unsigned int> rootWeights_;

    //some values we'll need:
    size_t nbSites_,         //the number of sites in the container
      nbClasses_,       //the number of rate classes
      nbStates_,        //the number of states in the alphabet
      nbDistinctSites_; //the number of distinct sites in the container

    mutable VVVdouble rootLikelihoods_;
    mutable VVdouble rootLikelihoodsS_;
    mutable Vdouble rootLikelihoodsSR_;
    mutable Vdouble dLikelihoods_;
    mutable Vdouble d2Likelihoods_;
    mutable VVdouble leafLikelihoods1_, leafLikelihoods2_;
  
    double minimumBrLen_;
    std::shared_ptr<Constraint> brLenConstraint_;
    double brLen_;

  public:
    TwoTreeLikelihood(
      const std::string& seq1, const std::string& seq2,  
      const AlignedValuesContainer& data,
      TransitionModel* model,
      DiscreteDistribution* rDist,
      bool verbose);

    TwoTreeLikelihood(const TwoTreeLikelihood& lik);
    
    TwoTreeLikelihood& operator=(const TwoTreeLikelihood& lik);

    TwoTreeLikelihood* clone() const { return new TwoTreeLikelihood(*this); } 

    virtual ~TwoTreeLikelihood();

  public:

    /**
     * @name The TreeLikelihood interface.
     *
     * Other methods are implemented in the AbstractTreeLikelihood class.
     *
     * @{
     */
    size_t getNumberOfStates() const { return model_->getNumberOfStates(); } 
    
    const std::vector<int>& getAlphabetStates() const { return model_->getAlphabetStates(); } 
    
    int getAlphabetStateAsInt(size_t i) const { return model_->getAlphabetStateAsInt(i); }
  
    std::string getAlphabetStateAsChar(size_t i) const { return model_->getAlphabetStateAsChar(i); }
    
    TreeLikelihoodData* getLikelihoodData()
    {
      throw NotImplementedException("TwoTreeLikelihood::getLikelihoodData.");
    }
    const TreeLikelihoodData* getLikelihoodData() const
    {
      throw NotImplementedException("TwoTreeLikelihood::getLikelihoodData.");
    }  
    double getLikelihood() const;
    double getLogLikelihood() const;
    double getLikelihoodForASite (size_t site) const;
    double getLogLikelihoodForASite(size_t site) const;
    ParameterList getBranchLengthsParameters() const;
    ParameterList getSubstitutionModelParameters() const;
  
    TransitionModel* getModelForSite(int nodeId, size_t siteIndex) { return model_; }
  
    const TransitionModel* getModelForSite(int nodeId, size_t siteIndex) const { return model_; }

    const std::vector<double>& getRootFrequencies(size_t siteIndex) const { return model_->getFrequencies(); }
    size_t getSiteIndex(size_t site) const { return rootPatternLinks_[site]; }
    /**
     * @brief This method is not applicable for this object.
     */
    VVVdouble getTransitionProbabilitiesPerRateClass(int nodeId, size_t siteIndex) const { return pxy_; }
    void setData(const AlignedValuesContainer& sites) {}
    void initialize();
    /** @} */

    /**
     * @name The DiscreteRatesAcrossSites interface implementation:
     *
     * @{
     */
    double getLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;
    double getLogLikelihoodForASiteForARateClass(size_t site, size_t rateClass) const;
    double getLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;
    double getLogLikelihoodForASiteForARateClassForAState(size_t site, size_t rateClass, int state) const;
    /** @} */

    const TransitionModel* getModel() const { return model_; }

    /**
     * @brief Get the substitution model used for the computation.
     *
     * @return A pointer toward the substitution model of this instance.
     */
    TransitionModel* getModel() { return model_; }

    ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const
    {
      throw NotImplementedException("TwoTreeLikelihood::getNewBranchSiteModelIterator. This class does not (yet) provide support for partition models.");
    }

    ConstSiteModelIterator* getNewSiteModelIterator(size_t siteIndex) const
    {
      throw NotImplementedException("TwoTreeLikelihood::getNewSiteModelIterator. This class is for inner use only and does not provide site model iterators.");
    }

    

    /**
     * @brief Implements the Function interface.
     *
     * Update the parameter list and call the applyParameters() method.
     * Then compute the likelihoods at each node (computeLikelihood() method)
     * and call the getLogLikelihood() method.
     *
     * If a subset of the whole parameter list is passed to the function,
     * only these parameters are updated and the other remain constant (i.e.
     * equal to their last value).
     *
     * @param parameters The parameter list to pass to the function.
     */
    void setParameters(const ParameterList& parameters);
    double getValue() const;
    
    /**
     * @name DerivableFirstOrder interface.
     *
     * @{
     */
    double getFirstOrderDerivative(const std::string& variable) const;
    /** @{ */

    /**
     * @name DerivableSecondOrder interface.
     *
     * @{
     */
    double getSecondOrderDerivative(const std::string& variable) const;
    double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const { return 0; } // Not implemented for now.
    /** @} */

    virtual void initBranchLengthsParameters();

    virtual void setMinimumBranchLength(double minimum)
    {
      minimumBrLen_ = minimum;
      brLenConstraint_ = std::make_shared<IntervalConstraint>(IntervalConstraint(1, minimumBrLen_, true));
      initBranchLengthsParameters();
    }

    virtual double getMinimumBranchLength() const { return minimumBrLen_; }

  protected:
    
    /**
     * @brief This method initializes the leaves according to a sequence container.
     *
     * Here the container shrunkData_ is used.
     * Likelihood is set to 1 for the state corresponding to the sequence site,
     * otherwise it is set to 0.
     *
     * The two likelihood arrays are initialized according to alphabet
     * size and sequences length, and filled with 1.
     *
     * NB: This method is recursive.
     *
     * @param sequences The sequence container to use.
     */
    virtual void initTreeLikelihoods(const SequencedValuesContainer & sequences);

    void fireParameterChanged(const ParameterList & params);
    virtual void computeTreeLikelihood();
    virtual void computeTreeDLikelihood();
    virtual void computeTreeD2Likelihood();
    /**
     * @brief This builds the <i>parameters</i> list from all parametrizable objects,
     * <i>i.e.</i> substitution model, rate distribution and tree.
     */
    virtual void initParameters();

    /**
     * @brief All parameters are stores in a parameter list.
     *
     * This function apply these parameters to the substitution model,
     * to the rate distribution and to the branch lengths.
     */
    virtual void applyParameters();  

  };

/**
 * @brief Estimate a distance matrix from sequence data, according to a given model.
 *
 * By default, the parameters of the model are fixed to there given values.
 * It is possible to estimate one or several parameters by setting them with the
 * setAdditionalParameters() method.
 * Parameters will be estimated separately for each pair of sequence.
 *
 * For now it is not possible to retrieve estimated values.
 * You'll have to specify a 'profiler' to the optimizer and then look at the file
 * if you want to do so.
 */
  class DistanceEstimation:
    public virtual Clonable
  {
  private:
    std::unique_ptr<TransitionModel> model_;
    std::unique_ptr<DiscreteDistribution> rateDist_;
    const AlignedValuesContainer* sites_;
    DistanceMatrix* dist_;
    Optimizer* optimizer_;
    MetaOptimizer* defaultOptimizer_;
    size_t verbose_;
    ParameterList parameters_;

  public:
  
    /**
     * @brief Create a new DistanceEstimation object according to a given substitution model and a rate distribution.
     *
     * This instance will own the model and distribution, and will take car of their recopy and destruction.
     *
     * @param model    The substitution model to use.
     * @param rateDist The discrete rate distribution to use.
     * @param verbose  The verbose level:
     *  - 0=Off,
     *  - 1=one * by row computation
     *  - 2=one * by row computation and one . by column computation
     *  - 3=2 + optimization verbose enabled
     *  - 4=3 + likelihood object verbose enabled
     */
    DistanceEstimation(
      TransitionModel* model,
      DiscreteDistribution* rateDist,
      size_t verbose = 1) :
      model_(model),
      rateDist_(rateDist),
      sites_(0),
      dist_(0),
      optimizer_(0),
      defaultOptimizer_(0),
      verbose_(verbose),
      parameters_()
    {
      init_();
    }
  
    /**
     * @brief Create a new DistanceEstimation object and compute distances
     * according to a given substitution model and a rate distribution.
     *
     * This instance will own the model and distribution, and will take car of their recopy and destruction.
     *
     * @param model    The substitution model to use.
     * @param rateDist The discrete rate distribution to use.
     * @param sites    The sequence data.
     * @param verbose  The verbose level:
     *  - 0=Off,
     *  - 1=one * by row computation
     *  - 2=one * by row computation and one . by column computation
     *  - 3=2 + optimization verbose enabled
     *  - 4=3 + likelihood object verbose enabled
     *  @param computeMat if true the computeMatrix() method is called.
     */
    DistanceEstimation(
      TransitionModel* model,
      DiscreteDistribution* rateDist,
      const AlignedValuesContainer* sites,
      size_t verbose = 1,
      bool computeMat = true) :
      model_(model),
      rateDist_(rateDist),
      sites_(sites),
      dist_(0),
      optimizer_(0),
      defaultOptimizer_(0),
      verbose_(verbose),
      parameters_()
    {
      init_();
      if(computeMat) computeMatrix();
    }
    
    /**
     * @brief Copy constructor.
     *
     * Only the distance matrix is hard-copied, if there is one.
     *
     * @param distanceEstimation The object to copy.
     */
    DistanceEstimation(const DistanceEstimation& distanceEstimation):
      model_(distanceEstimation.model_->clone()),
      rateDist_(distanceEstimation.rateDist_->clone()),
      sites_(distanceEstimation.sites_),
      dist_(0),
      optimizer_(dynamic_cast<Optimizer *>(distanceEstimation.optimizer_->clone())),
      defaultOptimizer_(dynamic_cast<MetaOptimizer *>(distanceEstimation.defaultOptimizer_->clone())),
      verbose_(distanceEstimation.verbose_),
      parameters_(distanceEstimation.parameters_)
    {
      if(distanceEstimation.dist_ != 0)
        dist_ = new DistanceMatrix(*distanceEstimation.dist_);
      else
        dist_ = 0;
    }

    /**
     * @brief Assigment operator.
     *
     * Only the distance matrix is hard-copied, if there is one.
     * 
     * @param distanceEstimation The object to copy.
     * @return A reference toward this object.
     */
    DistanceEstimation& operator=(const DistanceEstimation& distanceEstimation)
    {
      model_.reset(distanceEstimation.model_->clone());
      rateDist_.reset(distanceEstimation.rateDist_->clone());
      sites_      = distanceEstimation.sites_;
      if (distanceEstimation.dist_ != 0)
        dist_     = new DistanceMatrix(*distanceEstimation.dist_);
      else
        dist_     = 0;
      optimizer_  = dynamic_cast<Optimizer *>(distanceEstimation.optimizer_->clone());
      // _defaultOptimizer has already been initialized since the default constructor has been called.
      verbose_    = distanceEstimation.verbose_;
      parameters_ = distanceEstimation.parameters_;
      return *this;
    }

    virtual ~DistanceEstimation()
    {
      if (dist_) delete dist_;
      delete defaultOptimizer_;
      delete optimizer_;
    }

    DistanceEstimation* clone() const { return new DistanceEstimation(*this); }
    
  private:
    void init_()
    {
      MetaOptimizerInfos* desc = new MetaOptimizerInfos();
      std::vector<std::string> name;
      name.push_back("BrLen");
      desc->addOptimizer("Branch length", new PseudoNewtonOptimizer(0), name, 2, MetaOptimizerInfos::IT_TYPE_FULL);
      ParameterList tmp = model_->getParameters();
      tmp.addParameters(rateDist_->getParameters());
      desc->addOptimizer("substitution model and rate distribution", new SimpleMultiDimensions(0), tmp.getParameterNames(), 0, MetaOptimizerInfos::IT_TYPE_STEP);
      defaultOptimizer_ = new MetaOptimizer(0, desc);
      defaultOptimizer_->setMessageHandler(0);
      defaultOptimizer_->setProfiler(0);
      defaultOptimizer_->getStopCondition()->setTolerance(0.0001);
      optimizer_ = dynamic_cast<Optimizer*>(defaultOptimizer_->clone());
    }

  public:

    /**
     * @brief Perform the distance computation.
     *
     * Result can be called by the getMatrix() method.
     *
     * @throw NullPointerException if at least one of the model,
     * rate distribution or data are not initialized.
     */
    void computeMatrix();
    
    /**
     * @brief Get the distance matrix.
     *
     * @return A pointer toward the computed distance matrix.
     */
    DistanceMatrix* getMatrix() const { return dist_ == 0 ? 0 : new DistanceMatrix(*dist_); }

    bool hasModel() const { return model_.get(); }

    const TransitionModel& getModel() const
    {
      if (hasModel())
        return *model_;
      else
        throw Exception("DistanceEstimation::getSubstitutionModel(). No model assciated to this instance.");
    }

    void resetSubstitutionModel(TransitionModel* model = 0) { model_.reset(model); }

    bool hasRateDistribution() const { return rateDist_.get(); }

    const DiscreteDistribution& getRateDistribution() const {
      if (hasRateDistribution())
        return *rateDist_;
      else
        throw Exception("DistanceEstimation::getRateDistribution(). No rate distribution assciated to this instance.");
    }

    void resetRateDistribution(DiscreteDistribution* rateDist = 0) { rateDist_.reset(rateDist); }

    void setData(const AlignedValuesContainer* sites) { sites_ = sites; }
    const AlignedValuesContainer* getData() const { return sites_; }
    void resetData() { sites_ = 0; }
    
    void setOptimizer(const Optimizer * optimizer)
    { 
      if (optimizer_) delete optimizer_;
      optimizer_ = dynamic_cast<Optimizer *>(optimizer->clone());
    }
    const Optimizer* getOptimizer() const { return optimizer_; }
    Optimizer* getOptimizer() { return optimizer_; }
    void resetOptimizer() { optimizer_ = dynamic_cast<Optimizer*>(defaultOptimizer_->clone()); }

    /**
     * @brief Specify a list of parameters to be estimated.
     *
     * Parameters will be estimated separately for each distance.
     *
     * @param parameters A list of parameters to estimate.
     */
    void setAdditionalParameters(const ParameterList& parameters)
    {
      parameters_ = parameters;
    }

    /**
     * @brief Reset all additional parameters.
     */
    void resetAdditionalParameters()
    {
      parameters_.reset();
    }

    /**
     * @param verbose Verbose level.
     */
    void setVerbose(size_t verbose) { verbose_ = verbose; }
    /**
     * @return Verbose level.
     */
    size_t getVerbose() const { return verbose_; }
  };

} //end of namespace bpp.

#endif //_DISTANCEESTIMATION_H_

