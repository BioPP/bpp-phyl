//
// File: AbstractNonHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Tue Oct 09 16:07 2007
// From file: AbstractHomogeneousTreeLikelihood.h
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

#ifndef _ABSTRACTNONHOMOGENEOUSTREELIKELIHOOD_H_
#define _ABSTRACTNONHOMOGENEOUSTREELIKELIHOOD_H_

#include "NonHomogeneousTreeLikelihood.h"
#include "AbstractDiscreteRatesAcrossSitesTreeLikelihood.h"

//From the STL:
#include <memory>

namespace bpp
{

/**
 * @brief Partial implementation for branch non-homogeneous models of the TreeLikelihood interface.
 *
 * This class provides a pointer toward a single substitution model + several utilitary variables.
 */
class AbstractNonHomogeneousTreeLikelihood:
  public virtual NonHomogeneousTreeLikelihood,
  public AbstractDiscreteRatesAcrossSitesTreeLikelihood
{
  public:

    class ConstNonHomogeneousSiteModelIterator :
      public ConstSiteModelIterator
    {
      private:
        std::vector<ConstNoPartitionSiteModelDescription> siteModelDescriptions_;
        size_t index_;
        size_t nbModels_;

      public:
        ConstNonHomogeneousSiteModelIterator(const SubstitutionModelSet* modelSet) :
          siteModelDescriptions_(), index_(0), nbModels_(modelSet->getNumberOfModels())
        {
          for (size_t i = 0; i < nbModels_; ++i)
            siteModelDescriptions_.push_back(ConstNoPartitionSiteModelDescription(modelSet->getModel(i), modelSet->getNodesWithModel(i)));        
        }

      public:
        ConstSiteModelDescription* next()
        {
          if (!hasNext())
            throw Exception("AbstractNonHomogeneousTreeLikelihood::ConstHomogeneousSiteModelIterator::next(). No more site in the set.");
          return &siteModelDescriptions_[index_++];
        }

        bool hasNext() const { return index_ < nbModels_; }
    };

  protected:
    SubstitutionModelSet* modelSet_;
    ParameterList brLenParameters_;
    
    mutable std::map<int, VVVdouble> pxy_;

    mutable std::map<int, VVVdouble> dpxy_;

    mutable std::map<int, VVVdouble> d2pxy_;
        
    std::vector<double> rootFreqs_;
        
    /**
     * @brief Pointer toward all nodes in the tree.
     *
     * The position in the array is the number used in the parameter name.
     * This may be different from the node id, unless you used the resetNodeId method on the input tree.
     */
    std::vector<Node*> nodes_;

    /**
     * @brief An index linking nodes to their id, for faster access than the getNode() method.
     */
    mutable std::map<int, const Node*> idToNode_;
 
    //some values we'll need:
    size_t nbSites_,         //the number of sites in the container
           nbDistinctSites_, //the number of distinct sites
           nbClasses_,       //the number of rate classes
           nbStates_,        //the number of states in the alphabet
           nbNodes_;         //the number of nodes in the tree

    bool verbose_;

    double minimumBrLen_;
    double maximumBrLen_;
    std::shared_ptr<Constraint> brLenConstraint_;

    bool reparametrizeRoot_;
    int root1_, root2_;

  public:
    AbstractNonHomogeneousTreeLikelihood(
      const Tree& tree,
      SubstitutionModelSet* modelSet,
      DiscreteDistribution* rDist,
      bool verbose = true,
      bool reparametrizeRoot = true);

    /**
     * @brief Copy constructor
     *
     * This constructor is to be called by the derived class copy constructor.
     */
    AbstractNonHomogeneousTreeLikelihood(const AbstractNonHomogeneousTreeLikelihood& lik);
    
    /**
     * @brief Assignation operator
     *
     * This operator is to be called by the derived class operator.
     */
    AbstractNonHomogeneousTreeLikelihood& operator=(const AbstractNonHomogeneousTreeLikelihood& lik);
 
    virtual ~AbstractNonHomogeneousTreeLikelihood() {}
    
  private:

    /**
     * @brief Method called by constructor.
     */
    void init_(
        const Tree& tree,
        SubstitutionModelSet* modelSet,
        DiscreteDistribution* rDist,
        bool verbose);

  public:
    
    /**
     * @name The TreeLikelihood interface.
     *
     * Other methods are implemented in the AbstractTreeLikelihood class.
     *
     * @{
     */
    size_t getNumberOfStates() const { return modelSet_->getNumberOfStates(); } 
    
    const std::vector<int>& getAlphabetStates() const { return modelSet_->getAlphabetStates(); } 
    
    int getAlphabetStateAsInt(size_t i) const { return modelSet_->getAlphabetStateAsInt(i); }
  
    std::string getAlphabetStateAsChar(size_t i) const { return modelSet_->getAlphabetStateAsChar(i); }
     
    void initialize();
    
    ParameterList getBranchLengthsParameters() const;
    
    ParameterList getSubstitutionModelParameters() const;
    
    ParameterList getRateDistributionParameters() const
    {
      return AbstractDiscreteRatesAcrossSitesTreeLikelihood::getRateDistributionParameters();
    }

    const TransitionModel* getModelForNode(int nodeId) const 
    {
      return modelSet_->getModelForNode(nodeId);
    }

    TransitionModel* getModelForNode(int nodeId)
    {
      return modelSet_->getModelForNode(nodeId);
    }

    const std::vector<double>& getRootFrequencies(size_t siteIndex) const { return rootFreqs_; }
    
    VVVdouble getTransitionProbabilitiesPerRateClass(int nodeId, size_t siteIndex) const { return pxy_[nodeId]; }

    ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const
    {
      return new ConstNoPartitionBranchModelIterator(modelSet_->getModelForNode(nodeId), nbDistinctSites_);
    }

    ConstSiteModelIterator* getNewSiteModelIterator(size_t siteIndex) const
    {
      return new ConstNonHomogeneousSiteModelIterator(modelSet_);
    }
       
    /** @} */

    /**
     * @name The NonHomogeneousTreeLikelihood interface.
     *
     * Other methods are implemented in the AbstractTreeLikelihood class.
     *
     * @{
     */
    const SubstitutionModelSet* getSubstitutionModelSet() const { return modelSet_; }
    
    SubstitutionModelSet* getSubstitutionModelSet() { return modelSet_; }
    
    void setSubstitutionModelSet(SubstitutionModelSet* modelSet);

    ParameterList getRootFrequenciesParameters() const
    {
      return modelSet_->getRootFrequenciesParameters();
    }
    /** @} */
    
  public: //Specific methods:

    /**
     * @brief This builds the <i>parameters</i> list from all parametrizable objects,
     * <i>i.e.</i> substitution model, rate distribution and tree.
     */
    virtual void initParameters();

    /**
     * @brief All parameters are stored in a parameter list.
     * This function apply these parameters to the substitution model,
     * to the rate distribution and to the branch lengths.
     */
    virtual void applyParameters();  

    virtual void initBranchLengthsParameters(bool verbose = true);

    virtual void setMinimumBranchLength(double minimum)
    {
      if (minimum > maximumBrLen_)
        throw Exception("AbstractNonHomogeneousTreeLikelihood::setMinimumBranchLength. Minimum branch length sould be lower than the maximum one: " + TextTools::toString(maximumBrLen_));
      minimumBrLen_ = minimum;
      brLenConstraint_=std::make_shared<IntervalConstraint>(minimumBrLen_, maximumBrLen_, true, true);
      initBranchLengthsParameters();
    }

    virtual void setMaximumBranchLength(double maximum)
    {
      if (maximum < minimumBrLen_)
        throw Exception("AbstractNonHomogeneousTreeLikelihood::setMaximumBranchLength. Maximum branch length sould be higher than the minimum one: " + TextTools::toString(minimumBrLen_));
      maximumBrLen_ = maximum;
      brLenConstraint_=std::make_shared<IntervalConstraint>(minimumBrLen_, maximumBrLen_, true, true);
      initBranchLengthsParameters();
    }

    virtual double getMinimumBranchLength() const { return minimumBrLen_; }
    virtual double getMaximumBranchLength() const { return maximumBrLen_; }


  protected:
    /**
     * @brief Fill the pxy_, dpxy_ and d2pxy_ arrays for all nodes.
     */
    virtual void computeAllTransitionProbabilities();
    /**
     * @brief Fill the pxy_, dpxy_ and d2pxy_ arrays for one node.
     */
    virtual void computeTransitionProbabilitiesForNode(const Node * node);

};

} //end of namespace bpp.

#endif //_ABSTRACTNONHOMOGENEOUSTREELIKELIHOOD_H_

