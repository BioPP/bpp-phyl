//
// File: AbstractNonHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Tue Oct 09 16:07 2007
// From file: AbstractHomogeneousTreeLikelihood.h
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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
        unsigned int index_;
        unsigned int nbSites_;

      public:
        ConstNonHomogeneousSiteModelIterator(const SubstitutionModelSet* modelSet, unsigned int nbSites) :
          siteModelDescriptions_(), index_(0), nbSites_(nbSites)
        {
          for (unsigned int i = 0; i < modelSet->getNumberOfModels(); ++i)
            siteModelDescriptions_.push_back(ConstNoPartitionSiteModelDescription(modelSet->getModel(i), modelSet->getNodesWithModel(i)));        
        }

      public:
        ConstSiteModelDescription* next() throw (Exception)
        {
          if (!hasNext())
            throw Exception("AbstractNonHomogeneousTreeLikelihood::ConstHomogeneousSiteModelIterator::next(). No more site in the set.");
          return &siteModelDescriptions_[index_++];
        }

        bool hasNext() const { return index_ < nbSites_; }
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
    std::vector<Node *> nodes_;

    /**
     * @brief An index linking nodes to their id, for faster access than the getNode() method.
     */
    mutable std::map<int, const Node *> idToNode_;
 
    //some values we'll need:
    unsigned int nbSites_,         //the number of sites in the container
                 nbDistinctSites_, //the number of distinct sites
                 nbClasses_,       //the number of rate classes
                 nbStates_,        //the number of states in the alphabet
                 nbNodes_;         //the number of nodes in the tree

    bool _verbose;

    double _minimumBrLen;
    Constraint * brLenConstraint_;

    int root1_, root2_;

  public:
    AbstractNonHomogeneousTreeLikelihood(
      const Tree& tree,
      SubstitutionModelSet* modelSet,
      DiscreteDistribution* rDist,
      bool verbose = true)
      throw (Exception);

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
 
    virtual ~AbstractNonHomogeneousTreeLikelihood();
    
  private:

    /**
     * @brief Method called by constructor.
     */
    void init_(
        const Tree& tree,
        SubstitutionModelSet* modelSet,
        DiscreteDistribution* rDist,
        bool verbose) throw (Exception);

  public:
    
    /**
     * @name The TreeLikelihood interface.
     *
     * Other methods are implemented in the AbstractTreeLikelihood class.
     *
     * @{
     */
    void initialize() throw(Exception);
    
    ParameterList getBranchLengthsParameters() const;
    
    ParameterList getSubstitutionModelParameters() const;
    
    ParameterList getRateDistributionParameters() const
    {
      return AbstractDiscreteRatesAcrossSitesTreeLikelihood::getRateDistributionParameters();
    }

    const SubstitutionModel* getSubstitutionModelForNode(int nodeId) const throw (NodeNotFoundException) 
    {
      return modelSet_->getModelForNode(nodeId);
    }

    SubstitutionModel* getSubstitutionModelForNode(int nodeId) throw (NodeNotFoundException)
    {
      return modelSet_->getModelForNode(nodeId);
    }

    const std::vector<double>& getRootFrequencies(unsigned int siteIndex) const { return rootFreqs_; }
    
    VVVdouble getTransitionProbabilitiesPerRateClass(int nodeId, unsigned int siteIndex) const { return pxy_[nodeId]; }

    ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const
    {
      return new ConstNoPartitionBranchModelIterator(*tree_, modelSet_->getModelForNode(nodeId), nbDistinctSites_);
    }

    ConstSiteModelIterator* getNewSiteModelIterator(unsigned int siteIndex) const
    {
      return new ConstNonHomogeneousSiteModelIterator(modelSet_, nbDistinctSites_);
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
    
    void setSubstitutionModelSet(SubstitutionModelSet* modelSet) throw (Exception);

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
    virtual void applyParameters() throw (Exception);  

    virtual void initBranchLengthsParameters();

    virtual void setMinimumBranchLength(double minimum)
    {
      _minimumBrLen = minimum;
      if(brLenConstraint_ != NULL) delete brLenConstraint_;
      brLenConstraint_ = new IncludingPositiveReal(_minimumBrLen);
      initBranchLengthsParameters();
    }

    virtual double getMinimumBranchLength() const { return _minimumBrLen; }

  protected:
    /**
     * @brief Fill the pxy_, dpxy_ and d2pxy_ arrays for all nodes.
     */
    void computeAllTransitionProbabilities();
    /**
     * @brief Fill the pxy_, dpxy_ and d2pxy_ arrays for one node.
     */
    void computeTransitionProbabilitiesForNode(const Node * node);

};

} //end of namespace bpp.

#endif //_ABSTRACTNONHOMOGENEOUSTREELIKELIHOOD_H_

