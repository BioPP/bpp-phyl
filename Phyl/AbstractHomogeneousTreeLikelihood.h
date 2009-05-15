//
// File: AbstractHomogeneousTreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Thr Dec 23 12:03 2004
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

#ifndef _ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_
#define _ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_

#include "AbstractDiscreteRatesAcrossSitesTreeLikelihood.h"
#include "HomogeneousTreeLikelihood.h"

namespace bpp
{

/**
 * @brief Partial implementation for homogeneous model of the TreeLikelihood interface.
 *
 * This class provides a pointer toward a single substitution model + several utilitary variables.
 */
class AbstractHomogeneousTreeLikelihood:
  public virtual HomogeneousTreeLikelihood,
  public AbstractDiscreteRatesAcrossSitesTreeLikelihood
{
	protected:
		SubstitutionModel * model_;
		ParameterList brLenParameters_;
		
		mutable map<int, VVVdouble> pxy_;

		mutable map<int, VVVdouble> dpxy_;

		mutable map<int, VVVdouble> d2pxy_;

    vector<double> rootFreqs_;
				
		/**
		 * @brief Pointer toward all nodes in the tree.
     *
     * The position in the array is the number used in the parameter name.
     * This may be different from the node id, unless you used the resetNodeId method on the input tree.
 		 */
		vector<Node *> nodes_;

		//some values we'll need:
		unsigned int nbSites_,         //the number of sites in the container
                 nbDistinctSites_, //the number of distinct sites
		             nbClasses_,       //the number of rate classes
		             nbStates_,        //the number of states in the alphabet
		             nbNodes_;         //the number of nodes in the tree

    bool verbose_;

    double minimumBrLen_;
    Constraint * brLenConstraint_;

	public:
		AbstractHomogeneousTreeLikelihood(
			const Tree & tree,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
      bool checkRooted = true,
			bool verbose = true)
      throw (Exception);

    /**
     * @brief Copy constructor
     *
     * This constructor is to be called by the derived class copy constructor.
     */
    AbstractHomogeneousTreeLikelihood(const AbstractHomogeneousTreeLikelihood & lik);
    
    /**
     * @brief Assignation operator
     *
     * This operator is to be called by the derived class operator.
     */
    AbstractHomogeneousTreeLikelihood & operator=(const AbstractHomogeneousTreeLikelihood & lik);
 
		virtual ~AbstractHomogeneousTreeLikelihood();
		
	private:

    /**
     * @brief Method called by constructor.
     */
    void _init(const Tree & tree,
			SubstitutionModel * model,
			DiscreteDistribution * rDist,
      bool checkRooted,
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
    /**
     * @brief Get the substitution model associated to a given node.
     *
     * In the homogeneous case, this function is aliased by the getSubstitutionModel function.
     *
     * @param nodeId The id of the request node.
     * @return A pointer toward the corresponding model.
     * @throw NodeNotFoundException This exception may be thrown if the node is not found (depending on the implementation).
     */
    const SubstitutionModel * getSubstitutionModelForNode(int nodeId) const throw (NodeNotFoundException) { return model_; }

    /**
     * @brief Get the substitution model associated to a given node.
     *
     * In the homogeneous case, this function is aliased by the getSubstitutionModel function.
     *
     * @param nodeId The id of the request node.
     * @return A pointer toward the corresponding model.
     * @throw NodeNotFoundException This exception may be thrown if the node is not found (depending on the implementation).
     */
    SubstitutionModel * getSubstitutionModelForNode(int nodeId) throw (NodeNotFoundException) { return model_; }

    vector<double> getRootFrequencies() const { return model_->getFrequencies(); }
    
    const VVVdouble & getTransitionProbabilitiesForNode(const Node* node) const { return pxy_[node->getId()]; }
    /** @} */

		/**
		 * @name The HomogeneousTreeLikelihood interface.
		 *
		 * Other methods are implemented in the AbstractTreeLikelihood class.
		 *
		 * @{
		 */
		const SubstitutionModel * getSubstitutionModel() const { return model_; }
		
		SubstitutionModel * getSubstitutionModel() { return model_; }
		
    void setSubstitutionModel(SubstitutionModel * model) throw (Exception);
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
      minimumBrLen_ = minimum;
      if(brLenConstraint_ != NULL) delete brLenConstraint_;
      brLenConstraint_ = new IncludingPositiveReal(minimumBrLen_);
      initBranchLengthsParameters();
    }

    virtual double getMinimumBranchLength() const { return minimumBrLen_; }

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

#endif //_ABSTRACTHOMOGENEOUSTREELIKELIHOOD_H_

