//
// File: TreeLikelihood.h
// Created by: Julien Dutheil
// Created on: Fri Oct 17 17:36:44 2003
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

#ifndef _TREELIKELIHOOD_H_
#define _TREELIKELIHOOD_H_

#include "Node.h"
#include "Tree.h"
#include "SubstitutionModel.h"

// From NumCalc:
#include <NumCalc/ParameterList.h>
#include <NumCalc/Parametrizable.h>
#include <NumCalc/Functions.h>
#include <NumCalc/VectorTools.h>

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/SiteContainer.h>

namespace bpp
{

/**
 * @brief The TreeLikelihood interface.
 *
 * This interface defines the methods needed for computing the likelihood
 * of a phylogenetic tree, given a dataset.
 */ 
class TreeLikelihood:
  public virtual DerivableSecondOrder
{
  public:
    TreeLikelihood() {}
    virtual ~TreeLikelihood() {}

#ifndef NO_VIRTUAL_COV
    TreeLikelihood * clone() const = 0;
#endif

  public:

    /**
     * @brief Set the dataset for which the likelihood must be evaluated.
     *
     * @param sites The data set to use.
     */
    virtual void setData(const SiteContainer & sites) = 0;
    
    /**
     * @brief Get the dataset for which the likelihood must be evaluated.
     *
     * @return A pointer toward the site container where the sequences are stored.
     */
    virtual const SiteContainer * getData() const = 0;

    /**
     * @brief Init the likelihood object.
     *
     * This method is used to initialize all parameters.
     * It is typically called after the constructor and the setData method.
     * It contains virtual methods that can't be called in the constructor.
     * @throw Exception if something bad happened, for instance if no data are associated to the likelihood function.
     */
    virtual void initialize() throw (Exception) = 0;

    /**
     * @return 'true' is the likelihood function has been initialized.
     */
    virtual bool isInitialized() const = 0;

    /**
     * @brief Get the likelihood for a site.
     *
     * @param site The site index to analyse.
     * @return The likelihood for site <i>site</i>.
     */
    virtual double getLikelihoodForASite(unsigned int site) const = 0;

    /**
     * @brief Get the logarithm of the likelihood for a site.
     *
     * @param site The site index to analyse.
     * @return The logarithm of the likelihood for site <i>site</i>.
     */
    virtual double getLogLikelihoodForASite(unsigned int site) const = 0;

    /**
     * @brief Get the likelihood for a site and for a state.
     *
     * @param site The site index to analyse.
     * @param state The state to consider.
     * @return The likelihood for site <i>site</i> and state <i>state</i>.
     */
    virtual double getLikelihoodForASiteForAState(unsigned int site, int state) const = 0;

    /**
     * @brief Get the logarithm of the likelihood for a site and for a state.
     *
     * @param site The site index to analyse.
     * @param state The state to consider.
     * @return The logarithm of the likelihood for site <i>site</i> and state <i>state</i>.
     */
    virtual double getLogLikelihoodForASiteForAState(unsigned int site, int state) const = 0;

    /**
     * @brief Get the likelihood for each site.
     *
     * @return A vector with all likelihoods for each site.
     */
    virtual Vdouble getLikelihoodForEachSite() const = 0;

    /**
     * @brief Get the logarithm of the likelihood for each site.
     *
     * @return A vector with all log likelihoods for each site.
     */
    virtual Vdouble getLogLikelihoodForEachSite() const = 0;

    /**
     * @brief Get the likelihood for each site and for each state.
     *
     * @return A 2d vector with all likelihoods for each site and for each state.
     */
    virtual VVdouble getLikelihoodForEachSiteForEachState() const = 0;

    /**
     * @brief Get the logarithm of the likelihood for each site and for each state.
     *
     * @return A 2d vector with all log likelihoods for each site and for each state.
     */
    virtual VVdouble getLogLikelihoodForEachSiteForEachState() const = 0;
    
    /**
     * @brief Get the likelihood for the whole dataset.
     *
     * @return The likelihood of the dataset.
     */
    virtual double getLikelihood() const = 0;

    /**
     * @brief Get the logarithm of the likelihood for the whole dataset.
     *
     * @return The logarithm of the likelihood of the dataset.
     */
    virtual double getLogLikelihood() const = 0;
  
    /**
     * @brief Get the tree (topology and branch lengths).
     *
     * @return The tree of this TreeLikelihood object.
      */
    virtual const Tree * getTree() const = 0;

    /**
     * @brief Get the number of sites in the dataset.
     *
     * @return the number of sites in the dataset.
     */
    virtual unsigned int getNumberOfSites() const = 0;

    /**
     * @brief Get the number of states in the alphabet associated to the dataset.
     *
     * @return the number of states in the alphabet associated to the dataset.
     */    
    virtual unsigned int getNumberOfStates() const = 0;
    
    /**
     * @brief Get the alphabet associated to the dataset.
     *
     * @return the alphabet associated to the dataset.
     */    
    virtual const Alphabet * getAlphabet() const = 0;
    
    /**
     * @name Retrieve some particular parameters subsets.
     *
     * @{
     */
    
    /**
     * @brief Get the branch lengths parameters.
     *
     * @return A ParameterList with all branch lengths.
     */
    virtual ParameterList getBranchLengthsParameters() const = 0;
    
    /**
     * @brief Get the parameters associated to substitution model(s).
     *
     * @return A ParameterList.
     */
    virtual ParameterList getSubstitutionModelParameters() const = 0;

    /**
     * @brief Get the substitution model associated to a given node.
     *
     * @param nodeId The id of the request node.
     * @return A pointer toward the corresponding model.
     * @throw NodeNotFoundException This exception may be thrown if the node is not found (depending on the implementation).
     */
    virtual const SubstitutionModel * getSubstitutionModelForNode(int nodeId) const throw (NodeNotFoundException) = 0;

    /**
     * @brief Get the substitution model associated to a given node.
     *
     * @param nodeId The id of the request node.
     * @return A pointer toward the corresponding model.
     * @throw NodeNotFoundException This exception may be thrown if the node is not found (depending on the implementation).
     */
    virtual SubstitutionModel * getSubstitutionModelForNode(int nodeId) throw (NodeNotFoundException) = 0;

    /**
     * @brief Get the values of the frequencies for each state in the alphabet at the root node.
     *
     * For reversible models, these are the equilibrium frequencies.
     * For non-reversible models, these usually are distinct parameters.
     *
     * @return A vector with ancestral frequencies for each state in the alphabet;
     */
    virtual vector<double> getRootFrequencies() const = 0;
    
    /** @} */

    /**
     * @brief Tell if derivatives must be computed.
     *
     * This methods calls the enableFirstOrderDerivatives and enableSecondOrderDerivatives.
     *
     * @param yn Yes or no.
     */
    virtual void enableDerivatives(bool yn) = 0;

    /**
     * @brief All derivable parameters.
     *
     * Usually, this contains all branch lengths parameters.
     *
     * @return A ParameterList.
     */
    virtual ParameterList getDerivableParameters() const = 0;

    /**
     * @brief All non derivable parameters.
     *
     * Usually, this contains all substitution model parameters and rate distribution.
     *
     * @return A ParameterList.
     */
    virtual ParameterList getNonDerivableParameters() const = 0;

    /**
     * @brief Retrieves all Pij(t) for a particular branch, defined by the upper node.
     *
     * These intermediate results may be used by other methods.
     *
     * @param node The node defining the branch of interest.
     * @return An array of dimension 3, where a[c][x][y] is the probability of substituting from x to y while being in rate class c.
     */
    virtual const VVVdouble & getTransitionProbabilitiesForNode(const Node* node) const = 0;
};

} //end of namespace bpp.

#endif  //_TREELIKELIHOOD_H_

