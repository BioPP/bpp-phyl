//
// File: SubstitutionProcess.h
// Created by: Julien Dutheil
// Created on: Tue May 15 13:11 2012
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

#ifndef _SUBSTITUTIONPROCESS_H_
#define _SUBSTITUTIONPROCESS_H_

#include "ParametrizableTree.h"
#include "ComputingTree.h"
#include "../Model/SubstitutionModel.h"

//From bpp-core:
#include <Bpp/Numeric/ParameterAliasable.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

//From the STL:
#include <memory>

namespace bpp
{

/**
 * @brief This interface describes the substitution process along the tree and sites of the alignment.
 *
 * It main purpose is to provide the necessary calculus for each branch-site-model class combination,
 * such as Markov generator and transition probabilities.
 * These are typically provided by a SubstitutionModel class, applied in various combination along the
 * tree (eg non-homogeneous models) and alignment (eg partition models).
 * The so-called "model class" refers to mixture models.
 *
 * An instance of the SubstitutionProcess class is always associated to an instance of a ParametrizableTree.
 * The substitution process is in charge of computing the transition probabilities for each branch of the tree,
 * and to update them in case a parameter (including branch length) changes.
 *
 * As several branches and sites can share the same generator/transition probabilities, calling these values
 * and therefore performing the underlying calculation for each branch-site can result in a very unefficient
 * code. Therefore, "model iterators" are provided, to allow smarter loops over the model structure,
 * by minimizing the amount of computation to be done.
 * Such iterator allow to loop over branches and sites in a clever way, through both directions, but do not
 * perform any calculations.
 * These are achieved through calls to the corresponding SubstitutionProcess class.
 */
  class SubstitutionProcess :
    public virtual ParameterAliasable
  {
  public:
    virtual SubstitutionProcess* clone() const = 0;

  public:
    virtual bool isCompatibleWith(const SiteContainer& data) const = 0;

    virtual const TreeTemplate<Node>& getTree() const = 0;
  
    virtual const ParametrizableTree& getParametrizableTree() const = 0;

    virtual bool hasDerivableParameter(const std::string& name) const = 0;

    virtual size_t getNumberOfClasses() const = 0;
  
    virtual size_t getNumberOfStates() const = 0;

    /**
     * @brief Get the substitution model corresponding to a certain
     * branch, site pattern, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */

    virtual const SubstitutionModel& getSubstitutionModel(int nodeId, size_t classIndex) const = 0;

    /**
     * @brief Get a pointer to the rate distribution (or null if there
     * is no rate distribution).
     *
     */

    virtual const DiscreteDistribution* getRateDistribution() const = 0;

    /**
     * @brief Methods to retrieve the parameters of specific objects.
     *
     */
  
    virtual ParameterList getSubstitutionModelParameters(bool independent) const = 0;

    virtual ParameterList getRateDistributionParameters(bool independent) const = 0;

    virtual ParameterList getRootFrequenciesParameters(bool independent) const = 0;

    virtual ParameterList getBranchLengthParameters(bool independent) const = 0;

    virtual ParameterList getDerivableParameters() const = 0;

    virtual ParameterList getNonDerivableParameters() const = 0;

    /**
     * @brief Get the transition probabilities corresponding to a
     * certain branch, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */
  
    virtual const Matrix<double>& getTransitionProbabilities(int nodeId, size_t classIndex) const = 0;
 
    /**
     * @brief Get the first order derivatives of the transition
     * probabilities according to time, corresponding to a certain
     * branch, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */
    virtual const Matrix<double>& getTransitionProbabilitiesD1(int nodeId, size_t classIndex) const = 0;
 
    /**
     * @brief Get the second order derivatives of the transition
     * probabilities according to time, corresponding to a certain
     * branch, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */
    virtual const Matrix<double>& getTransitionProbabilitiesD2(int nodeId, size_t classIndex) const = 0;
 
    /**
     * @brief Get the generator corresponding to a certain branch, and
     * model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */
    virtual const Matrix<double>& getGenerator(int nodeId, size_t classIndex) const = 0;

    /**
     * @brief Get the values of the frequencies for each state in the
     * alphabet at the root node.
     *
     * For reversible models, these are the equilibrium frequencies.
     * For non-reversible models, these usually are distinct parameters.
     *
     * @return A vector with ancestral frequencies for each state in the alphabet;
     */
    virtual const std::vector<double>& getRootFrequencies() const = 0;

    virtual const FrequenciesSet* getRootFrequenciesSet() const = 0;

    /**
     * This method is used to initialize likelihoods in reccursions.
     * It typically sends 1 if i = state, 0 otherwise, where
     * i is one of the possible states of the alphabet allowed in the model
     * and state is the observed state in the considered sequence/site.
     *
     * @param i the index of the state in the model.
     * @param state An observed state in the sequence/site.
     * @return 1 or 0 depending if the two states are compatible.
     * @throw BadIntException if states are not allowed in the associated alphabet.
     * @see getStates();
     * @see SubstitutionModel
     */
    virtual double getInitValue(size_t i, int state) const throw (BadIntException) = 0;

    /**
     * @return The probability associated to the given model class.
     *
     * @param classIndex The model class index.
     */
    virtual double getProbabilityForModel(size_t classIndex) const = 0;

    /**
     * @return The vector of the probabilities associated to the
     * classes.
     *
     */
  
    virtual Vdouble getClassProbabilities() const = 0;

    /**
     * @return The substitution rate associated to the given model class.
     *
     * @param classIndex The model class index.
     */
    virtual double getRateForModel(size_t classIndex) const = 0;

    /**
     * @brief Tell if the transition probabilities have changed after the last call to setParameters().
     * @return True if transition probabilities have changed.
     */

    //virtual bool transitionProbabilitiesHaveChanged() const = 0; Not sure we need that anymore...

    /**
     * @brief A virtual method to retrieve the ComputingTree defined in
     * inheriting classes.
     *
     */
  
    virtual const ComputingTree& getComputingTree() const = 0;

    virtual ComputingTree& getComputingTree() = 0;

  };

} // end namespace bpp

#endif // _SUBSTITUTIONPROCESS_H_

