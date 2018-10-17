//
// File: JointLikelihoodFunction.cpp
// Created by: Keren Halabi
// Created on: Thu Aug 28 13:14 2018
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

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

// for bpp-core
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/App/BppApplication.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Model/MixedSubstitutionModelSet.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>


using namespace bpp;
using namespace std;


class JointLikelihoodFunction:
  public Function,
  public AbstractParametrizable
{

  public:
    enum Hypothesis
    {
        null = 0,
        alternative = 1
    };

    enum OptimizationScope 
    { 
      none=0, 
      onlyCharacter=1, 
      onlySequence=2, 
      both=3 
    };

  private:
    Hypothesis hypothesis_;
    BppApplication* bppml_;
    RHomogeneousTreeLikelihood* characterTreeLikelihood_; 
    RNonHomogeneousMixedTreeLikelihood* sequenceTreeLikelihood_; 
    StochasticMapping* stocMapping_;
    OptimizationScope optimizationScope_;
    double logl_;

  public:

    /* basic functions */

    // base constrcutor - let it get the character model and character data for convenience - but it will create the treeLikelihood instance nad the sotchasticMapping instnace from the constructor
    JointLikelihoodFunction(BppApplication* bppml, const Tree* tree, const VectorSiteContainer* characterData, TransitionModel* characterModel, const VectorSiteContainer* sequenceData, MixedSubstitutionModelSet* sequenceModel, DiscreteDistribution* rDist);

    // copy constructor
    JointLikelihoodFunction(const JointLikelihoodFunction& jlf):
      AbstractParametrizable(""), 
      hypothesis_(jlf.hypothesis_),
      bppml_(jlf.bppml_),
      characterTreeLikelihood_(jlf.characterTreeLikelihood_->clone()),
      sequenceTreeLikelihood_(dynamic_cast<RNonHomogeneousMixedTreeLikelihood*>(jlf.sequenceTreeLikelihood_->clone())),
      stocMapping_(jlf.stocMapping_->clone()),
      optimizationScope_(jlf.optimizationScope_),
      logl_(jlf.logl_)
    { }

    // assignment operator
    JointLikelihoodFunction& operator=(const JointLikelihoodFunction& jlf)
    {
      hypothesis_ = jlf.hypothesis_;
      bppml_ = jlf.bppml_;
      characterTreeLikelihood_ = jlf.characterTreeLikelihood_->clone(); 
      sequenceTreeLikelihood_ = dynamic_cast<RNonHomogeneousMixedTreeLikelihood*>(jlf.sequenceTreeLikelihood_->clone()); 
      stocMapping_ = jlf.stocMapping_->clone();
      optimizationScope_ = jlf.optimizationScope_;
      logl_ = jlf.logl_;
      return *this;
    }

    // destructor
    ~JointLikelihoodFunction();

    // clone operator
    JointLikelihoodFunction* clone() const { return new JointLikelihoodFunction(*this); }
    
    /* basic Function methods */

    /**
     * @brief Set the paramters of the joint likelihood function according ot a given list of parameters
     * 
     * @param pl  A list of parameters whose value should be set inside the joint likelihood function 
     */
    void setParameters(const ParameterList& pl)
    {
      matchParametersValues(pl);
    }
    
    /**
     * @brief Return the last computed log likelihood
     */
    double getValue() const { return logl_; }

    /**
     * @brief Return the pointer to the character likelihood function
     */
    const HomogeneousTreeLikelihood* getCharacterLikelihoodFunction() const { return characterTreeLikelihood_; }

    /**
     * @brief Return the pointer to the sequence likelihood function
     */
    const RNonHomogeneousMixedTreeLikelihood* getSequenceLikelihoodFunction() const { return sequenceTreeLikelihood_; }
 
    /**
     * @brief Computes the value of the joint likelihood function depending on hypothesis - calls either computeNullJointLikelihood or computeAlternativeJointLikelihood
     * 
     * @param pl  A list of parameters whose value should be set inside the joint likelihood function 
     */
    void fireParameterChanged(const ParameterList& pl);

    /* Auxiliary methods */

    /**
     * @brief Set the hypothesis data memeber (either null or alternative)
     */
    void setHypothesis(JointLikelihoodFunction::Hypothesis hypothesis) { hypothesis_ = hypothesis; }

    /**
     * @brief Set the hypothesis data memeber (none / onlyCharacter / onlySequence / both)
     */
    void setOptimizationScope(JointLikelihoodFunction::OptimizationScope optimizationScope) { optimizationScope_ = optimizationScope; }

    /**
     * @brief Defines the tree partition into the two sub-models of the sequence model according to a given history
     * 
     * @param history  A tree in which each node has a state propetly assinged to it, corresponding to its binary character state
     */
    void setPartitionByHistory(const Tree* history);

    /**
     * @brief replaced the instnace of sequwnce tree likelihood so that the updated instance will consider the new character history
     * 
     * @param history  A tree in which each node has a state propetly assinged to it, corresponding to its binary character state
     */
    void updatesequenceTreeLikelihood(const Tree* history);

    /**
     * @brief Computes the joint likelihood function given the null model (that is, no change in selective pressure along the phylogeny, regardless of the character history)
     */
    void computeNullJointLikelihood();
    
    /**
     * @brief Computes the joint likelihood function given the null model (that is, change in selective pressure along the phylogeny depending on the character history)
     */
    void computeAlternativeJointLikelihood();

    /**
     * @brief Optimizes the parameters of the character model
     */
    void optimizeCharacterModel();

    /**
     * @brief Optimizes the parameters of the sequence model
     */
    void optimizeSequenceModel();

    /**
     * @brief Reports the joint model paramerers to std and updates the logl_ data member
     */
    void reportResults();

};