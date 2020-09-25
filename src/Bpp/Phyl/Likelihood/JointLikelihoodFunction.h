//
// File: JointLikelihoodFunction.h
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
    map<string, double> previousParametersValues_;
    StochasticMapping* stocMapping_;
    OptimizationScope optimizationScope_;
    double logl_;
    bool characterChanged_; // if characterChanged_ == true -> need to compute sequence likelihood again (both null and alternative) + create new character history (only in alternative)
    bool sequenceChanged_;  // if sequenceChanged_ == true -> need to compute sequence likelihood again
    double previousK_;
    double origTreeLength_; // save the original tree length in order to report the induced sequence scaling factor when needed
    string debugDir_;
    bool debug_;
    unsigned int cycleNum_;
  
  protected:

    /**
     * @brief Internal builtin function by which stsrting points are sorted in the optimization procedure
     */
    static bool sortStartingPointsFunction(map<string,double> i,map<string,double> j) { return (i["Overall Log likelihood"]<j["Overall Log likelihood"]); }


  public:

    /* basic functions */

    // base constrcutor - let it get the character model and character data for convenience - but it will create the treeLikelihood instance nad the sotchasticMapping instnace from the constructor
    JointLikelihoodFunction(BppApplication* bppml, const Tree* tree, const VectorSiteContainer* characterData, TransitionModel* characterModel, const VectorSiteContainer* sequenceData, MixedSubstitutionModelSet* sequenceModel, DiscreteDistribution* rDist, bool debug=false);

    // copy constructor
    JointLikelihoodFunction(const JointLikelihoodFunction& jlf):
      AbstractParametrizable(""), 
      hypothesis_(jlf.hypothesis_),
      bppml_(jlf.bppml_),
      characterTreeLikelihood_(jlf.characterTreeLikelihood_->clone()),
      sequenceTreeLikelihood_(dynamic_cast<RNonHomogeneousMixedTreeLikelihood*>(jlf.sequenceTreeLikelihood_->clone())),
      previousParametersValues_(jlf.previousParametersValues_),
      stocMapping_(jlf.stocMapping_->clone()),
      optimizationScope_(jlf.optimizationScope_),
      logl_(jlf.logl_),
      characterChanged_(jlf.characterChanged_),
      sequenceChanged_(jlf.sequenceChanged_),
      previousK_(jlf.previousK_),
      origTreeLength_(jlf.origTreeLength_),
      debugDir_(jlf.debugDir_),
      debug_(jlf.debug_),
      cycleNum_(jlf.cycleNum_)
    {   
    }

    // assignment operator
    JointLikelihoodFunction& operator=(const JointLikelihoodFunction& jlf)
    {
      hypothesis_ = jlf.hypothesis_;
      bppml_ = jlf.bppml_;
      characterTreeLikelihood_ = jlf.characterTreeLikelihood_->clone(); 
      sequenceTreeLikelihood_ = dynamic_cast<RNonHomogeneousMixedTreeLikelihood*>(jlf.sequenceTreeLikelihood_->clone()); 
      previousParametersValues_ = jlf.previousParametersValues_;
      stocMapping_ = jlf.stocMapping_->clone();
      optimizationScope_ = jlf.optimizationScope_;
      logl_ = jlf.logl_;
      characterChanged_ = jlf.characterChanged_;
      sequenceChanged_ = jlf.sequenceChanged_;
      previousK_ = jlf.previousK_;
      origTreeLength_ = jlf.origTreeLength_;
      debugDir_ = jlf.debugDir_;
      debug_ = jlf.debug_;
      cycleNum_ = jlf.cycleNum_;
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
    HomogeneousTreeLikelihood* getCharacterLikelihoodFunction() const { return characterTreeLikelihood_; }

    /**
     * @brief Return the pointer to the sequence likelihood function
     */
    RNonHomogeneousMixedTreeLikelihood* getSequenceLikelihoodFunction() const { return sequenceTreeLikelihood_; }
 
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
    void setHypothesis(JointLikelihoodFunction::Hypothesis hypothesis);
    
    /**
     * @brief Set the hypothesis data memeber (none / onlyCharacter / onlySequence / both)
     */
    void setOptimizationScope(JointLikelihoodFunction::OptimizationScope optimizationScope) { optimizationScope_ = optimizationScope; }

    /**
     * @brief Incorporates the labels of the nodes into their names for debugging and writing purposes
     * 
     * @param mapping  A tree in which each node has a state propetly assinged to it, corresponding to its binary character state
     */
    void updateStatesInNodesNames(Tree* mapping);

    /**
     * @brief Defines the tree partition into the two sub-models of the sequence model according to a given history
     * 
     * @param history  A tree in which each node has a state propetly assinged to it, corresponding to its binary character state
     */
    void setPartitionByHistory(Tree* history);

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
     * @brief Returns a map of the names and values of the sequence model parameters
     * @bool verbose - indicator weather the scanned parameters should also be reported to stdout
     */
    map<string,double> getModelParameters(bool verbose=true);

    /**
     * @brief Optimizes the parameters of the sequence model
     */
    void optimizeSequenceModel();

    /**
     * @brief Computea and reports the scaling factor of the tree with respect to the sequence data and model
     */
    void reportSequenceScalingFactor();

    /**
     * @brief Reports the joint model paramerers to std
     */
    //void reportResults();

    /**
     * @brief updates the values of the joint model parameters in a map  and maintains the flags characterChanged_ and sequenceChanged_ 
     */
    void updatePreviousParametersValues();

    // functions to override if using reference ot character and sequence models parameters instead of copies
    // all these functions use parameters_ datamember which would be empty in this case

    bool hasParameter(const std::string& name) const { return (characterTreeLikelihood_->hasParameter(name) || sequenceTreeLikelihood_->hasParameter(name)); }
    
    const Parameter& getParameter(const std::string& name) const
    {
      if (characterTreeLikelihood_->hasParameter(name))
      {
        return characterTreeLikelihood_->getParameter(name);  
      }
      else if (sequenceTreeLikelihood_->hasParameter(name))
      {
        return sequenceTreeLikelihood_->getParameter(name); 
      }
      else
      {
        throw ParameterNotFoundException("JointLikelihoodFunction::getParameter", name);
      } 
    }

    const std::shared_ptr<Parameter>& getSharedParameter(const std::string& name) const
    {
      const std::shared_ptr<Parameter>& sharedParameter1  = characterTreeLikelihood_->getSharedParameter(name);
      if (!sharedParameter1)
      {
        const std::shared_ptr<Parameter>& sharedParameter2 = sequenceTreeLikelihood_->getSharedParameter(name);
        return sharedParameter2;
      }
      return sharedParameter1; 
    }

    double getParameterValue(const std::string& name) const
    { 
      double parameterValue;
      string strippedName;
      if (name.find("TwoParameterBinary") != std::string::npos)
      {
        strippedName = name.substr(19);
        parameterValue = characterTreeLikelihood_->getParameter(name).getValue();
      }
      else
      {
        parameterValue = sequenceTreeLikelihood_->getParameter(name).getValue();
      }
      return parameterValue;
    }

    void setAllParametersValues(const ParameterList & parameters)
    {
      characterTreeLikelihood_->setAllParametersValues(parameters);
      sequenceTreeLikelihood_->setAllParametersValues(parameters);
      ParameterList pl;
    }

    void setParameterValue(const std::string& name, double value)
    {
      if (name.find("TwoParameterBinary") != std::string::npos)
      {
        characterTreeLikelihood_->setParameterValue(name, value);
      }
      else
      {
        sequenceTreeLikelihood_->setParameterValue(name, value);
      }
      ParameterList pl;
    }

    void setParametersValues(const ParameterList& parameters)
    { 
      characterTreeLikelihood_->setParametersValues(parameters);
      sequenceTreeLikelihood_->setParametersValues(parameters);
      ParameterList pl;
      fireParameterChanged(pl);
    }

    bool matchParametersValues(const ParameterList& parameters)
    {
      characterTreeLikelihood_->matchParametersValues(parameters);
      sequenceTreeLikelihood_->matchParametersValues(parameters);
      ParameterList pl;
      fireParameterChanged(pl);
      return true;
    }

    size_t getNumberOfParameters() const { return (characterTreeLikelihood_->getParameters().size() + sequenceTreeLikelihood_->getParameters().size()); }
	
	double getSequenceScalingFactor(bool verbose = true);
	
	void scaleSequenceTree(double factor);

  double getLikelihood();

  vector<double> getLikelihoodForEachSite();
};