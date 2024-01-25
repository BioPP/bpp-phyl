//
// File: NonHomogeneousSequenceSimulator.h
// Created by: Julien Dutheil
//             Bastien Boussau
// Created on: Wed Aug  24 15:20 2005
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

#ifndef _NONHOMOGENEOUSSEQUENCESIMULATOR_H_
#define _NONHOMOGENEOUSSEQUENCESIMULATOR_H_

#include "../../Simulation/DetailedSiteSimulator.h"
#include "../../Simulation/SequenceSimulator.h"
#include "../../Tree/TreeTemplate.h"
#include "../../Tree/NodeTemplate.h"
#include "../../Model/SubstitutionModel.h"
#include "../../Likelihood/ParametrizablePhyloTree.h"

#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Site.h>
#include <Bpp/Seq/Container/SiteContainer.h>

// From the STL:
#include <map>
#include <vector>

#include "../Model/SubstitutionModelSet.h"

namespace bpp
{

class SimData
{
  public:
    size_t state;
    std::vector<size_t> states;
    VVVdouble cumpxy;
    std::shared_ptr<const TransitionModelInterface> model;

  public:
    SimData(): state(), states(), cumpxy(), model(0) {}
    SimData(const SimData& sd): state(sd.state), states(sd.states), cumpxy(), model(sd.model) {}
    SimData& operator=(const SimData& sd)
    {
      state  = sd.state;
      states = sd.states;
      cumpxy = sd.cumpxy;
      model  = sd.model;
      return *this;
    }
};

typedef NodeTemplate<SimData> SNode;

/**
 * @brief (Legacy) Site and sequences simulation under non-homogeneous models.
 *
 * Former class maintained for backward compatibility only. 
 * The old class has been ported to the new interfaces, but some methods are left unimplemented.
 * 
 * Rate across sites variation is supported, using a DiscreteDistribution object or by specifying explicitely the rate of the sites to simulate.
 */
class NonHomogeneousSequenceSimulator:
  public DetailedSiteSimulatorInterface,
  public virtual SequenceSimulatorInterface
{
  private:
    std::shared_ptr<const SubstitutionModelSet> modelSet_;
    mutable std::shared_ptr<const Alphabet> alphabet_;
    std::vector<int>            supportedStates_;
    std::shared_ptr<const DiscreteDistributionInterface> rate_;
    std::shared_ptr<const Tree> templateTree_;
    mutable TreeTemplate<SNode> tree_;
    mutable std::shared_ptr<const ParametrizablePhyloTree> phyloTree_;

    /**
     * @brief This stores once for all all leaves in a given order.
     * This order will be used during site creation.
     */
    std::vector<SNode*> leaves_;

    std::vector<std::string> seqNames_;

    size_t nbNodes_;
    size_t nbClasses_;
    size_t nbStates_;

    bool continuousRates_;

    // Should we ouptut internal sequences as well?
    bool outputInternalSequences_;

    /**
     * @name Stores intermediate results.
     *
     * @{
     */

  public:
    NonHomogeneousSequenceSimulator(
      std::shared_ptr<const SubstitutionModelSet> modelSet,
      std::shared_ptr<const DiscreteDistributionInterface> rate,
      std::shared_ptr<const Tree> tree
    );

    NonHomogeneousSequenceSimulator(
      std::shared_ptr<const TransitionModelInterface> model,
      std::shared_ptr<const DiscreteDistributionInterface> rate,
      std::shared_ptr<const Tree> tree
    );

    virtual ~NonHomogeneousSequenceSimulator() {}

    NonHomogeneousSequenceSimulator(const NonHomogeneousSequenceSimulator& nhss) = default;

    NonHomogeneousSequenceSimulator& operator=(const NonHomogeneousSequenceSimulator& nhss) = default;

    NonHomogeneousSequenceSimulator* clone() const override
    { 
      return new NonHomogeneousSequenceSimulator(*this);
    }

  private:
    /**
     * @brief Init all probabilities.
     *
     * Method called by constructors.
     */
    void init();

  public:

    /**
     * @name The SiteSimulator interface
     *
     * @{
     */
    std::unique_ptr<Site> simulateSite() const override;

    std::unique_ptr<Site> simulateSite(size_t ancestralStateIndex) const override;

    std::unique_ptr<Site> simulateSite(size_t ancestralStateIndex, double rate) const override;

    std::unique_ptr<Site> simulateSite(double rate) const override;

    std::vector<std::string> getSequenceNames() const override { return seqNames_; }
  
    /** @} */

    /**
     * @name The DetailedSiteSimulator interface.
     *
     * @{
     */
    std::unique_ptr<SiteSimulationResult> dSimulateSite() const override;

    std::unique_ptr<SiteSimulationResult> dSimulateSite(size_t ancestralStateIndex) const override;

    std::unique_ptr<SiteSimulationResult> dSimulateSite(size_t ancestralStateIndex, double rate) const override;

    std::unique_ptr<SiteSimulationResult> dSimulateSite(double rate) const override;
    /** @} */

    /**
     * @name The SequenceSimulator interface
     *
     * @{
     */
    std::unique_ptr<SiteContainerInterface> simulate(size_t numberOfSites) const override;
    
    const SiteSimulatorInterface& siteSimulator(size_t pos) const override 
    {
      throw Exception("NonHomogeneousSequenceSimulator::siteSimulator: not implemented.");
    }
    /** @} */

    /**
     * @name SiteSimulator and SequenceSimulator interface
     *
     * @{
     */
    std::shared_ptr<const Alphabet> getAlphabet() const override { return alphabet_; }
    
    const Alphabet& alphabet() const override { return *alphabet_; }
    /** @} */

    /**
     * @name Functions with rate classes instead of absolute rates.
     *
     * @{
     */
    virtual std::unique_ptr<Site> simulateSite(size_t ancestralStateIndex, size_t rateClass) const;
    virtual std::unique_ptr<SiteSimulationResult> dSimulateSite(size_t ancestralStateIndex, size_t rateClass) const;
    /** @} */

    /**
     * @brief Get the substitution model associated to this instance.
     *
     * @return The substitution model associated to this instance.
     */
    std::shared_ptr<const SubstitutionModelSet> getSubstitutionModelSet() const { return modelSet_; }

    /**
     * @brief Get the substitution model associated to this instance.
     *
     * @return The substitution model associated to this instance.
     */
    const SubstitutionModelSet& substitutionModelSet() const { return *modelSet_; }

    /**
     * @brief Get the rate distribution associated to this instance.
     *
     * @return The DiscreteDistribution object associated to this instance.
     */
    std::shared_ptr<const DiscreteDistributionInterface> getRateDistribution() const { return rate_; }

    /**
     * @brief Get the rate distribution associated to this instance.
     *
     * @return The DiscreteDistribution object associated to this instance.
     */
    const DiscreteDistributionInterface& rateDistribution() const { return *rate_; }

   /**
     * @brief Get the tree associated to this instance.
     *
     * @return The Tree object associated to this instance.
     */
    std::shared_ptr<const Tree> getTree() const { return templateTree_; }

    /**
     * @brief Get the tree associated to this instance.
     *
     * @return The Tree object associated to this instance.
     */
    const Tree& tree() const { return *templateTree_; }

    /**
     * @brief Enable the use of continuous rates instead of discrete rates.
     *
     * To work, the DiscreteDistribution object used should implement the randC method.
     *
     * @param yn Tell if we should use continuous rates.
     */
    void enableContinuousRates(bool yn) { continuousRates_ = yn; }

    /**
     * @brief Sets whether we will output the internal sequences or not.
     *
     *
     * @param yn Tell if we should output internal sequences.
     */
    void outputInternalSequences(bool yn) override;
    void outputInternalSites(bool yn) override
    {
      throw Exception("NonHomogeneousSequenceSimulator::outputInternalSites: not implemented.");
    }


  protected:

    /**
     * @brief Evolve from an initial state along a branch, knowing the evolutionary rate class.
     *
     * This method is fast since all pijt have been computed in the constructor of the class.
     * This method is used for the implementation of the SiteSimulator interface.
     */
    size_t evolve(const SNode* node, size_t initialStateIndex, size_t rateClass) const;

    /**
     * @brief Evolve from an initial state along a branch, knowing the evolutionary rate.
     *
     * This method is slower than the previous one since exponential terms must be computed.
     * This method is used for the implementation of the SiteSimulator interface.
     */
    size_t evolve(const SNode* node, size_t initialStateIndex, double rate) const;

    /**
     * @brief The same as the evolve(initialState, rateClass) function, but for several sites at a time.
     *
     * This method is used for the implementation of the SequenceSimulator interface.
     */
    void multipleEvolve(
        const SNode* node,
        const std::vector<size_t>& initialStateIndices,
        const std::vector<size_t>& rateClasses,
        std::vector<size_t>& finalStates) const;

    std::unique_ptr<SiteContainerInterface> multipleEvolve(
        const std::vector<size_t>& initialStates,
        const std::vector<size_t>& rateClasses) const;

    void dEvolve(size_t initialState, double rate, RASiteSimulationResult& rassr) const;

    /**
     * @name The 'Internal' methods.
     *
     * @{
     */

    /**
     * This method uses the states_ variable for saving ancestral states.
     */
    void evolveInternal(SNode* node, size_t rateClass) const;
    /**
     * This method uses the states_ variable for saving ancestral states.
     */
    void evolveInternal(SNode* node, double rate) const;
    /**
     * This method uses the multipleStates_ variable for saving ancestral states.
     */
     void multipleEvolveInternal(SNode* node, const std::vector<size_t>& rateClasses) const;

    /**
     * This method uses the states_ variable for saving ancestral states.
     */
    void dEvolveInternal(SNode * node, double rate, RASiteSimulationResult & rassr) const;
    /** @} */

};

} //end of namespace bpp.

#endif //_NONHOMOGENEOUSSEQUENCESIMULATOR_H_
