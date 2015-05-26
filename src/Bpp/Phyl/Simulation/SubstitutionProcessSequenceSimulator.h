//
// File: SubstitutionProcessSequenceSimulator.h
// Created by: Laurent Guéguen
// Created on: jeudi 9 avril 2015, à 17h 17
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

#ifndef _SUBSTITUTIONPROCESSSEQUENCESIMULATOR_H_
#define _SUBSTITUTIONPROCESSSEQUENCESIMULATOR_H_

#include "DetailedSiteSimulator.h"
#include "SequenceSimulator.h"
#include "../Tree/TreeTemplate.h"
#include "../Tree/NodeTemplate.h"
#include "../Model/SubstitutionModel.h"

#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Site.h>
#include <Bpp/Seq/Container/SiteContainer.h>

// From the STL:
#include <map>
#include <vector>

#include "../NewLikelihood/SubstitutionProcessCollection.h"
#include "../NewLikelihood/SubstitutionProcess.h"
#include "../NewLikelihood/SequenceEvolution.h"

namespace bpp
{
  
  class SimProcessData
  {
  private:
    size_t state;
    std::vector<size_t> states;
    VVVdouble cumpxy;
    const SubstitutionProcess* process_;
    
  public:
    SimProcessData(): state(), states(), cumpxy(), process_(0) {}
    SimProcessData(const SimProcessData& sd): state(sd.state), states(sd.states), cumpxy(), process_(sd.process_) {}
    SimProcessData& operator=(const SimProcessData& sd)
    {
      state  = sd.state;
      states = sd.states;
      cumpxy = sd.cumpxy;
      process_ = sd.process_;
      return *this;
    }

    friend class SimpleSubstitutionProcessSequenceSimulator;
    
  };
  
  typedef NodeTemplate<SimProcessData> SPNode;

/**
 * @brief Site and sequences simulation under a unique substitution process.
 *
 */

  class SimpleSubstitutionProcessSequenceSimulator:
    public DetailedSiteSimulator,
    public virtual SequenceSimulator
  {
  private:
    const SubstitutionProcess* process_;
    const Alphabet*            alphabet_;
    std::vector<int>           supportedStates_;
    const Tree*                templateTree_;
    mutable TreeTemplate<SPNode> tree_;
  
    /**
     * @brief This stores once for all all leaves in a given order.
     * This order will be used during site creation.
     */
    std::vector<SPNode*> leaves_;
  
    std::vector<std::string> seqNames_;

    size_t nbNodes_;
    size_t nbClasses_;
    size_t nbStates_;

    bool continuousRates_;
  
    /**
     * @name Stores intermediate results.
     *
     * @{
     */

  public:    
    SimpleSubstitutionProcessSequenceSimulator(
      const SubstitutionProcess& process) throw (Exception);

    virtual ~SimpleSubstitutionProcessSequenceSimulator()
    {
    }

    SimpleSubstitutionProcessSequenceSimulator(const SimpleSubstitutionProcessSequenceSimulator& nhss) :
      process_        (nhss.process_),
      alphabet_       (nhss.alphabet_),
      supportedStates_(nhss.supportedStates_),
      templateTree_   (nhss.templateTree_),
      tree_           (nhss.tree_),
      leaves_         (nhss.leaves_),
      seqNames_       (nhss.seqNames_),
      nbNodes_        (nhss.nbNodes_),
      nbClasses_      (nhss.nbClasses_),
      nbStates_       (nhss.nbStates_),
      continuousRates_(nhss.continuousRates_)
    {}

    SimpleSubstitutionProcessSequenceSimulator& operator=(const SimpleSubstitutionProcessSequenceSimulator& nhss)
    {
      process_        = nhss.process_;
      alphabet_        = nhss.alphabet_;
      supportedStates_ = nhss.supportedStates_;
      templateTree_    = nhss.templateTree_;
      tree_            = nhss.tree_;
      leaves_          = nhss.leaves_;
      seqNames_        = nhss.seqNames_;
      nbNodes_         = nhss.nbNodes_;
      nbClasses_       = nhss.nbClasses_;
      nbStates_        = nhss.nbStates_;
      continuousRates_ = nhss.continuousRates_;
      return *this;
    }

    SimpleSubstitutionProcessSequenceSimulator* clone() const { return new SimpleSubstitutionProcessSequenceSimulator(*this); }

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
    Site* simulateSite() const;

    Site* simulateSite(size_t ancestralStateIndex) const;

    Site* simulateSite(size_t ancestralStateIndex, double rate) const;
    
    Site* simulateSite(double rate) const;

    std::vector<std::string> getSequencesNames() const { return seqNames_; }
    /** @} */
    
    /**
     * @name The DetailedSiteSimulator interface.
     *
     * @{
     */
    SiteSimulationResult* dSimulateSite() const;
    
    SiteSimulationResult* dSimulateSite(size_t ancestralStateIndex) const;
    
    SiteSimulationResult* dSimulateSite(size_t ancestralStateIndex, double rate) const;
    
    SiteSimulationResult* dSimulateSite(double rate) const;
    
    /** @} */

    /**
     * @name The SequenceSimulator interface
     *
     * @{
     */
    SiteContainer* simulate(size_t numberOfSites) const;
    /** @} */
    
    /**
     * @name SiteSimulator and SequenceSimulator interface
     *
     * @{
     */
    const Alphabet* getAlphabet() const { return alphabet_; }
    /** @} */

    /**
     * @name Functions with rate classes instead of absolute rates.
     *
     * @{
     */

    // //virtual Site* simulateSite(size_t ancestralStateIndex, size_t rateClass) const;

    // virtual SiteSimulationResult* dSimulateSite(size_t ancestralStateIndex, size_t rateClass) const;
    
    /** @} */
  
    /**
     * @brief Get the substitution process associated to this instance.
     *
     * @return The substitution process associated to this instance.
     */
    const SubstitutionProcess* getSubstitutionProcess() const { return process_; }
    
    /**
     * @brief Get the rate distribution associated to this instance.
     *
     * @return The DiscreteDistribution object associated to this instance.
     */

//const DiscreteDistribution* getRateDistribution() const { return rate_; }

    /**
     * @brief Get the tree associated to this instance.
     *
     * @return The Tree object associated to this instance.
     */
    const Tree* getTree() const { return templateTree_; }

    /**
     * @brief Enable the use of continuous rates instead of discrete rates.
     *
     * To work, the DiscreteDistribution object used should implement
     * the randC method.
     *
     * In this case, sampling is done jointly on the process classes
     * and on the continuous rate.
     *
     * @param yn Tell if we should use continuous rates.
     */
    void enableContinuousRates(bool yn) { continuousRates_ = yn; }
  
  protected:
    
    /**
     * @brief Evolve from an initial state along a branch, knowing the evolutionary class.
     *
     * This method is fast since all pijt have been computed in the constructor of the class.
     * This method is used for the implementation of the SiteSimulator interface.
     */
    
    size_t evolve(const SPNode* node, size_t initialStateIndex, size_t rateClass) const;
    
    /**
     * @brief Evolve from an initial state along a branch, knowing the evolutionary class.
     *
     * This method is fast since all pijt have been computed in the constructor of the class.
     * This method is used for the implementation of the SiteSimulator interface.
     */
    
    size_t evolve(const SPNode* node, size_t initialStateIndex, size_t rateClass, double rate) const;
    
     /**
     * @brief The same as the evolve(initialState, rateClass)
     * function, but for several sites at a time.
     *
     * This method is used for the implementation of the SequenceSimulator interface.
     */
    
    void multipleEvolve(
        const SPNode* node,
        const std::vector<size_t>& initialStateIndices,
        const std::vector<size_t>& rateClasses,
        std::vector<size_t>& finalStates) const;
    
    SiteContainer* multipleEvolve(
      const std::vector<size_t>& initialStates,
      const std::vector<size_t>& rateClasses) const;

    void dEvolve(size_t initialState, size_t rateClass, SiteSimulationResult& ssr) const;

    void dEvolve(size_t initialState, size_t rateClass, double rate, SiteSimulationResult& ssr) const;

    /**
     * @name The 'Internal' methods.
     *
     * @{
     */

    /**
     * This method uses the states_ variable for saving ancestral states.
     */
    void evolveInternal(SPNode* node, size_t rateClass) const;

    /**
     * This method uses the states_ variable for saving ancestral states.
     */
    void evolveInternal(SPNode* node, size_t rateClass, double rate) const;

    /**
     * This method uses the multipleStates_ variable for saving ancestral states.
     */
     void multipleEvolveInternal(SPNode* node, const std::vector<size_t>& rateClasses) const;

    /**
     * This method uses the states_ variable for saving ancestral states.
     */

    void dEvolveInternal(SPNode * node, size_t rateClass, double rate, SiteSimulationResult & ssr) const;

    void dEvolveInternal(SPNode * node, size_t rateClass, SiteSimulationResult & ssr) const;
/** @} */

  };


  /**
   * @brief Sequences simulation under position specific substitution process.
   *
   */

  class SubstitutionProcessSequenceSimulator:
    public SequenceSimulator
  {
  protected:
    /**
     * @brief the vector of the process simulators.
     *
     */
    
    std::map<size_t, SimpleSubstitutionProcessSequenceSimulator*> mProcess_;

    /**
     * @brief The vector of the site specific process in mProcess_;
     * is mutable because can be changed for each simulation (for ex
     * in case of HMM).
     */
    
    mutable std::vector<size_t> vMap_;

    /**
     * @brief all processes trees must have at least the same sequence
     * names as the first process of the map.
     *
     */
    
    std::vector<std::string> seqNames_;

    /**
     * @brief correspondance map of seqNames positions of the several trees.
     * Reference is the tree of the first process of the map.
     *
     * mvPosNames[process id][i] is the position in the id_th tree
     * leaves names of the i_th name of seqName_.
     */

    std::map<size_t, std::vector<size_t> > mvPosNames_;
    
    
  public:
    SubstitutionProcessSequenceSimulator(const SequenceEvolution& evol);

    SubstitutionProcessSequenceSimulator(const std::map<size_t, const SubstitutionProcess&>& mSP);

    SubstitutionProcessSequenceSimulator(const SubstitutionProcessCollection& spc);

    SubstitutionProcessSequenceSimulator(const SubstitutionProcessSequenceSimulator&);
    
    SubstitutionProcessSequenceSimulator& operator=(const SubstitutionProcessSequenceSimulator&);

    SubstitutionProcessSequenceSimulator* clone() const { return new SubstitutionProcessSequenceSimulator(*this); }

    ~SubstitutionProcessSequenceSimulator();

    SimpleSubstitutionProcessSequenceSimulator& getSimpleProcessSimulator(size_t pos)
    {
      if (pos>vMap_.size())
        throw BadIntegerException("Out of range position for SubstitutionProcessSequenceSimulator", (int)pos);
      return *mProcess_[vMap_[pos]];
    }

    /**
     * @brief reset the set of processes.
     *
     */
    
    virtual void resetSiteSimulators(size_t numberOfSites) const
    {
    }

    void setMap(std::vector<size_t> vMap);

    SiteContainer* simulate(size_t numberOfSites) const;

    SiteContainer* simulate(const std::vector<double>& rates) const;

    SiteContainer* simulate(const std::vector<size_t>& states) const;

    SiteContainer* simulate(const std::vector<double>& rates, const std::vector<size_t>& states) const;

    const Alphabet* getAlphabet() const;

  };
  
    

} //end of namespace bpp.

#endif //_SUBSTITUTIONPROCESSSEQUENCESIMULATOR_H_

