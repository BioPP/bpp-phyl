//
// File: SimpleSubstitutionProcessSequenceSimulator.h
// Created by: Laurent Guéguen
// Created on: dimanche 24 mai 2020, à 07h 30
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

#ifndef _SIMPLE_SUBSTITUTION_PROCESS_SEQUENCE_SIMULATOR_H_
#define _SIMPLE_SUBSTITUTION_PROCESS_SEQUENCE_SIMULATOR_H_

#include "New_DetailedSiteSimulator.h"
#include "SequenceSimulator.h"
#include "../NewLikelihood/ParametrizablePhyloTree.h"
#include "../Model/SubstitutionModel.h"
#include <Bpp/Phyl/NewLikelihood/ProcessComputationTree.h>


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
#include "../Tree/AwareNode.h"
#include "../Tree/PhyloBranch.h"

namespace bpp
{

  class SimProcessNode:
    public ProcessComputationNode
  {
  private:
    
    // states during simulation
    size_t state_;

    // probabilities to choose in case of mixture node
    Vdouble cumProb_;

    // Sons in case of mixture node
    std::vector<std::shared_ptr<SimProcessNode>> sons_;
    
  public:
    SimProcessNode(const ProcessComputationNode& pcn):
      ProcessComputationNode(pcn), state_() {}

    friend class SimpleSubstitutionProcessSequenceSimulator;
  };
    
  class SimProcessEdge :
    public ProcessComputationEdge
  {
  private:
    // Cumulative pxy for all rates
    VVVdouble cumpxy_;
    
  public:
    SimProcessEdge(const ProcessComputationEdge& pce):
      ProcessComputationEdge(pce), cumpxy_() {}

    friend class SimpleSubstitutionProcessSequenceSimulator;
  };

  typedef AssociationTreeGlobalGraphObserver<SimProcessNode, SimProcessEdge>  SPTree;

/**
 * @brief Site and sequences simulation under a unique substitution process.
 *
 */

  class SimpleSubstitutionProcessSequenceSimulator:
    public New_DetailedSiteSimulator,
    public virtual SequenceSimulator
  {
  private:
    const SubstitutionProcess*     process_;
    const Alphabet*                alphabet_;
    const ParametrizablePhyloTree* phyloTree_;
    SPTree        tree_;

    
    /**
     * @brief Vector of indexes of sequenced output species
     *
     */

    std::vector<size_t> seqIndexes_;
    
    /**
     * @brief Vector of names of sequenced output species
     *
     */
    
    std::vector<std::string> seqNames_; 
    
    /**
     * @brief Map between species Indexes & used nodes, may change at
     * each simulation.
     *
     */

    mutable std::map<size_t, std::shared_ptr<SimProcessNode>> speciesNodes_;
    
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
    SimpleSubstitutionProcessSequenceSimulator(
      const SubstitutionProcess& process);

    virtual ~SimpleSubstitutionProcessSequenceSimulator()
    {
    }

    SimpleSubstitutionProcessSequenceSimulator(const SimpleSubstitutionProcessSequenceSimulator& nhss) :
      process_        (nhss.process_),
      alphabet_       (nhss.alphabet_),
      phyloTree_      (nhss.phyloTree_),
      tree_           (nhss.tree_),
      seqIndexes_     (nhss.seqIndexes_),
      seqNames_       (nhss.seqNames_),
      speciesNodes_  (nhss.speciesNodes_),
      nbNodes_        (nhss.nbNodes_),
      nbClasses_      (nhss.nbClasses_),
      nbStates_       (nhss.nbStates_),
      continuousRates_(nhss.continuousRates_),
      outputInternalSequences_(nhss.outputInternalSequences_)
    {}

    SimpleSubstitutionProcessSequenceSimulator& operator=(const SimpleSubstitutionProcessSequenceSimulator& nhss)
    {
      process_        = nhss.process_;
      alphabet_        = nhss.alphabet_;
      phyloTree_       = nhss.phyloTree_;
      tree_            = nhss.tree_;
      seqIndexes_      = nhss.seqIndexes_;
      seqNames_        = nhss.seqNames_;
      speciesNodes_   = nhss.speciesNodes_;
      nbNodes_         = nhss.nbNodes_;
      nbClasses_       = nhss.nbClasses_;
      nbStates_        = nhss.nbStates_;
      continuousRates_ = nhss.continuousRates_;
      outputInternalSequences_ = nhss.outputInternalSequences_;
      
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

    Site* simulateSite(double rate) const;
    
    Site* simulateSite(size_t ancestralStateIndex, double rate) const;

    std::vector<std::string> getSequencesNames() const { return seqNames_; }
    /** @} */
    
    /**
     * @name The DetailedSiteSimulator interface.
     *
     * @{
     */
    New_SiteSimulationResult* dSimulateSite() const;

    New_SiteSimulationResult* dSimulateSite(size_t ancestralStateIndex) const;

    New_SiteSimulationResult* dSimulateSite(double rate) const;
    
    New_SiteSimulationResult* dSimulateSite(size_t ancestralStateIndex, double rate) const;

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
     * @brief Get the substitution process associated to this instance.
     *
     * @return The substitution process associated to this instance.
     */
    const SubstitutionProcess* getSubstitutionProcess() const { return process_; }
    
    /**
     * @brief Get the tree associated to this instance.
     *
     * @return The Tree object associated to this instance.
     */
    const ParametrizablePhyloTree* getTree() const { return phyloTree_; }

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

    /**
     * @brief Sets whether we will output the internal sequences or not.
     *
     *
     * @param yn Tell if we should output internal sequences.
     */
    void outputInternalSequences(bool yn) ;
    
  protected:

    /**
     * @name The 'Internal' methods.
     *
     * @{
     */

    /**
     * This method uses the states_ variable for saving ancestral states.
     */
    void evolveInternal(std::shared_ptr<SimProcessNode> node, size_t rateClass, New_SiteSimulationResult * ssr = 0) const;

    /**
     * This method uses the states_ variable for saving ancestral states.
     */
    void evolveInternal(std::shared_ptr<SimProcessNode> node, double rate, New_SiteSimulationResult * ssr = 0) const;

     /** @} */

  };


  

} //end of namespace bpp.

#endif //_SIMPLESUBSTITUTIONPROCESSSEQUENCESIMULATOR_H_

