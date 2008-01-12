//
// File: DetailedSiteSimulator.h
// Created by: Julien Dutheil
// Created on: Tue Mar  14 10:51 2006
// from old file DetailedSequenceSimulator.h
// Created on: Wed Aug  24 15:20 2005
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

#ifndef _DETAILEDSITESIMULATOR_H_
#define _DETAILEDSITESIMULATOR_H_

#include "SiteSimulator.h"
#include "TreeTemplate.h"
#include "MutationProcess.h"

// From the STL:
#include <map>
#include <vector>

using namespace std;

namespace bpp
{

/**
 * @brief Data structure to store the result of a DetailedSiteSimulator.
 *
 * This data structure stores each transitional state, and the time when it occured.
 */
class SiteSimulationResult
{
  protected:
    mutable map<int, unsigned int> _indexes;
    unsigned int _currentIndex;
    vector<MutationPath> _paths;
    vector<int> _ancestralStates;
    const Tree * _tree;
    vector<int> _leavesId;
    const Alphabet * _alphabet;
    
  public:
    SiteSimulationResult(const Tree * tree, const Alphabet * alphabet, int ancestralState):
      _currentIndex(0)
    {
      _tree = tree;
      _alphabet = alphabet;
      _indexes[tree->getRootId()] = 0;
      _ancestralStates.push_back(ancestralState);
      _leavesId = tree->getLeavesId();
    }

    virtual ~SiteSimulationResult() {}
  
  public:
    /**
     * @return The alphabet associated to this simulation.
     */
    virtual const Alphabet * getAlphabet() const { return _alphabet; }
    
    virtual void addNode(int nodeId, MutationPath path)
    {
      _currentIndex++;
      _indexes[nodeId] = _currentIndex;
      _paths.push_back(path);
      _ancestralStates.push_back(path.getFinalState());
    }

    virtual int getAncestralState(unsigned int i)    const { return _ancestralStates[i]; }

    virtual int getAncestralState(int nodeId) const { return _ancestralStates[_indexes[nodeId]]; }

    virtual unsigned int getSubstitutionCount(unsigned int i) const { return _paths[i].getNumberOfEvents(); }
    
    virtual unsigned int getSubstitutionCount(int nodeId) const { return _paths[_indexes[nodeId]].getNumberOfEvents(); }
    
    virtual vector<double> getSubstitutionVector() const
    {
      unsigned int n = _paths.size();
      vector<double> counts(n);
      for(unsigned int i = 0; i < n; i++)
        counts[i] = (double)_paths[i].getNumberOfEvents();
      return counts;
    }

    /**
     * @return The states at the leaves.
     */
    virtual vector<int> getFinalStates() const
    {
      unsigned int n = _leavesId.size(); 
      vector<int> states(n);
      for(unsigned int i = 0; i < n; i++)
      {
        states[i] = _ancestralStates[_indexes[_leavesId[i]]];
      }
      return states;
    }

    /**
     * @return The site corresponding to this simulation.
     */
    virtual Site * getSite() const { return new Site(getFinalStates(), _alphabet); }

    /**
     * @return A vector with the leaves names.
     */
    virtual vector<string> getLeaveNames() const
    {
      unsigned int n = _leavesId.size(); 
      vector<string> names(n);
      for(unsigned int i = 0; i < n; i++)
      {
        names[i] = _tree->getNodeName(_leavesId[i]);
      }
      return names;
    }

};

//---------------------------------------------------------------------------

/**
 * @brief Data structure to store the result of a DetailedSiteSimulator.
 *
 * This sructure inherits from the SequenceSimulationResult class, and add support for
 * rate variation across sites.
 */
class RASiteSimulationResult:
  public SiteSimulationResult
{
  protected:
    double _rate;
    
  public:
    RASiteSimulationResult(const Tree* tree, const Alphabet * alphabet, int ancestralState, double rate):
      SiteSimulationResult(tree, alphabet, ancestralState),
      _rate(rate) {}

    virtual ~RASiteSimulationResult() {}
  
  public:
    /**
     * @return The rate of this simulation.
     */
    virtual double getRate() const { return _rate; }
};

//---------------------------------------------------------------------------

/**
 * @brief This interface adds the dSimulate method to the SiteSimulator interface.
 *
 * Instances of this class should be used when a detailed output of the simulation is needed.
 */
class DetailedSiteSimulator:
  public virtual SiteSimulator
{
  public:
    DetailedSiteSimulator() {}
    virtual ~DetailedSiteSimulator() {}
  
  public:
    /**
     * @brief Get a detailed simulation result for one site.
     *
     * @return A SiteSimulationResult object with all ancestral
     * states for all nodes and branches.
     */
    virtual SiteSimulationResult * dSimulate() const = 0;
    virtual SiteSimulationResult * dSimulate(int ancestralState) const = 0;
    virtual SiteSimulationResult * dSimulate(int ancestralState, double rate) const = 0;
    virtual SiteSimulationResult * dSimulate(double rate) const = 0;
    
};

} //end of namespace bpp.

#endif // _DETAILEDSITESIMULATOR_H_

