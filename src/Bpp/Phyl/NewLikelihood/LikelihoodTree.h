//
// File: LikelihoodTree.h
// Created by: Julien Dutheil, Laurent Guéguen
// Created on: mardi 23 juin 2015, à 14h 16
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

#ifndef _LIKELIHOOD_TREE_H_
#define _LIKELIHOOD_TREE_H_

#include "../Tree/Node.h"
#include "../Tree/TreeTemplate.h"
#include "SubstitutionProcess.h"
#include "AbstractLikelihoodNode.h"

//From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/SiteContainer.h>

namespace bpp
{

/**
 * @brief Interface LikelihoodTree data structure.
 *
 * The structure is initiated according to a tree topology, and 
 * data can be retrieved through node ids.
 *
 * @see LikelihoodNode
 */
  
    class LikelihoodTree:
      public virtual Clonable
    {
    public:
      LikelihoodTree() {}
      virtual ~LikelihoodTree() {}
    
#ifndef NO_VIRTUAL_COV
      LikelihoodTree* clone() const = 0;
#endif

    public:
      virtual const Alphabet* getAlphabet() const = 0;
      virtual size_t getRootArrayPosition(size_t site) const = 0;
      virtual std::vector<size_t>& getRootArrayPositions() = 0; 

      // virtual LikelihoodNode& getNodeData(int nodeId, size_t nClass) = 0;
      
      // virtual const LikelihoodNode& getNodeData(int nodeId, size_t nClass) const = 0;

      /**
       * @return The number of non redundant patterns.
       */
      virtual size_t getNumberOfDistinctSites() const = 0;
    
      /**
       * @return The total number of sites.
       */
      virtual size_t getNumberOfSites() const = 0;
    
      /**
       * @return Get the number of states used in the model.
       */
      virtual size_t getNumberOfStates() const = 0;
    
      /**
       * @return The number of classes used in the model.
       */
      virtual size_t getNumberOfClasses() const = 0;

      /**
       * @return The frequency of a given pattern.
       */
      virtual unsigned int getWeight(size_t pos) const = 0;
    
      /**
       * @return Frequencies for each pattern.
       */
      virtual const std::vector<unsigned int>& getWeights() const = 0;


      /**
       * @brief Resize and initialize all likelihood arrays according
       * to the given sizes, with values computed at leaves by the
       * process on the sites.
       */

      virtual void initLikelihoods(const SiteContainer& sites, const SubstitutionProcess& process) = 0;

      /**
       * @brief Resize and initialize all likelihood arrays at a given
       * node according to the given sizes, for a given derivation
       * class.
       *
       */

      virtual void resetLikelihoods(int nodeId, size_t nbSites, size_t nbStates, unsigned char DX) = 0;

      /**
       * @brief sets using log in all likelihood arrays.
       *
       */
    
      virtual void setAllUseLog(bool useLog) = 0;

    };

} //end of namespace bpp.

#endif //_LIKELIHOOD_TREE_H_

