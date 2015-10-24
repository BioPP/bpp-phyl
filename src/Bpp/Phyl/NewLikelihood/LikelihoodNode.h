//
// File: LikelihoodNode.h
// Created by: Laurent Guéguen
// Created on: mardi 23 juin 2015, à 00h 26
// From file HomogeneousTreeLikelihood.h
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

#ifndef _LIKELIHOOD_NODE_H_
#define _LIKELIHOOD_NODE_H_

#include "../Tree/Node.h"

#include <Bpp/Numeric/VectorTools.h>

// From the STL:
#include <map>

namespace bpp
{

/**
 * @brief Likelihood data structure for a node.
 * 
 */
  class LikelihoodNode : 
    public Node 
  {
  public:
    LikelihoodNode() :
      Node() {}
    
    LikelihoodNode(const LikelihoodNode& data) :
      Node(data) {}

    LikelihoodNode(const Node& np) :
      Node(np) {}

    LikelihoodNode(int num, std::string st) :
      Node(num, st) 
    {}

    LikelihoodNode& operator=(const LikelihoodNode& data)
    {
      Node::operator=(*this);

      return *this;
    }
    
  public:
    virtual bool usesLog() const = 0;
    
    virtual VVdouble& getLikelihoodArray(unsigned char DX) = 0;
    virtual const VVdouble& getLikelihoodArray(unsigned char DX) const = 0;
    
    virtual void resetLikelihoods(size_t nbSites, size_t nbStates, unsigned char DX) = 0;

    virtual double& operator()(size_t nSite, size_t nState) = 0;

    virtual double operator()(size_t nSite, size_t nState) const = 0;

  };

} //end of namespace bpp.

#endif //_LIKELIHOOD_NODE_H_

