//
// File: AbstractLikelihoodNode.h
// Authors:
//   Laurent Guéguen
//   Francois Gindraud (2017)
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

#ifndef _ABSTRACT_LIKELIHOOD_NODE_H_
#define _ABSTRACT_LIKELIHOOD_NODE_H_

#include <Bpp/Phyl/NewLikelihood/DataFlow.h>
#include <Bpp/Phyl/NewLikelihood/ForRange.h>

#include "LikelihoodNode.h"
#include "SpeciationComputingNode.h"

// From the STL:
#include <map>

namespace bpp
{
  /**
   * @brief Basic Likelihood data structure for a node.
   *
   * This class is for use with the AbstractLikelihoodTree class.
   *
   * Store all conditional likelihoods:
   * <pre>
   * x[i][s]
   *   |------> Site i
   *      |---> Ancestral state s
   * </pre>
   * We call this the <i>likelihood array</i> for each node.
   * In the same way, we store first and second order derivatives.
   *
   * @see AbstractLikelihoodTree
   */

  class AbstractLikelihoodNode : public LikelihoodNode
  {
  protected:
    /**
     * @brief Likelihoods results.
     */
    CachedValue<VVdouble> nodeLikelihoods_{};
    CachedValue<VVdouble> nodeDLikelihoods_{};
    CachedValue<VVdouble> nodeD2Likelihoods_{};

    /**
     * @brief vector of pattern pointers from this towards sons (used if patterns only!)
     */
    std::vector<const std::vector<size_t>*> vPatt_{};

    /* @brief says if the likelihood arrays store the log-likelihoods
     * (instead of the likelihoods).
     *
     * Note that it is not valable for the DXlikelihoods.
     *
     */
    bool usesLog_{false};

  public:
    AbstractLikelihoodNode() = default;

    AbstractLikelihoodNode(const PhyloNode& np)
      : LikelihoodNode(np)
    {
    }

    AbstractLikelihoodNode(const AbstractLikelihoodNode&) = default;
    AbstractLikelihoodNode& operator=(const AbstractLikelihoodNode&) = default;

    virtual AbstractLikelihoodNode* clone() const override { return new AbstractLikelihoodNode(*this); }

  public:
    bool usesLog() const { return usesLog_; }

    virtual void setUseLog(bool useLog);

    void setUseLogDownward(bool useLog);

    /**
     * @brief Several Likelihood Arrays
     *
     */

    VVdouble& getLikelihoodArray(unsigned char DX) override
    {
      switch (DX)
      {
        case ComputingNode::D0:
          return nodeLikelihoods_.get_mutable();
        case ComputingNode::D1:
          return nodeDLikelihoods_.get_mutable();
        case ComputingNode::D2:
          return nodeD2Likelihoods_.get_mutable();
        default:
          throw Exception("Unknown derivative " + TextTools::toString(DX));
      }
    }

    const VVdouble& getLikelihoodArray(unsigned char DX) const override
    {
      switch (DX)
      {
        case ComputingNode::D0:
          return nodeLikelihoods_.get();
        case ComputingNode::D1:
          return nodeDLikelihoods_.get();
        case ComputingNode::D2:
          return nodeD2Likelihoods_.get();
        default:
          throw Exception("Unknown derivative " + TextTools::toString(DX));
      }
    }

    bool isUp2date(unsigned char DX) const
    {
      switch (DX)
      {
        case ComputingNode::D0:
          return nodeLikelihoods_.is_valid();
        case ComputingNode::D1:
          return nodeDLikelihoods_.is_valid();
        case ComputingNode::D2:
          return nodeD2Likelihoods_.is_valid();
        default:
          return false;
      }
    }

    void update(bool check, unsigned char DX)
    {
      if (!check)
      {
        switch (DX)
        {
          case ComputingNode::D0:
            nodeLikelihoods_.invalidate();
          case ComputingNode::D1:
            nodeDLikelihoods_.invalidate();
          case ComputingNode::D2:
            nodeD2Likelihoods_.invalidate();
        }
      }
      else
      {
        switch (DX)
        {
          case ComputingNode::D0:
            nodeLikelihoods_.make_valid();
            break;
          case ComputingNode::D1:
            nodeDLikelihoods_.make_valid();
            break;
          case ComputingNode::D2:
            nodeD2Likelihoods_.make_valid();
            break;
        }
      }
    }

    void resetLikelihoods(size_t nbSites, size_t nbStates, unsigned char DX)
    {
      auto& array = getLikelihoodArray(DX);
      array.resize(nbSites);

      for (size_t i = 0; i < nbSites; i++)
        array[i].resize(nbStates, 0);

      update(false, DX);
    }

    double& operator()(size_t nSite, size_t nState) { return nodeLikelihoods_.get_mutable()[nSite][nState]; }

    double operator()(size_t nSite, size_t nState) const { return nodeLikelihoods_.get()[nSite][nState]; }

    /*
     * @brief recursively set likelihood arrays down from the node.
     *
     */

    void resetDownwardLikelihoods(size_t nbSites, size_t nbStates, unsigned char DX)
    {
      resetLikelihoods(nbSites, nbStates, DX);

      for (auto i : make_range(getNumberOfSons()))
      {
        static_cast<AbstractLikelihoodNode*>(getSon(i))->resetDownwardLikelihoods(nbSites, nbStates, DX);
      }

      update(false, DX);
    }

    /*
     * @brief recursively set patterns
     *
     */

    void setPatterns(const std::map<int, std::map<int, std::vector<size_t>>>& patterns)
    {
      const auto& patt = patterns.find(getId())->second;

      for (auto i : make_range(getNumberOfSons()))
      {
        static_cast<AbstractLikelihoodNode*>(getSon(i))->setPatterns(patterns);
        vPatt_.push_back(&(patt.find(getSon(i)->getId())->second));
      }
    }

    /**
     * @brief Compute the posterior probabilities for each state and
     * each distinct site.
     *
     * @param vPP The 2-dimension site X state of the a posteriori
     *        probabilities [in, out].
     */

    void getPosteriorProbabilitiesForEachState(VVdouble& vPP) const;
  };

} // end of namespace bpp.

#endif //_ABSTRACT_LIKELIHOOD_NODE_H_
