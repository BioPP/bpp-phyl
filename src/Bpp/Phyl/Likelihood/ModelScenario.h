//
// File: ModelScenario.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: jeudi 28 novembre 2019, ÃÂ  00h 59
//

/*
  Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef BPP_PHYL_LIKELIHOOD_MODELSCENARIO_H
#define BPP_PHYL_LIKELIHOOD_MODELSCENARIO_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/VectorTools.h>

#include "ModelPath.h"

namespace bpp
{
/**
 * @brief Organization of submodels in mixed substitution models as
 * paths.
 *
 * These sets are defined through an hypergraph, denoted as a
 * scenario, ie a list of hypernodes. HyperNodes are defined in class
 * ModelPath.
 *
 * For example, suppose there are 3 mixed models (M1,M2 and M3),
 * with 2, 3, 4 submodels (S1, S2, ...) each.
 *
 * If the sites are allowed to follow any combination of submodels
 * (12 combinations) the corresponding scenario has only one
 * model path: (<1,2>,<1,2,3>,<1,2,3,4>).
 *
 * The scenario with model paths
 * ((<1>,<1,2>,<1,2>),(<2>,<3>,<3,4>)) means that a site either
 * follows 6 combinations:
 *
 * M1:S1 and ( M2:S1 or M2:S2 ) and ( M3:S1 or M3:S2 ).
 *
 * or
 *
 * M1:S2 and M2:S3 and (M3:S3 or M3:S4).
 *
 *
 * If a mixed model is not represented in any model path, it is
 * considered as non-mixed, which means that his transition
 * probabilities are the mixture of the transition probabilities of
 * its submodels.
 *
 *
 * Dependency of the submodels entails constraints in the
 * probabilities of the submodels, and definition of the model paths
 * must be taken with care for the whole modelling to be possible.
 *
 *
 * In this implementation, for sake of simplification (and for reason
 * of time), each submodel must belong to exactly one given model
 * path, but in theory more complex dependencies are possible.
 *
 *
 * Inside each path model, for each set of submodels, the
 * probabilities of the submodels are conditional probabilities, which
 * means that they sum 1 and their ratio are unchanged.
 *
 * For instance, for scenario ((<1>,<1,2>,<1,2>),(<2>,<3>,<3,4>)),
 * the probabilities of model paths are the probabilities of M1:S1
 * and M1:S2. In the first model path, the probabilities of M2:S1 and
 * M2:S2 are P(M2:S1)/(P(M2:S1)+P(M2:S2)) and
 * P(M2:S2)/(P(M2:S1)+P(M2:S2)).
 *
 * We do not certify that the probability parameters of the mixed
 * models are all useful, and then identifiability problems may be
 * encountered.
 *
 * There is a method ("complete") that creates an additional model
 * path to ensure that all submodels belong to at least an model path.
 *
 * To do (perhaps) : compute the probabilities of the hypernodes from
 * a specific mixed model
 *
 *
 */

class ModelScenario
{
private:
  std::vector<std::shared_ptr<ModelPath> > vModelPaths_;

public:
  ModelScenario() :
    vModelPaths_() {}

  ~ModelScenario(){}

  ModelScenario(std::vector<std::shared_ptr<ModelPath> > vModelPaths) :
    vModelPaths_(vModelPaths)
  {}

  ModelScenario(const ModelScenario& set) :
    vModelPaths_(set.vModelPaths_)
  {}

  ModelScenario& operator=(const ModelScenario& set)
  {
    vModelPaths_ = set.vModelPaths_;
    return *this;
  }

//  ModelScenario* clone() const { return new ModelScenario(*this); }

  /**
   * @brief Resets the list of the ModelPaths
   */
  void clear()
  {
    vModelPaths_.clear();
  }

  /*
   *@brief adds the copy of an ModelPath to the end of the
   * ModelPaths list.
   */
  void addModelPath(std::shared_ptr<ModelPath> hn)
  {
    vModelPaths_.push_back(hn);
  }

  /*
   *@brief change from a MixedTransitionModel to a
   * MixedTransitionModel in all ModelPath.
   *
   */

  void changeModel(std::shared_ptr<MixedTransitionModel> m1,
                   std::shared_ptr<MixedTransitionModel> m2);

  /*
   *@brief If necessary, adds a new ModelPath such that all submodels
   *       of the declared mixture models are at least in an
   *       ModelPath.
   *
   * Returns true iff a new path has been built.
   *
   */

  bool complete();

  /*
   *@brief adds a submodel number to the nMth mixed model of the
   *  nHth ModelPath of the list (default nH=0). Checks if all the
   *  numbers are valid.
   *
   *@param nM number of the mixed model
   *@param vnS number of the submodel
   *@param nH number of the concerned ModelPath (default the last element of
   *     the list)
   */

//  void addToModelPath(size_t nM, const Vint& vnS, int nH = -1);

  size_t getNumberOfModelPaths() const { return vModelPaths_.size(); }

  std::shared_ptr<ModelPath> getModelPath(size_t i) {return vModelPaths_[i]; }

  std::shared_ptr<const ModelPath> getModelPath(size_t i) const {return vModelPaths_[i]; }

  /*
   * @brief return models found in several paths
   *
   */

  std::vector<std::shared_ptr<MixedTransitionModel> > getModels() const;

  /*
   *@brief Checks if all the path (ie hypernodes) are exclusive.
   *
   */

  bool hasExclusivePaths() const;

//  void fireParameterChanged(const ParameterList& parameters);

  /*
   *@brief compute the probabilities in all the ModelPaths
   *
   * The probabilities of the hypernodes are computed from the lead
   * mixed model of the FIRST ModelPath (and not of the others).
   *
   */

  void computeModelPathsProbabilities();

  /*
   *@brief computes the probability of an ModelPath, given the
   *     conditional probabilities of the submodels computed from the
   *     hypernodes of this ModelScenario object. If the
   *     ModelPath does not match the structure of allowed by this
   *     ModelScenario, an Exception is thrown.
   *
   *     The probability of an ModelPath is the product -- on the set
   *     of the mixed models -- of the sums of the conditional
   *     probabilities of the submodels that belong to this hypernode
   *     for each mixed model.
   *
   *@param hn the ModelPath which conditional probability is computed.
   */

  /**
   * @brief string description
   *
   */

  std::string to_string() const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_MODELSCENARIO_H
