//
// File: SimpleSubstitutionProcess.h
// Created by: Julien Dutheil
// Created on: Tue May 15 13:11 2012
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

#ifndef _SIMPLESUBSTITUTIONPROCESS_H_
#define _SIMPLESUBSTITUTIONPROCESS_H_

#include "SubstitutionProcess.h"

//From bpp-core:
#include <Bpp/Numeric/ParameterAliasable.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

//From the stl:
#include <memory>

namespace bpp
{

/**
 * @brief Space and time homogeneous substitution process, without mixture.
 */
class SimpleSubstitutionProcess :
  public AbstractSubstitutionProcess,
  public AbstractParameterAliasable
{
protected:
  std::auto_ptr<SubstitutionModel> model_;
  /**
   * @brief All transition probabilities, one set per node.
   */
  std::vector< RowMatrix<double> >probabilities_;
  std::vector<bool> computeProbability_;
  /**
   * @brief The hash table is used to index probability matrices and node ids.
   */
  std::map<int, size_t> nodeIndex_;

public:
  SimpleSubstitutionProcess(SubstitutionModel* model, ParametrizableTree* tree) :
    AbstractSubstitutionProcess(tree, new TrivialSitePartition()),
    AbstractParameterAliasable(model->getNamespace()),
    model_(model),
    probabilities_(),
    computeProbability_(),
    nodeIndex_()
  {
    if (!model)
      throw Exception("SimpleSubstitutionProcess. A model instance must be provided.");
    // Add parameters:
    addParameters_(tree->getParameters());  //Branch lengths
    addParameters_(model->getParameters()); //Substitution model
  }

  SimpleSubstitutionProcess(const SimpleSubstitutionProcess& ssp) :
    AbstractParameterAliasable(ssp),
    model_(ssp.model_->clone())
  {}

  SimpleSubstitutionProcess& operator=(const SimpleSubstitutionProcess& ssp)
  {
    AbstractParameterAliasable::operator=(ssp),
    model_.reset(ssp.model_->clone());
    return *this;
  }

public:
  SimpleSubstitutionProcess* clone() const { return new SimpleSubstitutionProcess(*this); }

  unsigned int getNumberOfClasses() const { return 1; }
  
  unsigned int getNumberOfStates() const { return model_->getNumberOfStates(); }

  ParameterList getTransitionProbabilitiesParameters() const {
    return model_->getParameters();
  }
  
  bool isCompatibleWith(const SiteContainer& data) const {
    return data.getAlphabet()->getAlphabetType() == model_->getAlphabet()->getAlphabetType();
  }
  
  const Matrix<double>& getTransitionProbabilities(int nodeId, unsigned int siteIndex, unsigned int classIndex) const
  {
    size_t i = nodeIndex_[nodeId];
    if (!computeProbability_[i]) {
      computeProbability_[i] = true; //From now on this matrix will be updated.
      //The transition matrix was never coputed before. We therefore have to compute it first:
      double l = pTree_->getBranchLengthParameter(nodeId).getValue();
      probabilities_[i] = model_->getPij_t(l);
    }
    return probabilities_[i];
  }

  const Matrix<double>& getGenerator(int nodeId, unsigned int siteIndex, unsigned int classIndex) const {
    return model_->getGenerator();
  }

  const std::vector<double>& getRootFrequencies(unsigned int siteIndex) const {
    return model_->getFrequencies();
  }
  
  double getInitValue(unsigned int i, int state) const throw (BadIntException) {
    return model_->getInitValue(i, state);
  }
  
  double getProbabilityForModel(unsigned int classIndex) const {
    if (classIndex != 0)
      throw IndexOutOfBoundsException("SimpleSubstitutionProcess::getProbabilityForModel.", classIndex, 0, 1);
    return 1;
  }

  ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const {
    return new ConstNoPartitionBranchModelIterator(model_.get(), sitePartition_->getNumberOfPatternsForPartition(0));
  }

  ConstSiteModelIterator* getNewSiteModelIterator(unsigned int siteIndex) const {
    return new ConstHomogeneousSiteModelIterator(*pTree_, model_.get());
  }

  //bool transitionProbabilitiesHaveChanged() const { return true; }
  protected:
    void fireParameterChanged(const ParameterList& pl); //Forward parameters and updates probabilities if needed.
};

//TODO: remove inheritence, has it contains more transition probabilities
/*
class RateAcrossSitesSubstitutionProcess :
  public SimpleSubstitutionProcess
{
private:
  std::auto_ptr<DiscreteDistribution> rDist_;

public:
  RateAcrossSitesSubstitutionProcess(SubstitutionModel* model, DiscreteDistribution* rdist, ParametrizableTree* tree) :
    SimpleSubstitutionProcess(model, tree),
    rDist_(rdist)
  {
    if (!rdist)
      throw Exception("RateAcrossSitesSubstitutionProcess. A rate distribution instance must be provided.");
    // Add parameters:
    addParameters_(rdist->getParameters());
  }

  RateAcrossSitesSubstitutionProcess(const RateAcrossSitesSubstitutionProcess& rassp) :
    SimpleSubstitutionProcess(rassp),
    rDist_(rassp.rDist_->clone())
  {}

  RateAcrossSitesSubstitutionProcess& operator=(const RateAcrossSitesSubstitutionProcess& rassp)
  {
    SimpleSubstitutionProcess::operator=(rassp),
    rDist_.reset(rassp.rDist_->clone());
    return *this;
  }

public:
  RateAcrossSitesSubstitutionProcess* clone() const { return new RateAcrossSitesSubstitutionProcess(*this); }

  unsigned int getNumberOfClasses() const { return rDist_->getNumberOfCategories(); }

  const Matrix<double>& getTransitionProbabilities(int nodeId, unsigned int siteIndex, unsigned int classIndex) const
  {
    double l = pTree_->getBranchLengthParameter(nodeId).getValue();
    double r = rDist_->getCategory(classIndex);
    return model_->getPij_t(l * r);
  }

  double getProbabilityForModel(unsigned int classIndex) const {
    if (classIndex >= rDist_->getNumberOfCategories())
      throw IndexOutOfBoundsException("RateAcrossSitesSubstitutionProcess::getProbabilityForModel.", classIndex, 0, rDist_->getNumberOfCategories());
    return rDist_->getProbability(classIndex);
  }

  // TODO: it actually depend on the distribution used, how it is parameterized. If classes are fixed and parameters after probabilities only, this should return false to save computational time!
  bool transitionProbabilitiesHaveChanged() const { return true; }
};
*/

} // end namespace bpp

#endif // _SUBSTITUTIONPROCESS_H_

