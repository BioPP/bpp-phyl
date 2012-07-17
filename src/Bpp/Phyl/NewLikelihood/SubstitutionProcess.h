//
// File: SubstitutionProcess.h
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

#ifndef _SUBSTITUTIONPROCESS_H_
#define _SUBSTITUTIONPROCESS_H_

#include "ParametrizableTree.h"
#include "SitePartition.h"
#include "ModelIterator.h"

#include <Bpp/Numeric/ParameterAliasable.h>
#include <Bpp/Numeric/AbstractParameterAliasable.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>

namespace bpp
{

/**
 * @brief This interface describes the substitution process along the tree and sites of the alignment.
 *
 * It main purpose is to provide the necessary calculus for each branch-site-model class combination,
 * such as Markov generator and transition probabilities.
 * These are typically provided by a SubstitutionModel class, applied in various combination along the tree (eg non-homogeneous models)
 * and alignment (eg partition models). The model class is there to account for putative mixture models.
 *
 * As several branches and sites can share the same generator/transition probabilities, calling these values
 * and therefore performing the underlying calculation for each branch-site can result in a very unefficient code.
 * Therefore, "model iterators" are provided, to allow smarter loops over the model structure, by minimizing
 * the amount of computation to be done.
 * Such iterator allow to loop over branches and sites in a clever way, through both direcions, but do not perform
 * any calculations. These are achieved through calls to the corresponding SubstitutionProcess class.
 */
class SubstitutionProcess :
  public virtual ParameterAliasable
{
public:
  virtual SubstitutionProcess* clone() const = 0;

public:
  virtual void registerWithTree(ParametrizableTree* pTree) = 0;
  virtual void registerWithSitePartition(SitePartition* sitePartition) = 0;
  virtual bool isRegistered() const = 0;

  virtual unsigned int getNumberOfClasses() const = 0;

  virtual const Matrix<double>& getTransitionProbabilities(int nodeId, unsigned int siteIndex, unsigned int classIndex) const = 0;
  virtual const Matrix<double>& getGenerator(int nodeId, unsigned int siteIndex, unsigned int classIndex) const = 0;

  virtual ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const = 0;
  virtual ConstSiteModelIterator* getNewSiteModelIterator(unsigned int siteIndex) const = 0;


  /**
   * @brief Tell if the transition probabilities have changed after the last call to setParameters().
   * @return True if trnasition probabilities have changed.
   */
  virtual bool transitionProbabilitiesHaveChanged() const = 0;
};


class AbstractSubstitutionProcess :
  public virtual SubstitutionProcess
{
protected:
  ParametrizableTree* pTree_;
  SitePartition* sitePartition_;

protected:
  AbstractSubstitutionProcess() : pTree_(0), sitePartition_(0) {}

  /**
   * When a process is copied, it becomes unregistered,
   * instead of being registered with the same tree and site partition. This will avoid some hidden
   * bugs...
   */
  AbstractSubstitutionProcess(const AbstractSubstitutionProcess& asp) :
    pTree_(0), sitePartition_(0) {}

  AbstractSubstitutionProcess& operator=(const AbstractSubstitutionProcess& asp)
  {
    pTree_ = 0;
    sitePartition_ = 0;
    return *this;
  }

public:
  void registerWithTree(ParametrizableTree* pTree) { pTree_ = pTree; }
  void registerWithSitePartition(SitePartition* sitePartition) { sitePartition_ = sitePartition; }
  bool isRegistered() const { return pTree_ != 0 && sitePartition_ != 0;}
};


class SimpleSubstitutionProcess :
  public AbstractSubstitutionProcess,
  public AbstractParameterAliasable
{
private:
  SubstitutionModel* model_;

public:
  SimpleSubstitutionProcess(SubstitutionModel* model) :
    AbstractParameterAliasable(model->getNamespace()),
    model_(model)
  {
    if (!model) throw Exception("SimpleSubstitutionProcess. A model instance must be provided.");
    // Add parameters:
    addParameters_(model->getParameters());
  }

  SimpleSubstitutionProcess(const SimpleSubstitutionProcess& ssp) :
    AbstractParameterAliasable(ssp),
    model_(ssp.model_->clone())
  {}

  SimpleSubstitutionProcess& operator=(const SimpleSubstitutionProcess& ssp)
  {
    AbstractParameterAliasable::operator=(ssp),
    delete model_;
    model_ = ssp.model_->clone();
    return *this;
  }

public:
  virtual SimpleSubstitutionProcess* clone() const { return new SimpleSubstitutionProcess(*this); }

  virtual unsigned int getNumberOfClasses() const { return 1; }

  virtual const Matrix<double>& getTransitionProbabilities(int nodeId, unsigned int siteIndex, unsigned int classIndex) const {
    double l = pTree_->getBranchLengthParameter(nodeId).getValue();
    return model_->getPij_t(l);
  }

  virtual const Matrix<double>& getGenerator(int nodeId, unsigned int siteIndex, unsigned int classIndex) const {
    return model_->getGenerator();
  }

  ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const {
    return new ConstNoPartitionBranchModelIterator(model_, sitePartition_->getNumberOfPatternsForPartition(0));
  }

  ConstSiteModelIterator* getNewSiteModelIterator(unsigned int siteIndex) const {
    return new ConstHomogeneousSiteModelIterator(*pTree_, model_);
  }

  bool transitionProbabilitiesHaveChanged() const { return true; }
};

class RateAcrossSitesSubstitutionProcess :
  public AbstractSubstitutionProcess,
  public AbstractParameterAliasable
{
private:
  SubstitutionModel* model_;
  DiscreteDistribution* rDist_;

public:
  RateAcrossSitesSubstitutionProcess(SubstitutionModel* model, DiscreteDistribution* rdist) :
    AbstractParameterAliasable(model->getNamespace()),
    model_(model),
    rDist_(rdist)
  {
    if (!model) throw Exception("RateAcrossSitesSubstitutionProcess. A model instance must be provided.");
    if (!rdist) throw Exception("RateAcrossSitesSubstitutionProcess. A rate distribution instance must be provided.");
    // Add parameters:
    addParameters_(model->getParameters());
    addParameters_(rdist->getParameters());
  }

  RateAcrossSitesSubstitutionProcess(const RateAcrossSitesSubstitutionProcess& rassp) :
    AbstractParameterAliasable(rassp),
    model_(rassp.model_->clone()),
    rDist_(rassp.rDist_->clone())
  {}

  RateAcrossSitesSubstitutionProcess& operator=(const RateAcrossSitesSubstitutionProcess& rassp)
  {
    AbstractParameterAliasable::operator=(rassp),
    delete model_;
    model_ = rassp.model_->clone();
    delete rDist_;
    rDist_ = rassp.rDist_->clone();
    return *this;
  }

public:
  virtual RateAcrossSitesSubstitutionProcess* clone() const { return new RateAcrossSitesSubstitutionProcess(*this); }

  virtual unsigned int getNumberOfClasses() const { return rDist_->getNumberOfCategories(); }

  virtual const Matrix<double>& getTransitionProbabilities(int nodeId, unsigned int siteIndex, unsigned int classIndex) const
  {
    double l = pTree_->getBranchLengthParameter(nodeId).getValue();
    double r = rDist_->getCategory(classIndex);
    return model_->getPij_t(l * r);
  }

  virtual const Matrix<double>& getGenerator(int nodeId, unsigned int siteIndex, unsigned int classIndex) const
  {
    return model_->getGenerator();
  }

  ConstBranchModelIterator* getNewBranchModelIterator(int nodeId) const {
    return new ConstNoPartitionBranchModelIterator(model_, sitePartition_->getNumberOfPatternsForPartition(0));
  }

  ConstSiteModelIterator* getNewSiteModelIterator(unsigned int siteIndex) const {
    return new ConstHomogeneousSiteModelIterator(*pTree_, model_);
  }

  // TODO: it actually depend on the distribution used, how it is parameterized. If classes are fixed and parameters after probabilities only, this should return false to save computational time!
  bool transitionProbabilitiesHaveChanged() const { return true; }
};
} // end namespace bpp

#endif // _SUBSTITUTIONPROCESS_H_

