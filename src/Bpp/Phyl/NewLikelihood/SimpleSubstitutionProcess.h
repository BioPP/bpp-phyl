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

#include "AbstractSubstitutionProcess.h"

//From the stl:
#include <memory>

namespace bpp
{

/**
 * @brief Space and time homogeneous substitution process, without mixture.
 */
class SimpleSubstitutionProcess :
    public AbstractSubstitutionProcess
{
protected:
  std::auto_ptr<SubstitutionModel> model_;

private:
  /**
   * @brief The related Computing Tree
   *
   */

  mutable std::auto_ptr<ComputingTree> computingTree_;

public:
  SimpleSubstitutionProcess(SubstitutionModel* model, ParametrizableTree* tree, bool checkRooted);

  SimpleSubstitutionProcess(const SimpleSubstitutionProcess& ssp);

  SimpleSubstitutionProcess& operator=(const SimpleSubstitutionProcess& ssp);

public:
  SimpleSubstitutionProcess* clone() const { return new SimpleSubstitutionProcess(*this); }

  size_t getNumberOfStates() const { return model_->getNumberOfStates(); }

  bool hasDerivableParameter(const std::string& name) const {
    return !model_->getIndependentParameters().hasParameter(name);
  }
  
  bool isCompatibleWith(const SiteContainer& data) const {
    return data.getAlphabet()->getAlphabetType() == model_->getAlphabet()->getAlphabetType();
  }
 
  const SubstitutionModel& getSubstitutionModel(int nodeId, size_t classIndex) const
  {
    return *model_;
  }

  const DiscreteDistribution* getRateDistribution() const
  {
    return 0;
  }
  
  const Matrix<double>& getGenerator(int nodeId, size_t classIndex) const {
    return model_->getGenerator();
  }

  const FrequenciesSet* getRootFrequenciesSet() const {
    return 0;
  }
  
  ParameterList getSubstitutionModelParameters(bool independent) const
  {
    return (independent?model_->getIndependentParameters():model_->getParameters());
  }

  ParameterList getRateDistributionParameters(bool independent) const {
    return ParameterList();
  }

  ParameterList getRootFrequenciesParameters(bool independent) const {
    return ParameterList();
  }

  ParameterList getBranchLengthParameters(bool independent) const {
    return getParametrizableTree().getParameters();
  }
  
  const std::vector<double>& getRootFrequencies() const {
    return model_->getFrequencies();
  }
  
  double getInitValue(size_t i, int state) const throw (BadIntException) {
    return model_->getInitValue(i, state);
  }
  
  double getProbabilityForModel(size_t classIndex) const {
    if (classIndex != 0)
      throw IndexOutOfBoundsException("SimpleSubstitutionProcess::getProbabilityForModel.", classIndex, 0, 1);
    return 1;
  }

  Vdouble getClassProbabilities() const
  {
    return Vdouble(1,1);
  }

  double getRateForModel(size_t classIndex) const {
    if (classIndex != 0)
      throw IndexOutOfBoundsException("SimpleSubstitutionProcess::getRateForModel.", classIndex, 0, 1);
    return 1;
  }

  const ComputingTree& getComputingTree() const
  {
    return *computingTree_.get();
  }

  ComputingTree& getComputingTree()
  {
    return *computingTree_.get();
  }

  //bool transitionProbabilitiesHaveChanged() const { return true; }
protected:
  void fireParameterChanged(const ParameterList& pl); //Forward parameters and updates probabilities if needed.

  bool modelChangesWithParameter_(size_t i, const ParameterList& pl) const {
    if (pl.getCommonParametersWith(model_->getParameters()).size() > 0)
      return true; 
    ParameterList plbl = pTree_->getBranchLengthParameters(i);
    if (plbl.getCommonParametersWith(pl).size() > 0)
      return true;
    return false;
  }

};

} // end namespace bpp

#endif // _SUBSTITUTIONPROCESS_H_

