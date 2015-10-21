//
// File: AbstractSubstitutionProcess.h
// Created by: Julien Dutheil
// Created on: Tue Marc 22 21:17 2013
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

#ifndef _ABSTRACTSUBSTITUTIONPROCESS_H_
#define _ABSTRACTSUBSTITUTIONPROCESS_H_

#include "SubstitutionProcess.h"
#include "ComputingTree.h"

//From the STL:
#include <memory>

//From bpp-core:
#include <Bpp/Numeric/AbstractParameterAliasable.h>

namespace bpp
{

/**
 * @brief A partial implementation of the SubstitutionProcess interface.
 *
 * This class handles a pointer toward a ParametrizableTree object, as well
 * as convenient arrays for storing previously computed probabilities.
 */
  class AbstractSubstitutionProcess :
    public virtual SubstitutionProcess,
    public virtual AbstractParameterAliasable
  {
  protected:
    std::auto_ptr<ParametrizableTree> pTree_;

    size_t nbClasses_;

  protected:
    AbstractSubstitutionProcess(ParametrizableTree* tree, size_t nbClasses, const std::string& prefix = "");

    AbstractSubstitutionProcess(const AbstractSubstitutionProcess& asp);

    AbstractSubstitutionProcess& operator=(const AbstractSubstitutionProcess& asp);

  public:

    const TreeTemplate<Node>& getTree() const {return pTree_->getTree(); }
  
    const ParametrizableTree& getParametrizableTree() const { return *pTree_; }

    size_t getNumberOfClasses() const { return nbClasses_; }

  protected:
    size_t getModelIndex_(int nodeId, size_t modelClass) const throw (NodeNotFoundException, IndexOutOfBoundsException);


  public:

    /**
     * @brief get (Non)Derivable INDEPENDENT parameters
     *
     **/

    ParameterList getDerivableParameters() const
    {
      return getBranchLengthParameters(true);
    }
  
    ParameterList getNonDerivableParameters() const;

    /**
     * @brief AbsractParametrizable interface
     *
     **/
  
    void fireParameterChanged(const ParameterList& pl);

    /**
     * @brief Get the transition probabilities corresponding to a certain branch, site pattern, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */
    const Matrix<double>& getTransitionProbabilities(int nodeId, size_t classIndex) const
    {
      return getComputingTree()[classIndex]->getNode(nodeId)->getTransitionProbabilities();
    }
 
    /**
     * @brief Get the first order derivatives of the transition probabilities according to time, corresponding to a certain branch, site pattern, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */
    const Matrix<double>& getTransitionProbabilitiesD1(int nodeId, size_t classIndex) const
    {
      return getComputingTree()[classIndex]->getNode(nodeId)->getTransitionProbabilitiesD1();
    }
 
    /**
     * @brief Get the second order derivatives of the transition probabilities according to time, corresponding to a certain branch, site pattern, and model class.
     *
     * @param nodeId The id of the node.
     * @param classIndex The model class index.
     */
    const Matrix<double>& getTransitionProbabilitiesD2(int nodeId, size_t classIndex) const
    {
      return getComputingTree()[classIndex]->getNode(nodeId)->getTransitionProbabilitiesD2();
    }

  };
  

} // end namespace bpp
 
#endif // _ABSTRACTSUBSTITUTIONPROCESS_H_

