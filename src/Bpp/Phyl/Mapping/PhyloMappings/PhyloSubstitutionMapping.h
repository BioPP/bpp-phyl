//
// File: PhyloSubstitutionMapping.h
// Created by: Laurent Guéguen
// Created on: mercredi 29 novembre 2017, à 22h 42
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

#ifndef _PHYLO_SUBSTITUTION_MAPPING_H_
#define _PHYLO_SUBSTITUTION_MAPPING_H_

// From bpp-seq:
#include <Bpp/Numeric/ParameterAliasable.h>
#include <Bpp/Numeric/ParametrizableCollection.h>
#include "../SubstitutionRegister.h"
#include "../../Model/SubstitutionModel.h"

#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
  class PhyloSubstitutionMapping :
    public Clonable
  {
  public:
    PhyloSubstitutionMapping() {}
    virtual ~PhyloSubstitutionMapping() {}

    PhyloSubstitutionMapping* clone() const = 0;

  public:


    /**
     * @brief Get the substitution model corresponding to a certain
     * branch and model class.
     *
     * @param branchId The id of the branch.
     * @param classIndex The model class index.
     */

    virtual const TransitionModel* getModel(unsigned int nodeId, size_t classIndex) const = 0;

    /**
     * @brief Get a list of nodes id for which the given model is associated.
     *
     * @param i The index of the model in the set.
     * @return A vector with the ids of the node associated to this model.
     * @throw IndexOutOfBoundsException If the index is not valid.
     */

//    virtual const std::vector<unsigned int> getNodesWithModel(size_t i) const = 0;

    /*
     * @brief Checks and sets the models with given parameters.
     *
     */
    
    virtual bool matchParametersValues(const ParameterList& nullParams) = 0;

    /*
     * @brief Gets the parameters.
     *
     */
    
    virtual const ParameterList& getParameters() const = 0;

    /*
     * @brief compute Normalizations
     * @param nullParams a list of null parameters
     * @param verbose  
     *
     */

    virtual void computeNormalizations(const ParameterList& nullParams, bool verbose = true) = 0;

    // virtual void computeNormalizationsForASite(size_t site, const ParameterList& nullParams, bool verbose = true) = 0;

    /*
     * @brief return if normalizations have been performed.
     *
     */
     
    virtual bool normalizationsPerformed() const = 0;

    /*
     * @brief return if counts have been performed.
     *
     */
     
    virtual bool countsPerformed() const = 0;
    
    /*
     * @brief ComputeCounts
     * @param threshold 
     * @param verbose  
     *
     */

    virtual void computeCounts(double threshold = -1, bool verbose=true) = 0;

    // virtual void computeCountsForASite(size_t site, double threshold = -1, bool verbose=true) = 0;

    /*
     * @brief Return the tree of counts
     *
     */

    virtual ProbabilisticSubstitutionMapping& getCounts() = 0;
    
    virtual const ProbabilisticSubstitutionMapping& getCounts() const = 0;
    
    /*
     * @brief Return the tree of factors
     *
     */
    
    virtual ProbabilisticSubstitutionMapping& getNormalizations() = 0;
    
    virtual const ProbabilisticSubstitutionMapping& getNormalizations() const = 0;

    /*
     * @brief change Distances
     *
     */

    virtual void setDistances(const AlphabetIndex2 & ndist) = 0;
    
  };

      
} //end of namespace bpp.

#endif  //_PHYLO_SUBSTITUTION_MAPPING_H_

