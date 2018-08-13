//
// File: AbstractCodonDistanceSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: jeudi 15 septembre 2011, à 21h 11
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

#ifndef _ABSTRACTCODON_BGC_SUBSTITUTIONMODEL_H_
#define _ABSTRACTCODON_BGC_SUBSTITUTIONMODEL_H_

#include "CodonSubstitutionModel.h"
#include <Bpp/Numeric/AbstractParameterAliasable.h>


// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Abstract class for modelling of non-synonymous and
 * synonymous substitution rates in codon models, with gBGC.
 *
 * @author Laurent Guéguen
 *
 * The non-synonymous substitution rate is multiplied with
 * @f$\frac{\epsilon B+S}{1-e^{-(\epsilon B+S)}}@f$.
 *
 * The synonymous substitution rate is multiplied with @f$\frac{\epsilon
 * B}{1-e^{-\epsilon B}}@f$.
 *
 * 
 * with positive parameter @f$S@f$ that stands for selection, and real
 * parameter @f$B@f$ for biased gene conversion. In the formula,
 * @f$\epsilon = 1@f$ for AT->GC substitutions, @f$\epsilon = -1@f$
 * for GC->AT substitution, and  @f$\epsilon = 0@f$ otherwise.
 *
 *
 * References:
 * - Galtier N, Duret L, Glémin S, Ranwez V (2009) GC-biased gene
 * conversion promotes the fixation of deleterious amino acid changes
 * in primates, Trends in Genetics, vol. 25(1) pp.1-5.
 *
 */

  class AbstractCodonBGCSubstitutionModel :
    public virtual CoreCodonSubstitutionModel,
    public virtual AbstractParameterAliasable
  {
  private:
    const GeneticCode* pgencode_;
  
    double B_, S_;
    
    std::shared_ptr<StateMap> stateMap_;
    
  public:
    /**
     * @brief Build a new AbstractCodonBGCSubstitutionModel object.
     *
     * @param pgencode the genetic code
     * @param prefix the Namespace
     */
    AbstractCodonBGCSubstitutionModel(
      const GeneticCode* pgencode,
      const std::string& prefix);

    AbstractCodonBGCSubstitutionModel(const AbstractCodonBGCSubstitutionModel& model) :
      AbstractParameterAliasable(model),
      pgencode_(model.pgencode_),
      B_(model.B_),
      S_(model.S_),
      stateMap_(model.stateMap_)
    {}

    AbstractCodonBGCSubstitutionModel& operator=(
      const AbstractCodonBGCSubstitutionModel& model)
    {
      AbstractParameterAliasable::operator=(model);
      pgencode_ = model.pgencode_;
      B_ = model.B_;
      S_ = model.S_;
      stateMap_ = model.stateMap_;

      return *this;
    }

    AbstractCodonBGCSubstitutionModel* clone() const
    {
      return new AbstractCodonBGCSubstitutionModel(*this);
    }
  
    virtual ~AbstractCodonBGCSubstitutionModel() {}

  public:
    void fireParameterChanged(const ParameterList& parameters);

    double getCodonsMulRate(size_t i, size_t j) const;

    const FrequenciesSet* getFrequenciesSet() const 
    {
      return 0;
    }

  };

} // end of namespace bpp.

#endif

