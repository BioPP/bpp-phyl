//
// File: AbstractCodonDistanceSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: jeudi 15 septembre 2011, à 21h 11
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _ABSTRACTCODONDISTANCESUBSTITUTIONMODEL_H_
#define _ABSTRACTCODONDISTANCESUBSTITUTIONMODEL_H_

#include "CodonSubstitutionModel.h"
#include <Bpp/Numeric/AbstractParameterAliasable.h>


// From SeqLib:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/StateProperties/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Abstract class for modelling of non-synonymous/synonymous
 * ratios of substitution rates in codon models.
 *
 * @author Laurent Guéguen
 *
 * If a distance @f$d@f$ between amino-acids is defined, the ratio
 * between non-synonymous and synonymous substitutions rates is, if
 * the coded amino-acids are @f$x@f$ and @f$y@f$,
 * @f$\beta*\exp(-\alpha.d(x,y))@f$ with non-negative parameter
 * @f$\alpha@f$ and positive parameter @f$\beta@f$.
 *
 * If such a distance is not defined, the ratio between non-synonymous
 * and synonymous substitutions rates is @f$\beta@f$ with positive
 * parameter @f$\beta@f$ (ie @f$d=0@f$).
 */

class AbstractCodonDistanceSubstitutionModel :
    virtual public CodonSubstitutionModel,
    virtual public AbstractParameterAliasable
{
private:
  const GeneticCode* geneticCode_;
  const AlphabetIndex2<double>* pdistance_;

  double alpha_, beta_;
public:
  /**
   *@brief Build a new AbstractCodonDistanceSubstitutionModel object from
   * a pointer to NucleotideSubstitutionModel.
   *
   *@param palph pointer to a GeneticCode
   *@param pdist optional pointer to a distance between amino-acids
   *@param prefix the Namespace
   */
  
  AbstractCodonDistanceSubstitutionModel(
    const GeneticCode* palph,
    const AlphabetIndex2<double>* pdist,
    const std::string& prefix);

  
  AbstractCodonDistanceSubstitutionModel(
    const AbstractCodonDistanceSubstitutionModel& model) :
    AbstractParameterAliasable(model.getNamespace()),
    geneticCode_(model.geneticCode_),
    pdistance_(model.pdistance_),
    alpha_(model.alpha_),
    beta_(model.beta_)
  {}

  AbstractCodonDistanceSubstitutionModel & operator=(
    const AbstractCodonDistanceSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    geneticCode_ = model.geneticCode_;
    pdistance_ = model.pdistance_;
    alpha_ = model.alpha_;
    beta_ = model.beta_;
    return *this;
  }

  ~AbstractCodonDistanceSubstitutionModel() {}

public:
  void fireParameterChanged(const ParameterList& parameters);

  const GeneticCode* getGeneticCode() const { return geneticCode_; }

public:
  double getCodonsMulRate(unsigned int, unsigned int) const;


};
} // end of namespace bpp.

#endif

