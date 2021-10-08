//
// File: AbstractCodonAARateSubstitutionModel.h
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

#ifndef _ABSTRACTCODON_AARATE_SUBSTITUTIONMODEL_H_
#define _ABSTRACTCODON_AARATE_SUBSTITUTIONMODEL_H_

#include "CodonSubstitutionModel.h"
#include "../Protein/ProteinSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Abstract class for modelling of non-synonymous and
 *  synonymous substitution rates in codon models, given an amino acid
 *  rate matrix (from a shared_ptr model).
 *
 * @author Laurent Guéguen
 *
 * From the generator @f$g@f$ between amino-acids, the non-synonymous
 *  rate is multiplied with, if the coded amino-acids are @f$x@f$ and
 *  @f$y@f$, @f$\beta*g(x,y)@f$ with positive parameter \c "beta".
 *
 * If paramSynRate is true, the synonymous substitution rate is
 *  multiplied with @f$\gamma@f$ (with optional positive parameter \c
 *  "gamma"), else it is multiplied with 1.
 *
 *
 */

class AbstractCodonAARateSubstitutionModel :
  public virtual CoreCodonSubstitutionModel,
  public virtual AbstractParameterAliasable
{
private:
  std::shared_ptr<ProteinSubstitutionModel> pAAmodel_;

  const GeneticCode* pgencode_;

  double beta_;

  double gamma_;

  std::shared_ptr<const StateMap> stateMap_;

public:
  /**
   * @brief Build a new AbstractCodonAARateSubstitutionModel object from
   *  a pointer to NucleotideSubstitutionModel.
   *
   * @param pmodel shared_ptr to an amino_acid generator
   * @param pgencode the genetic code
   * @param prefix the Namespace
   * @param paramSynRate is true iff synonymous rate is parameterised
   *       (default=false).
   */

  AbstractCodonAARateSubstitutionModel(
    std::shared_ptr<ProteinSubstitutionModel> pmodel,
    const GeneticCode* pgencode,
    const std::string& prefix,
    bool paramSynRate = false);


  AbstractCodonAARateSubstitutionModel(const AbstractCodonAARateSubstitutionModel& model) :
    AbstractParameterAliasable(model),
    pAAmodel_(model.pAAmodel_),
    pgencode_(model.pgencode_),
    beta_(model.beta_),
    gamma_(model.gamma_),
    stateMap_(model.stateMap_)
  {}

  AbstractCodonAARateSubstitutionModel& operator=(
    const AbstractCodonAARateSubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    pAAmodel_ = model.pAAmodel_;
    pgencode_ = model.pgencode_;
    beta_ = model.beta_;
    gamma_ = model.gamma_;
    stateMap_ = model.stateMap_;

    return *this;
  }

  AbstractCodonAARateSubstitutionModel* clone() const
  {
    return new AbstractCodonAARateSubstitutionModel(*this);
  }

  virtual ~AbstractCodonAARateSubstitutionModel() {}

public:
  void fireParameterChanged(const ParameterList& parameters);

  double getCodonsMulRate(size_t i, size_t j) const;

  void setNamespace(const std::string& prefix)
  {
    AbstractParameterAliasable::setNamespace(prefix);
    pAAmodel_->setNamespace(prefix + pAAmodel_->getNamespace());
  }

  /*
   * @brief links to a new AA model
   *
   */
  void setAAModel(std::shared_ptr<ProteinSubstitutionModel> model)
  {
    pAAmodel_ = model;
  }

  const std::shared_ptr<ProteinSubstitutionModel>  getAAModel() const
  {
    return pAAmodel_;
  }

  const std::shared_ptr<FrequencySet> getFrequencySet() const
  {
    return 0;
  }

  void setFreq(std::map<int, double>& frequencies){}
};
} // end of namespace bpp.

#endif// _ABSTRACTCODON_AARATE_SUBSTITUTIONMODEL_H_
