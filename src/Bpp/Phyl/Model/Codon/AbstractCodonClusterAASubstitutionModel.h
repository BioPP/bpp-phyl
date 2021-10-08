//
// File: AbstractCodonClusterAASubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: mercredi 22 août 2018, à 14h 11
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

#ifndef _ABSTRACT_CODON_CLUSTER_AA_SUBSTITUTION_MODEL_H_
#define _ABSTRACT_CODON_CLUSTER_AA_SUBSTITUTION_MODEL_H_

#include "CodonSubstitutionModel.h"
#include <Bpp/Numeric/AbstractParameterAliasable.h>


// From bpp-seq:
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Abstract class for modelling of non-synonymous and
 *  synonymous substitution rates in codon models, with AA clustered.
 *
 * @author Laurent Guéguen
 *
 * Non-synonymous rates between amino-acids in the same cluster are
 *  multiplied with @f$\omega_C@f$, and non-synonymous rates between
 *  amino-acids in different clusters are multiplied with @f$\omega_R@f$.
 *
 *  Parameters names are \c "omegaR" and "omegaC".
 *
 *  Clusters can be defined as a vector of assignations, amino-acids
 *  ordered as 3-letters abbreviates ("Ala", "Arg", ...). Default
 *  cluster splits according to polarity and size: "AGPV", "RQEHKWY",
 *  "NDCST", "ILMF", which produces vector:
 *
 * (1,2,3,3,3,2,2,1,2,4,4,2,4,4,1,3,3,2,2,1)
 *
 * References:
 *
 * Sainudiin et al., 2005, Detecting Site-Specific Physicochemical
 *    Selective Pressures: Applications to the Class I HLA of the
 *    Human Major Histocompatibility Complex and the SRK of the Plant
 *    Sporophytic Self-Incompatibility System, Journal of Molecular
 *    Evolution 60(3):315-26
 *
 *
 * Claudia C Weber, Simon Whelan, 2019, Physicochemical Amino Acid
 *    Properties Better Describe Substitution Rates in Large
 *    Populations, Molecular Biology and Evolution, Volume 36, Issue
 *    4, April 2019, Pages 679–690,
 *    https://doi.org/10.1093/molbev/msz003
 *
 */

class AbstractCodonClusterAASubstitutionModel :
  public virtual CoreCodonSubstitutionModel,
  public virtual AbstractParameterAliasable
{
private:
  const GeneticCode* pgencode_;

  double omegaR_, omegaC_;

  std::vector<uint> assign_;

  std::shared_ptr<const StateMap> stateMap_;

public:
  /**
   * @brief Build a new AbstractCodonClusterAASubstitutionModel object.
   *
   * @param pgencode the genetic code
   * @param prefix the Namespace
   * @param assign an paramSynRate is true iff synonymous rate is parametrised
   *       (default categories:   "AGPV", "RQEHKWY", "NDCST", "ILMF")
   */

  AbstractCodonClusterAASubstitutionModel(
    const GeneticCode* pgencode,
    const std::string& prefix,
    const std::vector<uint>& assign = {1, 2, 3, 3, 3, 2, 2, 1, 2, 4, 4, 2, 4, 4, 1, 3, 3, 2, 2, 1});

  AbstractCodonClusterAASubstitutionModel(const AbstractCodonClusterAASubstitutionModel& model) :
    AbstractParameterAliasable(model),
    pgencode_(model.pgencode_),
    omegaR_(model.omegaR_),
    omegaC_(model.omegaC_),
    assign_(model.assign_),
    stateMap_(model.stateMap_)
  {}

  AbstractCodonClusterAASubstitutionModel& operator=(
    const AbstractCodonClusterAASubstitutionModel& model)
  {
    AbstractParameterAliasable::operator=(model);
    pgencode_ = model.pgencode_;
    omegaR_ = model.omegaR_;
    omegaC_ = model.omegaC_;
    assign_ = model.assign_;
    stateMap_ = model.stateMap_;

    return *this;
  }

  AbstractCodonClusterAASubstitutionModel* clone() const
  {
    return new AbstractCodonClusterAASubstitutionModel(*this);
  }

  virtual ~AbstractCodonClusterAASubstitutionModel() {}

public:
  void fireParameterChanged(const ParameterList& parameters);

  double getCodonsMulRate(size_t i, size_t j) const;

  const std::shared_ptr<FrequencySet> getFrequencySet() const
  {
    return 0;
  }

  const std::vector<uint>& getAssign() const
  {
    return assign_;
  }

  void setFreq(std::map<int, double>& frequencies){}
};
} // end of namespace bpp.

#endif
