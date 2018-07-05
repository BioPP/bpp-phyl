//
// File: AbstractCodonCpGSubstitutionModel.h
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

#ifndef _ABSTRACTCODONCPGSUBSTITUTIONMODEL_H_
#define _ABSTRACTCODONCPGSUBSTITUTIONMODEL_H_

#include "CodonSubstitutionModel.h"
#include <Bpp/Numeric/AbstractParameterAliasable.h>


namespace bpp
{
/**
 * @brief Abstract class for modelling of CpG -> CpA or TpG (symetric)
 *  hypermutability substitution rate inside codons. Note that the
 *  neihbouring effects between codons are not considered.
 *
 * @author Laurent Guéguen
 *
 * Substitution rate from C to T (resp. from G to A) is multiplied by
 *  a factor @f$\rho@f$ if C is followed by a G (resp. if G is
 *  following a C).
 *
 * Hypermutability parameter is named \c "rho".
 *
 */

  class AbstractCodonCpGSubstitutionModel :
    public virtual CoreCodonSubstitutionModel,
    public virtual AbstractParameterAliasable
  {
  private:
    double rho_;
    
    std::shared_ptr<StateMap> stateMap_;

  public:
    /**
     * @brief Build a new AbstractCodonCpGSubstitutionModel object from
     *  a pointer to NucleotideSubstitutionModel.
     *
     * @param alphabet the codon alphabet
     * @param prefix the Namespace
     */
    
    AbstractCodonCpGSubstitutionModel(
      const CodonAlphabet& alphabet,
      const std::string& prefix);

    AbstractCodonCpGSubstitutionModel(const AbstractCodonCpGSubstitutionModel& model) :
      AbstractParameterAliasable(model),
      rho_(model.rho_),
      stateMap_(model.stateMap_)
    {}

    AbstractCodonCpGSubstitutionModel& operator=(
      const AbstractCodonCpGSubstitutionModel& model)
    {
      AbstractParameterAliasable::operator=(model);
      rho_ = model.rho_;
      stateMap_ = model.stateMap_;
      return *this;
    }

    AbstractCodonCpGSubstitutionModel* clone() const
    {
      return new AbstractCodonCpGSubstitutionModel(*this);
    }

    virtual ~AbstractCodonCpGSubstitutionModel() {}

  public:
    void fireParameterChanged(const ParameterList& parameters);

    double getCodonsMulRate(size_t i, size_t j) const;

    const FrequenciesSet* getFrequenciesSet() const 
    {
      return 0;
    }

    void setFreq(std::map<int, double>& frequencies){};
  };

} // end of namespace bpp.

#endif

