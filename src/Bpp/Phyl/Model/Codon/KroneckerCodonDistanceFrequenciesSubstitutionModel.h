//
// File: KroneckerCodonDistanceFrequenciesSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: mardi 26 juillet 2016, à 23h 43
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

#ifndef _KRONECKER_CODONDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H_
#define _KRONECKER_CODONDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H_

#include "AbstractKroneckerCodonSubstitutionModel.h"
#include "AbstractCodonDistanceSubstitutionModel.h"
#include "AbstractCodonFrequenciesSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Class for non-synonymous substitution models on codons with
 * parameterized equilibrium frequencies and nucleotidic models,
 * with allowed multiple substitutions.
 *
 * Rates of multiple substitutions equal the product of single
 * substitutions involved, before mgmt of selection and removing stop
 * codons.
 *
 * @author Laurent Guéguen
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, ans does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * See description in AbstractKroneckerCodonDistanceSubstitutionModel
 * and AbstractCodonFrequenciesSubstitutionModel class.
 *
 * The additional parameters to
 * AbstractCodonFrequenciesSubstitutionModel are the rates of
 * nonsynonymous over synonymous substitutions.
 *
 * If a distance @f$d@f$ between amino-acids is defined, the
 *  non-synonymous rate is multiplied with, if the coded amino-acids
 *  are @f$x@f$ and @f$y@f$, @f$\beta*\exp(-\alpha.d(x,y))@f$ with
 *  non-negative parameter \c "alpha" and positive parameter \c
 *  "beta".
 *
 * If such a distance is not defined, the non-synonymous substitution
 *  rate is multiplied with @f$\beta@f$ with positive parameter \c
 *  "beta" (ie @f$d=0@f$).
 *
 * If such a distance is not defined, the ratio between non-synonymous
 * and synonymous substitutions rates is @f$\beta@f$ with positive
 * parameter \c "beta".
 */
  
class KroneckerCodonDistanceFrequenciesSubstitutionModel :
    public AbstractKroneckerCodonSubstitutionModel,
    public AbstractCodonDistanceSubstitutionModel,
    public AbstractCodonFrequenciesSubstitutionModel
{
public:
  /**
   * @brief Build a new KroneckerCodonDistanceFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * Any number of any positions can change simultaneously
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod pointer to the NucleotideSubstitutionModel to use in
   *        the three positions. It is owned by the instance.
   * @param pfreq pointer to the FrequenciesSet* equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids
   *        distance object.
   */

  KroneckerCodonDistanceFrequenciesSubstitutionModel(
      const GeneticCode* gCode,
      NucleotideSubstitutionModel* pmod,
      FrequenciesSet* pfreq,
      const AlphabetIndex2* pdist = 0);

  /**
   * @brief Build a new KroneckerCodonDistanceFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod pointer to the NucleotideSubstitutionModel to use in
   *        the three positions. It is owned by the instance.
   * @param vPos a vector of sets of simultaneously changing
   *   positions.
   * @param pfreq pointer to the FrequenciesSet* equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids
   *        distance object.
   */

  KroneckerCodonDistanceFrequenciesSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod,
    const std::vector<std::set< size_t> >& vPos,
    FrequenciesSet* pfreq,
    const AlphabetIndex2* pdist = 0);

  /**
   * @brief Build a new KroneckerCodonDistanceFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * Any number of any positions can change simultaneously
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod1, pmod2, pmod3 are pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid redundant
   *   parameters.  They are owned by the instance.
   * @param pfreq pointer to the FrequenciesSet* equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids
   *   distance object.
   */

  KroneckerCodonDistanceFrequenciesSubstitutionModel(
      const GeneticCode* gCode,
      NucleotideSubstitutionModel* pmod1,
      NucleotideSubstitutionModel* pmod2,
      NucleotideSubstitutionModel* pmod3,
      FrequenciesSet* pfreq,
      const AlphabetIndex2* pdist = 0);

  /**
   * @brief Build a new KroneckerCodonDistanceFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * @param gCode pointer to a GeneticCode
   * @param pmod1, pmod2, pmod3 are pointers to the
   *   NucleotideSubstitutionModel to use in the three positions.
   *   All the models must be different objects to avoid redundant
   *   parameters.  They are owned by the instance.
   * @param vPos a vector of sets of simultaneously changing
   *   positions.
   * @param pfreq pointer to the FrequenciesSet* equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids
   *   distance object.
   */

  KroneckerCodonDistanceFrequenciesSubstitutionModel(
    const GeneticCode* gCode,
    NucleotideSubstitutionModel* pmod1,
    NucleotideSubstitutionModel* pmod2,
    NucleotideSubstitutionModel* pmod3,
    const std::vector<std::set< size_t> >& vPos,
    FrequenciesSet* pfreq,
    const AlphabetIndex2* pdist = 0);

  virtual ~KroneckerCodonDistanceFrequenciesSubstitutionModel() {}

  KroneckerCodonDistanceFrequenciesSubstitutionModel* clone() const
  {
    return new KroneckerCodonDistanceFrequenciesSubstitutionModel(*this);
  }

public:
  void fireParameterChanged(const ParameterList& parameterlist);

  std::string getName() const;

  double getCodonsMulRate(size_t i, size_t j) const;

  void setNamespace(const std::string&);

  void setFreq(std::map<int,double>& frequencies);

  const FrequenciesSet* getFrequenciesSet() const {
    return AbstractCodonFrequenciesSubstitutionModel::getFrequenciesSet();
  }


};

} // end of namespace bpp.

#endif // _KRONECKER_CODONDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H_

