//
// File: DFPDistanceFrequenciesSubstitutionModel.h
// Authors:
//   Laurent Gueguen
// Created: mercredi 4 novembre 2020, ÃÂ  13h 12
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef BPP_PHYL_MODEL_CODON_DFPDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H
#define BPP_PHYL_MODEL_CODON_DFPDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H


#include "AbstractCodonDistanceSubstitutionModel.h"
#include "AbstractCodonFrequenciesSubstitutionModel.h"
#include "AbstractDFPSubstitutionModel.h"

namespace bpp
{
/**
 * @brief Class for non-synonymous substitution models on codons with
 * parameterized equilibrium frequencies and nucleotidic models,
 * with allowed multiple substitutions as parameterized in DFP model.
 *
 * See description in class AbstractDFPSubstitutionModel
 *
 * Reference: Adi Doron-Faigenboim, Tal Pupko, 2007, A Combined
 * Empirical and Mechanistic Codon Model, Molecular Biology and
 * Evolution, Volume 24, Issue 2, Pages 388Ã¢ÂÂ397,
 * https://doi.org/10.1093/molbev/msl175
 *
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, ans does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * See description in AbstractDFPDistanceSubstitutionModel
 * and AbstractCodonFrequenciesSubstitutionModel class.
 *
 * The additional parameters to AbstractDFP and
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

class DFPDistanceFrequenciesSubstitutionModel :
  public AbstractDFPSubstitutionModel,
  public AbstractCodonDistanceSubstitutionModel,
  public AbstractCodonFrequenciesSubstitutionModel
{
public:
  /**
   * @brief Build a new DFPDistanceFrequenciesSubstitutionModel object
   * from three pointers to AbstractSubstitutionModels. NEW
   * AbstractSubstitutionModels are copied from the given ones.
   *
   * Any number of any positions can change simultaneously
   *
   * @param gCode pointer to a GeneticCode
   * @param pfreq pointer to the std::shared_ptr<FrequencySet> equilibrium frequencies
   * @param pdist optional pointer to the AlphabetIndex2 amino-acids
   *        distance object.
   */
  DFPDistanceFrequenciesSubstitutionModel(
    std::shared_ptr<const GeneticCode> gCode,
    std::unique_ptr<CodonFrequencySetInterface> pfreq,
    std::shared_ptr<const AlphabetIndex2> pdist = nullptr);

  virtual ~DFPDistanceFrequenciesSubstitutionModel() {}

  DFPDistanceFrequenciesSubstitutionModel* clone() const override
  {
    return new DFPDistanceFrequenciesSubstitutionModel(*this);
  }

public:
  
  void fireParameterChanged(const ParameterList& parameterlist) override;

  std::string getName() const override;

  double getCodonsMulRate(size_t i, size_t j) const override;

  void setNamespace(const std::string&) override;

  size_t getNumberOfStates() const override
  {
    return 64;
  }

  void setFreq(std::map<int, double>& frequencies) override;

  const CodonFrequencySetInterface& codonFrequencySet() const override
  {
    return AbstractCodonFrequenciesSubstitutionModel::codonFrequencySet();
  }

  bool hasCodonFrequencySet() const override
  {
    return AbstractCodonFrequenciesSubstitutionModel::hasCodonFrequencySet();
  }

};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_CODON_DFPDISTANCEFREQUENCIESSUBSTITUTIONMODEL_H
