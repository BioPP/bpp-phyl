//
// File: CodonAsynonymousReversibleSubstitutionModel.h
// Created by: Laurent Gueguen
// Created on: Tue Dec 24 11:03:53 2003
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

#ifndef _CODONASYNONYMOUSREVERSIBLESUBSTITUTIONMODEL_H_
#define _CODONASYNONYMOUSREVERSIBLESUBSTITUTIONMODEL_H_

#include "AbstractCodonReversibleSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>
#include <Bpp/Seq/StateProperties/AlphabetIndex2.h>

namespace bpp
{
/**
 * @brief Class for asynonymous substitution models on codons.
 * @author Laurent Guéguen
 *
 * Objects of this class are built from three reversible substitution
 * models of NucleicAlphabets. No model is directly accessible. </p>
 *
 * Only substitutions with one letter changed are accepted. </p>
 *
 *
 * If a distance @f$d@f$ between amino-acids is defined, the ratio between
 * non-synonymous and synonymous substitutions rates is, if the codied
 * amino-acids are @f$x@f$ and @f$y@f$, @f$\beta*\exp(-\alpha.d(x,y))@f$ with
 * non-negative parameter @f$\alpha@f$ and positive parameter @f$\beta@f$.
 *
 * If such a distance is not defined, the ratio between non-synonymous
 * and synonymous substitutions rates is @f$\beta@f$ with positive
 * parameter @f$\beta@f$.
 */

class CodonAsynonymousReversibleSubstitutionModel :
  public AbstractCodonReversibleSubstitutionModel
{
private:
  const GeneticCode* geneticCode_;
  const AlphabetIndex2<double>* pdistance_;

public:
  /**
   * @brief Build a new CodonNeutralReversibleSubstitutionModel object from
   * a pointer to NucleotideSubstitutionModel.
   *
   * @param palph pointer to a GeneticCode
   * @param pmod  pointer to the NucleotideSubstitutionModel to use in the three positions.
   * The instance will then own this substitution model.
   * @param pdist optional pointer to a distance between amino-acids
   */
  CodonAsynonymousReversibleSubstitutionModel(
    const GeneticCode* palph,
    NucleotideSubstitutionModel* pmod,
    const AlphabetIndex2<double>* pdist = 0);

  /**
   * @brief Build a new CodonNeutralReversibleSubstitutionModel object
   * from three pointers to NucleotideSubstitutionModels.
   *
   * @param palph pointer to a GeneticCode
   * @param pmod1, pmod2, pmod3 pointers to the
   * NucleotideSubstitutionModels to use in the three
   * positions. Either all the models are different objects to avoid
   * parameters redondancy, or only the first model is used in every
   * position. The used models are owned by the instance.
   * @param pdist optional pointer to the AlphabetIndex2<double> amino-acids distance object.
   */
  CodonAsynonymousReversibleSubstitutionModel(
    const GeneticCode* palph,
    NucleotideSubstitutionModel* pmod1,
    NucleotideSubstitutionModel* pmod2,
    NucleotideSubstitutionModel* pmod3,
    const AlphabetIndex2<double>* pdist = 0);

  CodonAsynonymousReversibleSubstitutionModel(
    const CodonAsynonymousReversibleSubstitutionModel& cm) :
    AbstractCodonReversibleSubstitutionModel(cm),
    geneticCode_(cm.geneticCode_),
    pdistance_(cm.pdistance_)
  {}

  CodonAsynonymousReversibleSubstitutionModel & operator=(
    const CodonAsynonymousReversibleSubstitutionModel& cm)
  {
    AbstractCodonReversibleSubstitutionModel::operator=(cm);
    geneticCode_ = cm.geneticCode_;
    pdistance_ = cm.pdistance_;
    return *this;
  }

  ~CodonAsynonymousReversibleSubstitutionModel() {}

  CodonAsynonymousReversibleSubstitutionModel* clone() const
  {
    return new CodonAsynonymousReversibleSubstitutionModel(*this);
  }

public:
  void completeMatrices();

  std::string getName() const;

  const GeneticCode* getGeneticCode() const { return geneticCode_; }
};
} // end of namespace bpp.

#endif

