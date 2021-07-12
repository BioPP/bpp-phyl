//
// File: KCM.h
// Created by: Laurent Gueguen
// Created on: mardi 26 juillet 2016, à 16h 46
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

#ifndef _KCM_H_
#define _KCM_H_

#include "../AbstractBiblioSubstitutionModel.h"
#include "KroneckerCodonDistanceSubstitutionModel.h"

namespace bpp
{
/**
 * @brief The general multiple substitution model for codons, from
 * Zaheri & al, 2014.
 *
 * @author Laurent Guéguen
 *
 * This model is built from one or several nucleotide substitution
 * models. It also allows distinct equilibrium frequencies between
 * codons. A multiplicative factor accounts for the selective
 * restraints at the amino acid level, depending on the synonymy of
 * the amino acids.
 *
 *
 * This model includes :
 *
 * - parameters of the nucleotide models, 
 * - @f$\omega@f$ for the ratio of non-synonymous vs synonymous
 * substitution rates,
 * - parameters of the equilibrium frequencies 
 
 * The codon frequencies @f$\pi_j@f$ are either observed or infered.
 *
 * Reference:
 * -  Zaheri, M. and Dib, L. and Salamin, N. A Generalized
 * Mechanistic Codon Model,  Mol Biol Evol. 2014 Sep;31(9):2528-41.
 * doi: 10.1093/molbev/msu196. Epub 2014 Jun 23.
 *
 */
  
  class KCM :
    public AbstractBiblioSubstitutionModel,
    public virtual CodonReversibleSubstitutionModel
  {
  private:
    std::unique_ptr<KroneckerCodonDistanceSubstitutionModel> pmodel_;
    bool oneModel_;
    
  public:
    /**
     * @brief constructor.
     *
     * If onemod, a unique GTR model is used, otherwise three
     * different GTR models are used.
     *
     **/
    
    KCM(const GeneticCode* gc, bool oneModel);

    KCM(const KCM& kcm);

    KCM& operator=(const KCM&);

    virtual ~KCM() {}

    KCM* clone() const { return new KCM(*this); }

  public:
    std::string getName() const { return "KCM"+std::string(oneModel_?"7":"19")+".";}

    const SubstitutionModel& getSubstitutionModel() const { return *pmodel_.get(); }

    const GeneticCode* getGeneticCode() const { return pmodel_->getGeneticCode(); }
  
    double getCodonsMulRate(size_t i, size_t j) const { return pmodel_->getCodonsMulRate(i, j); }

  protected:
    SubstitutionModel& getSubstitutionModel() { return *pmodel_.get(); }

    const std::shared_ptr<FrequencySet> getFrequencySet() const {
      return 0;
    }

    void setFreq(std::map<int, double>& frequencies)
    {
      AbstractBiblioSubstitutionModel::setFreq(frequencies);
    };

  };

} // end of namespace bpp.

#endif  // _KCM_H_

