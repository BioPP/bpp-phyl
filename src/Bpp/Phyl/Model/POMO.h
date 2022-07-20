//
// File: POMO.h
// Authors: Laurent Guéguen
// Date:  jeudi 23 décembre 2021, à 14h 47
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

#ifndef BPP_PHYL_MODEL_POMO_H
#define BPP_PHYL_MODEL_POMO_H


#include "FrequencySet/FrequencySet.h"
#include "AbstractSubstitutionModel.h"

/**
 * @brief POMO-like Model, based on AllelicAlphabet.
 *
 *
 *
 * In the alphabet from which the AllelicAlphabet is built, this model
 * needs a given mutation model with generator (aka SubstitutionModel
 * in Bio++), and allele selection coefficients from a finess vector
 * (aka FrequencySet in Bio++).
 *
 * If we denote @f$F@f$ this fitness, the drift term towards the
 * fixation of allele @f$i@f$ is proportional @f$F_i@f$, and follows
 * the Moran law.
 *
 * See:
 *
 * De Maio, N., Schrempf, D., and Kosiol, C. (2015). PoMo: An allele
 * frequency-based approach for 279 species tree estimation. Systematic
 * Biology, 64(6):1018–1031.
 *
 * Borges R, Szöllősi GJ, Kosiol C. Quantifying GC-Biased Gene Conversion
 * in Great Ape Genomes Using Polymorphism-Aware Models.
 * Genetics. 2019 Aug;212(4):1321-1336.
 *
 * Borges, R. and Kosiol, C. (2020). Consistency and identifiability
 * of the polymorphism-aware phylo-genetic models. Journal of
 * Theoretical Biology, 486:110074.
 */

#include <Bpp/Seq/Alphabet/AllelicAlphabet.h>

namespace bpp
{
  class POMO :
    public AbstractSubstitutionModel
  {
  private:
    uint nbAlleles_;
  
    std::shared_ptr<SubstitutionModel> pmodel_;
    std::shared_ptr<FrequencySet> pfitness_;

  public:
    /**
     * @brief Build a POMO instance
     *
     */

    POMO(const AllelicAlphabet* allAlph,
         std::shared_ptr<SubstitutionModel> pmodel,
         std::shared_ptr<FrequencySet> pfitness);

    POMO(const POMO& model) :
      AbstractParameterAliasable(model),
      AbstractSubstitutionModel(model),
      nbAlleles_(model.nbAlleles_),
      pmodel_(model.pmodel_->clone()),
      pfitness_(model.pfitness_->clone())
    {}

    POMO& operator=(const POMO& model)
    {
      AbstractParameterAliasable::operator=(model);
      nbAlleles_ = model.nbAlleles_;
      pmodel_ = std::shared_ptr<SubstitutionModel>(model.pmodel_->clone());
      pfitness_ = std::shared_ptr<FrequencySet>(model.pfitness_->clone());
      return *this;
    }

    POMO* clone() const
    {
      return new POMO(*this);
    }

    void fireParameterChanged(const ParameterList& parameters);

    void setFreq(std::map<int, double>& frequencies);

    void setNamespace(const std::string& prefix)
    {
      AbstractParameterAliasable::setNamespace(prefix);
      pmodel_->setNamespace(prefix+pmodel_->getNamespace());
      pfitness_->setNamespace(prefix+pfitness_->getNamespace());
    }

    uint getNbAlleles() const
    {
      return nbAlleles_;
    }

    const std::shared_ptr<FrequencySet> getFitness() const
    {
      return pfitness_;
    }

    const std::shared_ptr<SubstitutionModel> getMutationModel() const
    {
      return pmodel_;
    }

    std::string getName() const
    {
      return "POMO";
    }
    
  protected:
    void updateMatrices();
  
  };
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_POMO_H
