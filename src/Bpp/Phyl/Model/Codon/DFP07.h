//
// File: DFP07.h
// Created by: Laurent Gueguen
// Created on: mardi 26 septembre 2017, à 23h 06
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

#ifndef _DFP07_H_
#define _DFP07_H_

#include "../AbstractBiblioMixedTransitionModel.h"
#include "../FrequencySet/CodonFrequencySet.h"
#include "../MixtureOfASubstitutionModel.h"
#include "../Protein/ProteinSubstitutionModel.h"

#include <Bpp/Seq/GeneticCode/GeneticCode.h>

namespace bpp
{

/**
 * @brief Class for non-synonymous substitution models on codons with
 * parameterized equilibrium frequencies and nucleotidic models, with
 * allowed multiple substitutions as parameterized in DFP model, with
 * correction to mimic AA substitution rates from a given protein
 * substitution model.
 *
 * Reference: Adi Doron-Faigenboim, Tal Pupko, 2007, A Combined
 * Empirical and Mechanistic Codon Model, Molecular Biology and
 * Evolution, Volume 24, Issue 2, Pages 388–397,
 * https://doi.org/10.1093/molbev/msl175
 *
 *
 * This class should be used with models which equilibrium
 * distribution is fixed, ans does not depend on the parameters.
 * Otherwise there may be problems of identifiability of the
 * parameters.
 *
 * See description in AbstractDFPSubstitutionModel and
 * AbstractCodonFrequenciesSubstitutionModel and
 * CodonSameAARateSubstitutionModel class.
 *
 * The additional parameters to AbstractDFPSubstitutionModel and
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

  class DFP07:
    public AbstractBiblioMixedTransitionModel,
    virtual public TransitionModel
  {
  protected:

    /*
     * redefined mixed model pointer
     *
     */
    
    const MixtureOfASubstitutionModel*  pmixsubmodel_;
  
    /**
     * @brief indexes of 2 codons states between which the substitution is
     * synonymous, to set a basis to the homogeneization of the rates.
     *
     */
    size_t synfrom_, synto_;
  
  public:
    DFP07(const GeneticCode* gCode, std::shared_ptr<ProteinSubstitutionModel> pAAmodel, std::shared_ptr<CodonFrequencySet> codonFreqs) ;
    
    DFP07(const DFP07& mod2) :
      AbstractBiblioMixedTransitionModel(mod2),
      pmixsubmodel_(),
      synfrom_(mod2.synfrom_),
      synto_(mod2.synto_)
    {
      pmixsubmodel_=dynamic_cast<const MixtureOfASubstitutionModel*>(&getMixedModel());      
    }

    virtual DFP07* clone() const {
      return new DFP07(*this);
    }
    
    DFP07& operator=(const DFP07& mod2)
    {
      const auto& eq = AbstractBiblioMixedTransitionModel::operator=(mod2);
      
      synfrom_ = mod2.synfrom_;
      synto_ = mod2.synto_;
      const auto& mm=eq.getMixedModel();
      
      pmixsubmodel_ = &dynamic_cast<const MixtureOfASubstitutionModel&>(mm);
        
      return *this;
    }
    
  protected:
    void updateMatrices();

  public:
    std::string getName() const { return "DFP07"; }

  };

} //end of namespace bpp.

#endif	//_DFP07_H_

