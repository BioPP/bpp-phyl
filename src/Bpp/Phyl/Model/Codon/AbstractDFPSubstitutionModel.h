//
// File: AbstractDFP.h
// Created by: Laurent Gueguen
// Created on: jeudi 29 octobre 2020, à 15h 59
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

#ifndef _ABSTRACT_DFP_SUBSTITUTIONMODEL_H_
#define _ABSTRACT_DFP_SUBSTITUTIONMODEL_H_

#include "../AbstractSubstitutionModel.h"
#include "CodonSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>

namespace bpp
{

/**
 * @brief Class for neutral substitution models on triplets, following
 * the mutation process proposed in Doron-Fagenboim & Pupko, 2006, but
 * without equilibrium frequencies. This model is an extension of
 * Kimura 2-rates substitution model to codons.
 *
 * There are 5 five free parameters:
 *
 * \c "tr" : rate between codons that differ by 1 transition
 * \c "tv" : rate between codons that differ by 1 transversion. This
 *           rate is set to 1, and other parameters are ratio on this one (so
 *           'tr' is similar to 'kappa' in YN98 model).
 * \c "trr" : rate between codons that differ by 2 transversions
 * \c "tvv" : rate between codons that differ by 2 transversions
 * \c "trv" : rate between codons that differ by 1 transition & 1 transversion
 * \c "tsub" : rate between codons that differ by 3 substitutions
 *
 * Reference: Adi Doron-Faigenboim, Tal Pupko, 2007, A Combined
 * Empirical and Mechanistic Codon Model, Molecular Biology and
 * Evolution, Volume 24, Issue 2, Pages 388–397,
 * https://doi.org/10.1093/molbev/msl175
 *
 */
  
  class AbstractDFPSubstitutionModel :
    public virtual CodonSubstitutionModel,
    public AbstractSubstitutionModel
  {
  private:
    const GeneticCode* gCode_;

    double tr_, trr_, tvv_, trv_, tsub_;
  public:

    /**
     *@brief Build a new AbstractDFPSubstitutionModel object
     *
     *@param palph pointer to a CodonAlphabet
     */
  
    AbstractDFPSubstitutionModel(const GeneticCode* gCode,
                                 const std::string& prefix = "AbstractDFP. ");

    AbstractDFPSubstitutionModel(const AbstractDFPSubstitutionModel& mod) :
      AbstractParameterAliasable(mod),
      AbstractSubstitutionModel(mod),
      gCode_(mod.gCode_),
      tr_(mod.tr_), trr_(mod.trr_), tvv_(mod.tvv_), trv_(mod.trv_), tsub_(mod.tsub_)
    {}

    AbstractDFPSubstitutionModel& operator=(const AbstractDFPSubstitutionModel& mod)
    {
      AbstractSubstitutionModel::operator=(mod);
      gCode_=mod.gCode_;
      tr_ = mod.tr_;
      trr_ = mod.trr_;
      tvv_ =  mod.tvv_;
      trv_ = mod.trv_;
      tsub_ = mod.tsub_;

      return *this;
    }

    ~AbstractDFPSubstitutionModel() {};

    AbstractDFPSubstitutionModel* clone() const = 0;

  public:
    const GeneticCode* getGeneticCode() const { return gCode_; }

    void fireParameterChanged(const ParameterList& parameters);

    using BranchModel::getNumberOfStates;
    size_t getNumberOfStates() { return 64;}
    
    /**
     * @brief Method inherited from AbstractSubstitutionModel
     *
     * This method sets the rates to/from stop codons to zero and
     * set the generator given parameters.
     *
     *
     */
    void updateMatrices();

    /*
     * Calls  the multiplication by the specific codon-codon rate.
     *
     */
    
    double getCodonsMulRate(size_t i, size_t j) const;
    
  };

} //end of namespace bpp.

#endif	

