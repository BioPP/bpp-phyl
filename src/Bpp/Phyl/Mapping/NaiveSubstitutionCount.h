//
// File: NaiveSubstitutionCount.h
// Created by: Julien Dutheil
// Created on: Wed Apr 5 11:08 2006
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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

#ifndef _NAIVESUBSTITUTIONCOUNT_H_
#define _NAIVESUBSTITUTIONCOUNT_H_

#include "WeightedSubstitutionCount.h"

#include <Bpp/Numeric/Matrix/Matrix.h>

namespace bpp
{

  /**
   * @brief Naive substitution count.
   *
   * This substitution count is defined as follow:
   * - 0 if @f$i = j@f$,
   * - 1 if @f$i \neq j @f$.
   *
   * Reference (for instance):
   * Tufféry P, Darlu P.
   * Exploring a phylogenetic approach for the detection of correlated substitutions in proteins.
   * Mol Biol Evol. 2000 Nov;17(11):1753-9
   *
   * @author Julien Dutheil
   */
  class NaiveSubstitutionCount:
    public AbstractSubstitutionCount,
    public AbstractWeightedSubstitutionCount
  {
  private:
    bool allowSelf_;
    std::vector<int> supportedChars_;

  public:
    /**
     * @brief Build a new simple substitution count.
     *
     * @param model The substitution model for which this substitution count is parametrized.
     *              The model is not used in the calculation, only for specifying the modeled states.
     * @param reg A pointer toward a substitution register object which discribes the type of substitutions to map.
     * @param allowSelf Tells if "self" mutations, from X to X should be counted together with the ones of type X to Y where X and Y are in the same category, if relevent.
     * The default is "no", to be consistent with other types of substitution counts which account for multiple substitutions, in which case it does not make sense to count "X to X".
     * @param weights the weights of the counts
     */
    NaiveSubstitutionCount(const SubstitutionModel* model, SubstitutionRegister* reg, bool allowSelf = false, const AlphabetIndex2* weights = 0) :
      AbstractSubstitutionCount(reg),
      AbstractWeightedSubstitutionCount(weights, true),
      allowSelf_(allowSelf),
      supportedChars_(model->getAlphabetStates()) {}				
		
    virtual ~NaiveSubstitutionCount() {}
	
    NaiveSubstitutionCount* clone() const { return new NaiveSubstitutionCount(*this); }

  public:
    double getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type = 1) const
    {
      if (initialState == finalState && !allowSelf_)
        return 0;
      else {
        int alphabetState1 = supportedChars_[initialState]; 
        int alphabetState2 = supportedChars_[finalState]; 
        return (register_->getType(initialState, finalState) == type ? (weights_ ? weights_->getIndex(alphabetState1, alphabetState2) : 1.) : 0.);
      }
    }

    Matrix<double>* getAllNumbersOfSubstitutions(double length, size_t type = 1) const;
    
    std::vector<double> getNumberOfSubstitutionsForEachType(size_t initialState, size_t finalState, double length) const
    {
      std::vector<double> v(getNumberOfSubstitutionTypes());
      for (size_t t = 1; t <= getNumberOfSubstitutionTypes(); ++t) {
        v[t - 1] = getNumberOfSubstitutions(initialState, finalState, length, t);
      }
      return v;
    }
    
    void setSubstitutionModel(const SubstitutionModel* model)
    {
      supportedChars_ = model->getAlphabetStates();
    }

  private:
    void substitutionRegisterHasChanged() {}
    void weightsHaveChanged() {}

  };

  /**
   * @brief Labelling substitution count.
   *
   * This substitution count return a distinct number for each possible mutation.
   * - 0 if @f$i = j@f$,
   * - @f$a(i,j)@f$ if @f$i \neq j @f$, where 'a' is an index giving a unique value for each combination of i and j.
   */
  class LabelSubstitutionCount:
    public AbstractSubstitutionCount
  {
  private:
    LinearMatrix<double> label_;
    std::vector<int> supportedChars_;
    
  public:
    LabelSubstitutionCount(const SubstitutionModel* model);

    virtual ~LabelSubstitutionCount() {}

    LabelSubstitutionCount* clone() const { return new LabelSubstitutionCount(*this); }
			
  public:
    double getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type = 1) const
    {
      return label_(initialState, finalState);
    }

    Matrix<double>* getAllNumbersOfSubstitutions(double length, size_t type = 1) const
    {
      return dynamic_cast<Matrix<double>*>(label_.clone());
    }
    
    std::vector<double> getNumberOfSubstitutionsForEachType(size_t initialState, size_t finalState, double length) const
    {
      std::vector<double> v(1);
      v[0] = label_(initialState, finalState);
      return v;
    }

    void setSubstitutionModel(const SubstitutionModel* model)
    {
      supportedChars_ = model->getAlphabetStates();
    }

    void setSubstitutionRegister(SubstitutionRegister* reg) throw (Exception) {
      throw Exception("OneJumpSubstitutionCount::setSubstitutionRegister. This SubstitutionsCount only works with a TotalSubstitutionRegister.");
    }

  private:
    void substitutionRegisterHasChanged() {}

  };

} //end of namespace bpp.

#endif // _SIMPLESUBSTITUTIONCOUNT_H_

