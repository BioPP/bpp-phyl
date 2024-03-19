// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MAPPING_SUBSTITUTIONCOUNT_H
#define BPP_PHYL_MAPPING_SUBSTITUTIONCOUNT_H

#include <Bpp/Numeric/Matrix/Matrix.h>

#include "../Model/SubstitutionModel.h"
#include "CategorySubstitutionRegister.h"
#include "SubstitutionRegister.h"

// From the STL:
#include <vector>

namespace bpp
{
/**
 * @brief The SubstitutionsCount interface.
 *
 * Provide a method to compute the @f$n_{x,y}(t)@f$ function, namely
 * the number of substitutions on a branch of length @f$t@f$, with
 * initial state @f$x@f$ and final state @f$y@f$.
 *
 * The new implementation offers to perform several counts
 * simultaneously, distinguishing between different types of
 * substitutions. Therefore substitution count object takes as input
 * a SubstitutionRegister, which describes all types of
 * substitutions and associate them with an index. All counts can be
 * retrieved in one go as a vector, the type serving as an indice.
 *
 * @author Julien Dutheil
 *
 * See:
 * Dutheil J, Pupko T, Jean-Marie A, Galtier N.
 * A model-based approach for detecting coevolving positions in a molecule.
 * Mol Biol Evol. 2005 Sep;22(9):1919-28. Epub 2005 Jun 8.
 */

class SubstitutionCountInterface :
  public virtual Clonable
{
public:
  SubstitutionCountInterface() {}
  virtual ~SubstitutionCountInterface() {}
  virtual SubstitutionCountInterface* clone() const = 0;

public:
  /**
   * @return Tell if a substitution register has been attached to this class.
   */
  virtual bool hasSubstitutionRegister() const = 0;

  /**
   * @return The SubstitutionRegister object associated to this instance. The register contains the description of the various substitutions types that are mapped.
   */
  virtual std::shared_ptr<const SubstitutionRegisterInterface> getSubstitutionRegister() const = 0;

  /**
   * @param reg The new SubstitutionRegister object to be associated to this instance.
   * The register contains the description of the various substitutions types that are mapped.
   */
  virtual void setSubstitutionRegister(std::shared_ptr<const SubstitutionRegisterInterface> reg) = 0;

  /**
   * @brief Short cut function, equivalent to getSubstitutionRegister().getNumberOfSubstitutionTypes().
   *
   * @return The number of substitution types supported by this instance.
   */
  virtual size_t getNumberOfSubstitutionTypes() const { return getSubstitutionRegister()->getNumberOfSubstitutionTypes(); }

  /**
   * @brief Short cut function, equivalent to getSubstitutionRegister()->getAlphabet().
   *
   * @return The alphabet associated to this substitution count.
   */
  virtual std::shared_ptr<const Alphabet> getAlphabet() const { return getSubstitutionRegister()->getAlphabet(); }

  /**
   * @brief Short cut function, equivalent to getSubstitutionRegister()->getAlphabet()->getSize().
   *
   * @return The number of states in the model/alphabet.
   */
  virtual size_t getNumberOfStates() const { return getSubstitutionRegister()->getAlphabet()->getSize(); }


  /**
   * @brief Get the number of susbstitutions on a branch, given the initial and final states, and the branch length.
   *
   * @param initialState The initial state.
   * @param finalState   The final state.
   * @param length       The length of the branch.
   * @param type         The type of substitution to count.
   * @return The number of substitutions on a branch of specified length and
   * according to initial and final states.
   */
  virtual double getNumberOfSubstitutions(size_t initialState, size_t finalState, double length, size_t type) const = 0;

  /**
   * @brief Get the numbers of susbstitutions on a branch, for each initial and final states, and given the branch length.
   *
   * @param length       The length of the branch.
   * @param type         The type of susbstitution to count.
   * @return A matrix with all numbers of substitutions for each initial and final states.
   */
  virtual std::unique_ptr< Matrix<double>> getAllNumbersOfSubstitutions(double length, size_t type) const = 0;

  /**
   * @brief Stores the numbers of susbstitutions on a branch, for
   * each initial and final states, and given the branch length.
   *
   * @param length       The length of the branch.
   * @param type         The type of susbstitution to count.
   * @param mat          The matrix filled with all numbers of substitutions
   * for each initial and final states.
   */

  virtual void storeAllNumbersOfSubstitutions(double length, size_t type, Eigen::MatrixXd& mat) const = 0;

  /**
   * @brief Get the numbers of susbstitutions on a branch for all types, for an initial and final states, given the branch length.
   *
   * @param initialState The initial state.
   * @param finalState   The final state.
   * @param length       The length of the branch.
   * @return A matrix with all numbers of substitutions for each initial and final states.
   */
  virtual std::vector<double> getNumberOfSubstitutionsPerType(size_t initialState, size_t finalState, double length) const = 0;

  /**
   * @brief Set the substitution model associated with this count, if relevent.
   *
   * @param model The substitution model to use with this count.
   */
  virtual void setSubstitutionModel(std::shared_ptr<const SubstitutionModelInterface> model) = 0;
};

/**
 * @brief Partial implementation of the SubstitutionCount interface.
 */
class AbstractSubstitutionCount :
  public virtual SubstitutionCountInterface
{
protected:
  std::shared_ptr<const SubstitutionRegisterInterface> register_;

public:
  AbstractSubstitutionCount(std::shared_ptr<const SubstitutionRegisterInterface> reg) :
    register_(reg)
  {}

  virtual ~AbstractSubstitutionCount() {}

public:
  bool hasSubstitutionRegister() const { return register_.get() != nullptr; }

  /**
   * @brief attribution of a SubstitutionRegister
   *
   * @param reg pointer to a SubstitutionRegister
   *
   */
  void setSubstitutionRegister(std::shared_ptr<const SubstitutionRegisterInterface> reg)
  {
    register_ = reg;
    substitutionRegisterHasChanged();
  }

  std::shared_ptr<const SubstitutionRegisterInterface> getSubstitutionRegister() const { return register_; }

protected:
  virtual void substitutionRegisterHasChanged() = 0;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MAPPING_SUBSTITUTIONCOUNT_H
