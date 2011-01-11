//
// File: SubstitutionRegister.h
// Created by: Julien Dutheil
// Created on: Mon Dec 6 16:32 2010
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

#ifndef _SUBSTITUTIONREGISTER_H_
#define _SUBSTITUTIONREGISTER_H_

//From bpp-core:
#include <Bpp/Clonable.h>

//From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

//From the STL:
#include <vector>

namespace bpp
{

/**
 * @brief The SubstitutionRegister interface.
 *
 * Substitution registers are simple classes that defines categories of substitutions, and assign them an index.
 *
 * @author Julien Dutheil
 */
class SubstitutionRegister:
  public virtual Clonable
{
	public:
		SubstitutionRegister() {}
		virtual ~SubstitutionRegister() {}

#ifndef NO_VIRTUAL_COV
    virtual SubstitutionRegister* clone() const = 0;
#endif

	public:
    /**
     * @return The alphabet associated to this instance.
     */
    virtual const Alphabet* getAlphabet() const = 0;

    /**
     * @return The number of substitution types supported by this class.
     */
    virtual unsigned int getNumberOfSubstitutionTypes() const = 0;

    /**
     * @brief Get the substitution type far a given pair of states.
     *
     * @param fromState Initial state (should be a state supported by the specified alphabet).
     * @param toState   Final state (should be a state supported by the specified alphabet).
     * @return The index of the corresponding substitution type, ranging from 0 to 'getNumberOfSubstitutionTypes' + 1,
     * as non-substitution (that is when fromState == toState) will always return 0.
     */
    virtual unsigned int getType(int fromState, int toState) const = 0;
};

class AbstractSubstitutionRegister:
  public SubstitutionRegister
{
  protected:
    const Alphabet* alphabet_;

  public:
    AbstractSubstitutionRegister(const Alphabet* alphabet):
      alphabet_(alphabet)
    {}

    AbstractSubstitutionRegister(const AbstractSubstitutionRegister&asr):
      alphabet_(asr.alphabet_)
    {}

    AbstractSubstitutionRegister& operator=(const AbstractSubstitutionRegister&asr) {
      alphabet_ = asr.alphabet_;
      return *this;
    }

    virtual ~AbstractSubstitutionRegister() {}

  public:
    const Alphabet* getAlphabet() const { return alphabet_; }

};

/**
 * @brief Count all substitutions.
 *
 * This register has only 1 substitution type, mapped as:
 * - 0 not a substitution
 * - 1 a substitution
 */
class TotalSubstitutionRegister:
  public AbstractSubstitutionRegister
{
  public:
    TotalSubstitutionRegister(const Alphabet* alphabet):
      AbstractSubstitutionRegister(alphabet)
    {}

    TotalSubstitutionRegister* clone() const { return new TotalSubstitutionRegister(*this); }

  public:
    unsigned int getNumberOfSubstitutionTypes() const { return 1; }

    unsigned int getType(int fromState, int toState) const {
      return (fromState == toState ? 0 : 1);
    }

};

/**
 * @brief Distinguishes all types of substitutions.
 *
 * This register has only n * (n-1) substitution type, where n is the size of the alphabet, mapped as:
 * - 0 not a substitution
 * - x in [1, n(n-1)] a substitution
 */
class ExhaustiveSubstitutionRegister:
  public AbstractSubstitutionRegister
{
  private:
    unsigned int nbTypes_;
    std::vector< std::vector<unsigned int> > index_;

  public:
    ExhaustiveSubstitutionRegister(const Alphabet* alphabet):
      AbstractSubstitutionRegister(alphabet),
      nbTypes_(0),
      index_(alphabet->getSize())
    {
      nbTypes_ = alphabet->getSize() * (alphabet->getSize() - 1);
      unsigned int count = 1;
      for (size_t i = 0; i < index_.size(); ++i) {
        index_[i].resize(alphabet_->getSize());
        for (size_t j = 0; j < index_.size(); ++j) {
          if (j != i)
            index_[i][j] = count++;
        }
      }
    }
    
    ExhaustiveSubstitutionRegister* clone() const { return new ExhaustiveSubstitutionRegister(*this); }

  public:
    unsigned int getNumberOfSubstitutionTypes() const { return nbTypes_; }

    unsigned int getType(int fromState, int toState) const {
      return (fromState == toState ? 0 : index_[static_cast<int>(fromState)][static_cast<int>(toState)]);
    }


};

/**
 * @brief Distinguishes AT<->GC from GC<->AT.
 *
 * This register has two substitution types, mapped as:
 * - 0 not a substitution
 * - 1 a AT->GC substitution
 * - 2 a GC->AT substitution
 * If "includeZeroCounts" is set to true, then two additionnal types are mapped:
 * - 3 AT->AT
 * - 4 GC->GC
 * so that 0 is never returned.
 */
class GCSubstitutionRegister:
  public AbstractSubstitutionRegister
{
  private:
    bool includeZeroCounts_;

  public:
    GCSubstitutionRegister(const NucleicAlphabet* alphabet, bool includeZeroCounts = false):
      AbstractSubstitutionRegister(alphabet),
      includeZeroCounts_(includeZeroCounts)
    {}
    
    GCSubstitutionRegister* clone() const { return new GCSubstitutionRegister(*this); }

  public:
    unsigned int getNumberOfSubstitutionTypes() const { return includeZeroCounts_ ? 4 : 2; }

    unsigned int getType(int fromState, int toState) const
    {
      if ((fromState == 0 || fromState == 3) && (toState == 1 || toState == 2))
        return 1;
      if ((fromState == 1 || fromState == 2) && (toState == 0 || toState == 3))
        return 2;
      if (includeZeroCounts_) {
        if ((fromState == 0 || fromState == 3) && (toState == 0 || toState == 3))
          return 3;
        if ((fromState == 1 || fromState == 2) && (toState == 1 || toState == 2))
          return 4;
        throw Exception("GCSubstitutionRegister::getType. Unvalid substitution type?");
      } else  {
        return 0;
      }
    }
};

/**
 * @brief Distinguishes transitions from transversions.
 *
 * This register has two substitution types, mapped as:
 * - 0 not a substitution
 * - 1 a transition
 * - 2 a transversion
 */
class TsTvSubstitutionRegister:
  public AbstractSubstitutionRegister
{
  public:
    TsTvSubstitutionRegister(const NucleicAlphabet* alphabet):
      AbstractSubstitutionRegister(alphabet)
    {}
    
    TsTvSubstitutionRegister* clone() const { return new TsTvSubstitutionRegister(*this); }

  public:
    unsigned int getNumberOfSubstitutionTypes() const { return 2; }

    unsigned int getType(int fromState, int toState) const
    {
      if (fromState == toState)
        return 0; //nothing happens
      if ((fromState == 0 && toState == 2)
       || (fromState == 2 && toState == 0)
       || (fromState == 1 && toState == 3)
       || (fromState == 3 && toState == 1))
        return 1; //This is a transition
      return 2; //This is a transversion
    }
};

/**
 * @brief Distinguishes synonymous from non-synonymous substitutions.
 *
 * This register has two substitution types, mapped as:
 * - 0 not a substitution
 * - 1 a synonymous substitution
 * - 2 a non-synonymous substitution
 */
class DnDsSubstitutionRegister:
  public AbstractSubstitutionRegister
{
  private:
    const GeneticCode* code_;

  public:
    DnDsSubstitutionRegister(const GeneticCode* gc):
      AbstractSubstitutionRegister(code_->getSourceAlphabet()),
      code_(gc)
    {}

    DnDsSubstitutionRegister(const DnDsSubstitutionRegister& reg):
      AbstractSubstitutionRegister(reg),
      code_(reg.code_)
    {}
 
    DnDsSubstitutionRegister& operator=(const DnDsSubstitutionRegister& reg)
    {
      AbstractSubstitutionRegister::operator=(reg);
      code_ = reg.code_;
      return *this;
    }
   
    DnDsSubstitutionRegister* clone() const { return new DnDsSubstitutionRegister(*this); }

  public:
    unsigned int getNumberOfSubstitutionTypes() const { return 2; }

    unsigned int getType(int fromState, int toState) const
    {
      if (fromState == toState)
        return 0; //nothing happens
      return code_->areSynonymous(fromState, toState) ? 1 : 2;
    }
};




} //end of namespace bpp.

#endif //_SUBSTITUTIONREGISTER_H_

