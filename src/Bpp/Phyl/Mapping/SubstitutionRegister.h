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
 * Substitution registers are simple classes that define categories of substitutions, and assign them an index.
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
  public virtual SubstitutionRegister
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
 * @brief Gather states into defined categories, and count the changes between categories.
 *
 * Optionally allows for within categories substitutions.
 */
class CategorySubstitutionRegister:
  public AbstractSubstitutionRegister
{
  protected:
    bool within_;
    unsigned int nbCategories_;
    mutable std::map<int, unsigned int> categories_;
    std::vector< std::vector<unsigned int> > index_;

  public:
    /**
     * @brief Build a new substitution register with categories. This class is mean to be inherited.
     *
     * @param alphabet The input alphabet.
     * @param within Specifies if within categories substitutions should be counted as well.
     */
    CategorySubstitutionRegister(const Alphabet* alphabet, bool within = false):
      AbstractSubstitutionRegister(alphabet),
      within_(within), nbCategories_(0), categories_(), index_()
    {}

  protected:
    template<class T>
    void setCategories(const std::map<int, T>& categories) {
      //First index categories:
      nbCategories_ = 0;
      std::map<T, unsigned int> cats;
      for (typename std::map<int, T>::const_iterator it = categories.begin(); it != categories.end(); ++it) {
        if (cats.find(it->second) == cats.end()) {
          ++nbCategories_;
          cats[it->second] = nbCategories_;
        }
      }

      //Now creates categories:
      std::vector<int> types = alphabet_->getSupportedInts();
      for (size_t i = 0; i < types.size(); ++i) {
        typename std::map<int, T>::const_iterator it = categories.find(types[i]);
        if (it != categories.end()) {
          categories_[types[i]] = cats[it->second];
        } else {
          categories_[types[i]] = 0;
        }
      }
      
      unsigned int count = 1;
      index_.resize(nbCategories_);
      for (size_t i = 0; i < index_.size(); ++i) {
        index_[i].resize(nbCategories_);
        for (size_t j = 0; j < index_.size(); ++j) {
          if (j != i)
            index_[i][j] = count++;
        }
      }
      if (within_) {
        for (size_t i = 0; i < index_.size(); ++i) {
          index_[i][i] = count++;
        }
      }
 
    } 

  public:

    virtual unsigned int getCategory(int type) const {
      if (!alphabet_->isIntInAlphabet(type))
        throw Exception("CategorySubstitutionRegister::getCategory(). Type is not supported by alphabet.");
      return categories_[type];
    }

    virtual bool allowWithin() const { return within_; }
    
    unsigned int getNumberOfCategories() const { return nbCategories_; }

    unsigned int getNumberOfSubstitutionTypes() const { return static_cast<double>(nbCategories_ * (nbCategories_ - 1)) + (within_ ? nbCategories_ : 0); }

    unsigned int getType(int fromState, int toState) const {
      unsigned int fromCat = categories_[fromState];
      unsigned int toCat   = categories_[toState];
      if (fromCat > 0 && toCat > 0)
        return index_[fromCat - 1][toCat - 1];
      else
        return 0;
    }

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
  public CategorySubstitutionRegister
{
  public:
    ExhaustiveSubstitutionRegister(const Alphabet* alphabet, bool within = false):
      CategorySubstitutionRegister(alphabet, within)
    {
      std::map<int, int> categories;
      for (int i = 0; i < static_cast<int>(alphabet->getSize()); ++i) {
        categories[i] = i;
      }
      setCategories<int>(categories);
    }
    
    ExhaustiveSubstitutionRegister* clone() const { return new ExhaustiveSubstitutionRegister(*this); }

};

/**
 * @brief Distinguishes AT<->GC from GC<->AT.
 *
 * This register has two substitution types, mapped as:
 * - 0 not a substitution
 * - 1 a AT->GC substitution
 * - 2 a GC->AT substitution
 */
class GCSubstitutionRegister:
  public CategorySubstitutionRegister
{
  public:
    GCSubstitutionRegister(const NucleicAlphabet* alphabet, bool within = false):
      CategorySubstitutionRegister(alphabet, within)
    {
      std::map<int, int> categories;
      categories[0] = 1;
      categories[1] = 2;
      categories[2] = 2;
      categories[3] = 1;
      setCategories<int>(categories);
    }
    
    GCSubstitutionRegister* clone() const { return new GCSubstitutionRegister(*this); }

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
    bool countMultiple_;

  public:
    DnDsSubstitutionRegister(const GeneticCode* gc, bool countMultiple = false):
      AbstractSubstitutionRegister(gc->getSourceAlphabet()),
      code_(gc),
      countMultiple_(countMultiple)
    {}

    DnDsSubstitutionRegister(const DnDsSubstitutionRegister& reg):
      AbstractSubstitutionRegister(reg),
      code_(reg.code_),
      countMultiple_(reg.countMultiple_)
    {}
 
    DnDsSubstitutionRegister& operator=(const DnDsSubstitutionRegister& reg)
    {
      AbstractSubstitutionRegister::operator=(reg);
      code_ = reg.code_;
      countMultiple_ = reg.countMultiple_;
      return *this;
    }
   
    DnDsSubstitutionRegister* clone() const { return new DnDsSubstitutionRegister(*this); }

  public:
    unsigned int getNumberOfSubstitutionTypes() const { return 2; }

    unsigned int getType(int fromState, int toState) const
    {
      const CodonAlphabet* cAlpha = dynamic_cast<const CodonAlphabet*>(alphabet_);
      if (cAlpha->isStop(fromState) || cAlpha->isStop(toState))
        return 0;
      if (fromState == toState)
        return 0; //nothing happens
      if (!countMultiple_) {
        unsigned int countPos = 0;
        if (cAlpha->getFirstPosition(fromState) != cAlpha->getFirstPosition(toState))
          countPos++;
        if (cAlpha->getSecondPosition(fromState) != cAlpha->getSecondPosition(toState))
          countPos++;
        if (cAlpha->getThirdPosition(fromState) != cAlpha->getThirdPosition(toState))
          countPos++;
        if (countPos > 1)
          return 0;
      }
      return code_->areSynonymous(fromState, toState) ? 1 : 2;
    }
};




} //end of namespace bpp.

#endif //_SUBSTITUTIONREGISTER_H_

