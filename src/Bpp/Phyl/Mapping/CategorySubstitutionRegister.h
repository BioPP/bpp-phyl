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


#ifndef _CATEGORY_SUBSTITUTIONREGISTER_H_
#define _CATEGORY_SUBSTITUTIONREGISTER_H_

#include "SubstitutionRegister.h"

namespace bpp
{
/**
 * @brief The CategorySubstitutionRegisters
 *
 * @author Julien Dutheil
 */

  
/**
 * @brief Gather states into defined categories, and count the changes between categories.
 *
 * Optionally allows for within categories substitutions.
 */
  class CategorySubstitutionRegister :
    public AbstractSubstitutionRegister
  {
  protected:
    bool within_;
    size_t nbCategories_;
    mutable std::map<size_t, size_t> categories_;
    std::vector<std::string> categoryNames_;
    std::vector< std::vector<size_t> > index_;
    std::vector< std::vector<size_t> > revIndex_;

    bool stationarity_;
    
  public:
    /**
     * @brief Build a new substitution register with categories. This class is meant to be inherited.
     *
     * @param model The model defining the states.
     * @param within Specifies if within categories substitutions should be counted as well.
     */
    CategorySubstitutionRegister(const SubstitutionModel* model, bool within = false) :
      AbstractSubstitutionRegister(model),
      within_(within),
      nbCategories_(0),
      categories_(),
      categoryNames_(),
      index_(),
      revIndex_(),
      stationarity_()
    {}

  protected:
    template<class T>
    void setAlphabetCategories(const std::map<int, T>& categories)
    {
      //We need to convert alphabet states into model states.
      std::map<size_t, T> modelCategories;
      for (typename std::map<int, T>::const_iterator it = categories.begin(); it != categories.end(); ++it)
      {  
        std::vector<size_t> states = model_->getModelStates(it->first);
        for (size_t i = 0; i < states.size(); ++i) {
          modelCategories[states[i]] = it->second;
        }
      }
      //Then we forward the model categories:
      setModelCategories<T>(modelCategories);
    }

    template<class T>
    void setModelCategories(const std::map<size_t, T>& categories)
    {
      // First index categories:
      nbCategories_ = 0;
      std::map<T, size_t> cats;
      for (typename std::map<size_t, T>::const_iterator it = categories.begin(); it != categories.end(); ++it)
      {
        if (cats.find(it->second) == cats.end())
        {
          ++nbCategories_;
          cats[it->second] = nbCategories_;
        }
      }

      // Now creates categories:
      categories_.clear();
      categoryNames_.resize(nbCategories_);
      for (size_t i = 0; i < model_->getNumberOfStates(); ++i)
      {
        typename std::map<size_t, T>::const_iterator it = categories.find(i);
        if (it != categories.end())
        {
          categories_[i] = cats[it->second];
          categoryNames_[cats[it->second] - 1] += model_->getStateMap().getStateDescription(i);
        }
        else
        {
          categories_[i] = 0;
        }
      }

      size_t count = 1;
      index_.resize(nbCategories_);
      for (size_t i = 0; i < index_.size(); ++i)
      {
        index_[i].resize(nbCategories_);
        for (size_t j = 0; j < index_.size(); ++j)
        {
          if (j != i)
          {
            index_[i][j] = count++;
            std::vector<size_t> pos(2);
            pos[0] = i;
            pos[1] = j;
            revIndex_.push_back(pos);
          }
        }
      }
      if (within_)
      {
        for (size_t i = 0; i < index_.size(); ++i)
        {
          index_[i][i] = count++;
          std::vector<size_t> pos(2);
          pos[0] = i;
          pos[1] = i;
          revIndex_.push_back(pos);
        }
      }
    }

  public:
    virtual size_t getCategory(size_t state) const
    {
      return categories_[state];
    }

    virtual size_t getCategoryFrom(size_t type) const
    {
      if (type <= nbCategories_ * (nbCategories_ - 1))
      {
        return revIndex_[type - 1][0] + 1;
      }
      else
      {
        if (within_)
          return revIndex_[type - 1][0] + 1;
        else
          throw Exception("CategorySubstitutionRegister::getCategoryFrom. Bad substitution type.");
      }
    }

    virtual size_t getCategoryTo(size_t type) const
    {
      if (type <= nbCategories_ * (nbCategories_ - 1))
      {
        return revIndex_[type - 1][1] + 1;
      }
      else
      {
        if (within_)
          return revIndex_[type - 1][1] + 1;
        else
          throw Exception("CategorySubstitutionRegister::getCategoryTo. Bad substitution type.");
      }
    }

    virtual std::string getCategoryName(size_t category) const
    {
      return categoryNames_[category - 1];
    }

    virtual bool allowWithin() const { return within_; }

    bool isStationary() const
    {
      return stationarity_;
    }

    void setStationarity(bool stat)
    {
      stationarity_=stat;
    }
    
    size_t getNumberOfCategories() const { return nbCategories_; }

    size_t getNumberOfSubstitutionTypes() const { return nbCategories_ * (nbCategories_ - 1) + (within_ ? nbCategories_ : 0); }

    virtual size_t getType(size_t fromState, size_t toState) const
    {
      size_t fromCat = categories_[fromState];
      size_t toCat   = categories_[toState];
      if (fromCat > 0 && toCat > 0)
        return index_[fromCat - 1][toCat - 1];
      else
        return 0;
    }

    std::string getTypeName(size_t type) const
    {
      return getCategoryName(getCategoryFrom(type)) + "->" +  getCategoryName(getCategoryTo(type));
    }
  };


/**
 * @brief Distinguishes all types of substitutions.
 *
 * This register has all n * (n-1) substitution type, where n is the size of the alphabet, mapped as:
 * - 0 not a substitution
 * - x in [1, n(n-1)] a substitution
 */
  class ComprehensiveSubstitutionRegister :
    public CategorySubstitutionRegister
  {
  public:
    ComprehensiveSubstitutionRegister(const SubstitutionModel* model, bool within = false) :
      CategorySubstitutionRegister(model, within)
    {
      std::map<int, int> categories;
      for (int i = 0; i < static_cast<int>(model->getAlphabet()->getSize()); ++i)
      {
        categories[i] = i;
      }
      setAlphabetCategories<int>(categories);
    }

    ComprehensiveSubstitutionRegister* clone() const { return new ComprehensiveSubstitutionRegister(*this); }
  };

/**
 * @brief Distinguishes AT<->GC from GC<->AT.
 *
 * This register has two substitution types, mapped as:
 * - 0 not a substitution
 * - 1 a AT->GC substitution
 * - 2 a GC->AT substitution
 */
  class GCSubstitutionRegister :
    public CategorySubstitutionRegister
  {
  public:
    GCSubstitutionRegister(const NucleotideSubstitutionModel* model, bool within = false) :
      CategorySubstitutionRegister(model, within)
    {
      std::map<int, int> categories;
      categories[0] = 1;
      categories[1] = 2;
      categories[2] = 2;
      categories[3] = 1;
      setAlphabetCategories<int>(categories);
    }

    GCSubstitutionRegister* clone() const { return new GCSubstitutionRegister(*this); }
  };


/**
 * @brief Distinguishes AT->GC vs GC->AT inside synonymous
 * substitutions on third codon position.
 *
 * This register has two substitution types, mapped as:
 * - 0 not a counted substitution
 * - 1 a AT->GC synonymous substitution
 * - 2 a GC->AT synonymous substitution
 *
 * Multiple substitutions are forbidden.
 *
 */

  class GCSynonymousSubstitutionRegister :
    public CategorySubstitutionRegister
  {
  private:
    const GeneticCode* code_;

  public:
    GCSynonymousSubstitutionRegister(const CodonSubstitutionModel* model, bool within = false) :
      CategorySubstitutionRegister(model, within),
      code_(model->getGeneticCode())
    {
      const CodonAlphabet* pCA = dynamic_cast<const CodonAlphabet*>(code_->getSourceAlphabet());

      std::map<int, int> categories;
      for (int i = 0; i < static_cast<int>(pCA->getSize()); i++)
      {
        int n = pCA->getThirdPosition(i);
        switch (n)
        {
        case 0:
        case 3:
          categories[i] = 1;
          break;
        case 1:
        case 2:
          categories[i] = 2;
          break;
        }
      }
      setAlphabetCategories<int>(categories);
    }

    GCSynonymousSubstitutionRegister(const GCSynonymousSubstitutionRegister& reg) :
      CategorySubstitutionRegister(reg),
      code_(reg.code_)
    {}

    GCSynonymousSubstitutionRegister& operator=(const GCSynonymousSubstitutionRegister& reg)
    {
      CategorySubstitutionRegister::operator=(reg);
      code_ = reg.code_;
      return *this;
    }

    GCSynonymousSubstitutionRegister* clone() const { return new GCSynonymousSubstitutionRegister(*this); }

  public:
    size_t getNumberOfSubstitutionTypes() const { return 2; }

    size_t getType(size_t fromState, size_t toState) const
    {
      int x = model_->getAlphabetStateAsInt(fromState);
      int y = model_->getAlphabetStateAsInt(toState);
      const CodonAlphabet* pCA = dynamic_cast<const CodonAlphabet*>(code_->getSourceAlphabet());
      if (code_->isStop(x) || code_->isStop( y) || !code_->areSynonymous(x, y))
        return 0;

      // only substitutions between 3rd positions

      if ((pCA->getFirstPosition(x) != pCA->getFirstPosition(y)) ||
          (pCA->getSecondPosition(x) != pCA->getSecondPosition(y)))
        return 0;

      size_t fromCat = categories_[fromState];
      size_t toCat   = categories_[toState];

      if (fromCat > 0 && toCat > 0)
        return index_[fromCat - 1][toCat - 1];
      else
        return 0;
    }

    std::string getTypeName (size_t type) const
    {
      if (type == 0)
      {
        return "no AT<->GC substitution or non-synonymous substitution";
      }
      else if (type == 1)
      {
        return "AT->GC synonymous";
      }
      else if (type == 2)
      {
        return "GC->AT synonymous";
      }
      else
      {
        throw Exception("GCSynonymousSubstitutionRegister::getTypeName. Bad substitution type.");
      }
    }
  };


} // end of namespace bpp.

#endif // _CATEGORY_SUBSTITUTIONREGISTER_H_
