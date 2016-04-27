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

#include "../Model/SubstitutionModel.h"
#include "../Model/Nucleotide/NucleotideSubstitutionModel.h"
#include "../Model/Protein/ProteinSubstitutionModel.h"
#include "../Model/Codon/CodonSubstitutionModel.h"

// From bpp-core:
#include <Bpp/Clonable.h>
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Text/StringTokenizer.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/NucleicAlphabet.h>
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/GeneticCode/GeneticCode.h>

// From the STL:
#include <vector>
#include <string>
#include <algorithm>

namespace bpp
{
/**
 * @brief The SubstitutionRegister interface.
 *
 * Substitution registers are simple classes that define categories of substitutions, and assign an index to them.
 * Substitution registers are defined according to a given substitution model.
 *
 * @author Julien Dutheil
 */
  class SubstitutionRegister :
    public virtual Clonable
  {
  public:
    SubstitutionRegister() {}
    virtual ~SubstitutionRegister() {}

    virtual SubstitutionRegister* clone() const = 0;

  public:
    /**
     * @return The alphabet associated to this instance.
     */
    virtual const Alphabet* getAlphabet() const = 0;

    /**
     * @return The substitution model associated to this instance.
     */
    virtual const SubstitutionModel* getSubstitutionModel() const = 0;

    /**
     * @return The number of substitution types supported by this class.
     */
    virtual size_t getNumberOfSubstitutionTypes() const = 0;

    /**
     * @brief Get the substitution type far a given pair of model states.
     *
     * @param fromState Initial state (should be a state supported by the specified alphabet).
     * @param toState   Final state (should be a state supported by the specified alphabet).
     * @return The index of the corresponding substitution type, ranging from 0 to 'getNumberOfSubstitutionTypes' + 1,
     * as non-substitution (that is when fromState == toState) will always return 0.
     */
    virtual size_t getType(size_t fromState, size_t toState) const = 0;

    /**
     * @brief Get the name of a given substitution type.
     *
     * This method is only used for user-friendlyness purposes, not computational goal.
     * I can therefore be left unimplemented in some cases.
     *
     * @param type Index of the substitution (should be an size_t contained in the register).
     * @return A string describing the substitution type.
     */
    virtual std::string getTypeName(size_t type) const = 0;
  };

  class AbstractSubstitutionRegister :
    public virtual SubstitutionRegister
  {
  protected:
    const SubstitutionModel* model_;

  public:
    AbstractSubstitutionRegister(const SubstitutionModel* model) :
      model_(model)
    {}

    AbstractSubstitutionRegister(const AbstractSubstitutionRegister& asr) :
      model_(asr.model_)
    {}

    AbstractSubstitutionRegister& operator=(const AbstractSubstitutionRegister& asr)
    {
      model_ = asr.model_;
      return *this;
    }

    virtual ~AbstractSubstitutionRegister() {}

  public:
    const SubstitutionModel* getSubstitutionModel() const { return model_; }

    const Alphabet* getAlphabet() const { return model_->getAlphabet(); }
  };

/**
 * @brief Count all substitutions.
 *
 * This register has only 1 substitution type, mapped as:
 * - 0 not a substitution
 * - 1 a substitution
 */
  class TotalSubstitutionRegister :
    public AbstractSubstitutionRegister
  {
  public:
    TotalSubstitutionRegister(const SubstitutionModel* model) :
      AbstractSubstitutionRegister(model)
    {}

    TotalSubstitutionRegister* clone() const { return new TotalSubstitutionRegister(*this); }

  public:
    size_t getNumberOfSubstitutionTypes() const { return 1; }

    size_t getType(size_t fromState, size_t toState) const
    {
      return fromState == toState ? 0 : 1;
    }

    std::string getTypeName(size_t type) const
    {
      if (type == 0)
      {
        return "no substitution";
      }
      else if (type == 1)
      {
        return "substitution";
      }
      else
      {
        throw Exception("TotalSubstitutionRegister::getTypeName. Bad substitution type.");
      }
    }
  };

/**
 * @brief Completion of a given substitution register to consider
 * all substitutions. The new substitutions are considered in an
 * additional type.
 */
  class CompleteSubstitutionRegister :
    public AbstractSubstitutionRegister
  {
  private:
    const SubstitutionRegister* preg_;

    bool isRegComplete_;

  public:
    CompleteSubstitutionRegister(const SubstitutionRegister& reg) :
      AbstractSubstitutionRegister(reg.getSubstitutionModel()),
      preg_(reg.clone()),
      isRegComplete_(true)
    {
      size_t size = reg.getAlphabet()->getSize();
      for (size_t i = 0; i < size; i++)
      {
        for (size_t j = 0; j < size; j++)
        {
          if ((i != j) && reg.getType(i, j) == 0)
          {
            isRegComplete_ = false;
            return;
          }
        }
      }
    }

    CompleteSubstitutionRegister* clone() const { return new CompleteSubstitutionRegister(*this); }

    CompleteSubstitutionRegister(const CompleteSubstitutionRegister& csr) :
      AbstractSubstitutionRegister(csr),
      preg_(csr.preg_->clone()),
      isRegComplete_(csr.isRegComplete_)
    {}

    CompleteSubstitutionRegister& operator=(const CompleteSubstitutionRegister& csr)
    {
      AbstractSubstitutionRegister::operator=(csr);
      preg_ = csr.preg_->clone();
      isRegComplete_ = csr.isRegComplete_;
      return *this;
    }

    ~CompleteSubstitutionRegister()
    {
      if (preg_)
        delete preg_;
      preg_ = 0;
    }

  public:
    size_t getNumberOfSubstitutionTypes() const
    {
      return preg_->getNumberOfSubstitutionTypes() + (isRegComplete_ ? 0 : 1);
    }

    size_t getType(size_t fromState, size_t toState) const
    {
      if (fromState==toState)
        return 0;
      
      size_t t = preg_->getType(fromState, toState);
      if (t == 0)
        return getNumberOfSubstitutionTypes();
      else
        return t;
    }

    std::string getTypeName(size_t type) const
    {
      try
      {
        return preg_->getTypeName(type);
      }
      catch (Exception& e)
      {
        if (type == getNumberOfSubstitutionTypes())
          return "Completion substitution";
        else
          throw Exception("CompleteSubstitutionRegister::getTypeName. Bad substitution type.");
      }
    }
  };


  /**
   * @brief Sets a Register based on a vector of Registers. The
   * categories are intersection of categories of those Registers.
   *
   * @author Laurent Guéguen
   *
   */

  class VectorOfSubstitionRegisters :
    public AbstractSubstitutionRegister
  {
  private:
    /**
     * @brief the vector of pointers to SubstitutionRegisters.
     *
     *  These SubstitutionRegisters belong to the object.
     */

    std::vector<SubstitutionRegister*> vSubReg_;

  public:

    VectorOfSubstitionRegisters(const SubstitutionModel* model) :
      AbstractSubstitutionRegister(model),
      vSubReg_()
    {}

    VectorOfSubstitionRegisters(const VectorOfSubstitionRegisters& vosr) :
      AbstractSubstitutionRegister(vosr),
      vSubReg_()
    {
      for (size_t i = 0; i < vosr.vSubReg_.size(); i++)
      {
        vSubReg_.push_back(vosr.vSubReg_[i]->clone());
      }
    }
    
    VectorOfSubstitionRegisters& operator=(const VectorOfSubstitionRegisters& vosr)
    {
      AbstractSubstitutionRegister::operator=(vosr);

      for (size_t i = 0; i < vSubReg_.size(); i++)
        delete vSubReg_[i];

      vSubReg_.clear();
      
      for (size_t i = 0; i < vosr.vSubReg_.size(); i++)
        vSubReg_.push_back(vosr.vSubReg_[i]->clone());

      return *this;
    }

    VectorOfSubstitionRegisters* clone() const
    {
      return new VectorOfSubstitionRegisters(*this);
    }
    
    ~VectorOfSubstitionRegisters()
    {
      for (size_t i = 0; i < vSubReg_.size(); i++)
        delete vSubReg_[i];

      vSubReg_.clear();
    }

    void addRegister(SubstitutionRegister* reg)
    {
      if (reg){
        if (reg->getSubstitutionModel()!=getSubstitutionModel())
          throw Exception("VectorOfSubstitionRegisters::addRegister : mismatch between models");
      
        vSubReg_.push_back(reg);
      }
    }
    
    size_t getType(size_t i, size_t j) const
    {
      if (i==j)
        return 0;

      size_t x=0;

      for (size_t p=0;p<vSubReg_.size();p++)
      {
        x*=vSubReg_[p]->getNumberOfSubstitutionTypes();
        size_t z=vSubReg_[p]->getType(i,j);
        if (z==0)
          return 0;
        x+=z-1;
      }

      return x+1;
    }

    size_t getNumberOfSubstitutionTypes() const
    {
      if (vSubReg_.size()==0)
        return 0;
      
      size_t n=1;
      
      for (size_t i=0;i<vSubReg_.size();i++)
        n*=vSubReg_[i]->getNumberOfSubstitutionTypes();

      return n;
    }

    /**
     * @brief names of the types are the list of their types in the
     * registers (separated with _).
     */
    
    std::string getTypeName(size_t type) const
    {
      if (type == 0)
        return "no substitution";

      std::string res="";
      size_t ty=type-1;
      
      for (size_t p=vSubReg_.size(); p>0; p--)
      {
        res+=vSubReg_[p-1]->getTypeName((ty%vSubReg_[p-1]->getNumberOfSubstitutionTypes())+1);

        ty/=vSubReg_[p-1]->getNumberOfSubstitutionTypes();
        
        if (p!=1)
          res+="_";
      }

      return res;
    }

  };
  
    
/**
 * @brief Sets a Register based on a matrix of integers. If M is the
 *  matrix, M[i,j] is the number of the substitution type from i to j,
 *  or 0 if there is no substitution type from i to j.
 *
 * @author Laurent Guéguen
 */
  class GeneralSubstitutionRegister :
    public AbstractSubstitutionRegister
  {
  protected:
    /**
     * @brief The size of the matrix, i.e. the number of states.
     */
    size_t size_;

    /**
     * @brief The matrix of the substitution register.
     */

    RowMatrix<size_t> matrix_;

    /**
     * @brief The map from substitution types to the map of from states
     * to the vector of target states.
     *
     * This is the reverse information of matrix_
     *
     */

    std::map<size_t, std::map<size_t, std::vector<size_t> > > types_;

  public:
    GeneralSubstitutionRegister(const SubstitutionModel* model) :
      AbstractSubstitutionRegister(model),
      size_(model->getNumberOfStates()),
      matrix_(size_, size_),
      types_()
    {}

    GeneralSubstitutionRegister(const SubstitutionModel* model, const RowMatrix<size_t>& matrix) :
      AbstractSubstitutionRegister(model),
      size_(model->getNumberOfStates()),
      matrix_(matrix),
      types_()
    {
      if (matrix_.getNumberOfRows() != size_)
        throw DimensionException("GeneralSubstitutionRegister", size_, matrix_.getNumberOfRows());
      if (matrix_.getNumberOfColumns() != size_)
        throw DimensionException("GeneralSubstitutionRegister", size_, matrix_.getNumberOfColumns());
      updateTypes_();
    }

    GeneralSubstitutionRegister(const GeneralSubstitutionRegister& gsr) :
      AbstractSubstitutionRegister(gsr),
      size_(gsr.size_),
      matrix_(gsr.matrix_),
      types_(gsr.types_)
    {}

    GeneralSubstitutionRegister& operator=(const GeneralSubstitutionRegister& gsr)
    {
      AbstractSubstitutionRegister::operator=(gsr);
      size_   = gsr.size_;
      matrix_ = gsr.matrix_;
      types_  = gsr.types_;
      return *this;
    }

    GeneralSubstitutionRegister* clone() const { return new GeneralSubstitutionRegister(*this); }

    virtual ~GeneralSubstitutionRegister() {}

    size_t getType(size_t i, size_t j) const
    {
      return matrix_(i, j);
    }

    size_t getNumberOfSubstitutionTypes() const
    {
      return types_.find(0) == types_.end() ? types_.size() : types_.size() - 1;
    }

    /**
     * @brief names of the types are their number.
     */
    std::string getTypeName(size_t type) const
    {
      if (types_.find(type) != types_.end())
        return TextTools::toString(type);

      throw Exception("Bad type number " + TextTools::toString(type) + " in GeneralSubstitutionRegister::getTypeName.");
    }

  protected:
    void updateTypes_();
  };


/**
 *  @brief Class inheriting from GeneralSubstitutionRegister, this one
 *  uses a special constructor which allows it to build a substitution
 *  matrix from string input specifying a desired substitutions.
 *
 *  @author Juraj Michalik
 */

  class SelectedSubstitutionRegister :
    public GeneralSubstitutionRegister
  {
    std::map<size_t, std::string> categoryNames_;

  public:
    SelectedSubstitutionRegister (const SubstitutionModel* model, std::string selectedSubstitutions) :
      GeneralSubstitutionRegister(model),
      categoryNames_()
    {
      selectedSubstitutions.erase(std::remove(selectedSubstitutions.begin(), selectedSubstitutions.end(), ' '), selectedSubstitutions.end());
      /**
       *  @brief This constructor creates an empty square matrix (nrow = ncol
       *  = length of alphabet) and takes a string with specific
       *  syntax to mark a substitutions with a certain index
       *  depending on the string entered.
       *
       *  The same group of substitution is delimited by parentheses.
       *  The name, if entered, is entered at the start of a string
       *  and followed by ";". Substitutions are delimited by ",", and
       *  each substitution is defined with a "->" symbol.
       */
      size_t typeSubs = 0;
      size_t coord1 = 0;
      size_t coord2 = 0;
      std::string codon1 = "";
      std::string codon2 = "";
      StringTokenizer subsGroup(selectedSubstitutions, "()");
      while (subsGroup.hasMoreToken())
      {
        typeSubs++;
        StringTokenizer namesSubs(subsGroup.nextToken(), ":");
        if (namesSubs.numberOfRemainingTokens() == 2)
        {
          categoryNames_[typeSubs] = namesSubs.nextToken();
        }
        else if (namesSubs.numberOfRemainingTokens() == 1)
        {
          categoryNames_[typeSubs] = TextTools::toString(typeSubs);
        }
        StringTokenizer substitutions(namesSubs.nextToken(), ",");
        while (substitutions.hasMoreToken())
        {
          StringTokenizer coordinates(substitutions.nextToken(), "->");
          codon1 = coordinates.nextToken();
          codon2 = coordinates.nextToken();
          coord1 = model_->getAlphabet()->getStateIndex(codon1);
          coord2 = model_->getAlphabet()->getStateIndex(codon2);
          this->matrix_(coord1, coord2) = typeSubs;
        }
      }
      updateTypes_();
    }

    SelectedSubstitutionRegister* clone() const { return new SelectedSubstitutionRegister(*this); }

    ~SelectedSubstitutionRegister() {}

    std::string getTypeName(size_t type) const
    {
      if (types_.find(type) != types_.end())
        return TextTools::toString(categoryNames_.find(type)->second);

      throw Exception("Bad type number " + TextTools::toString(type) + " in GeneralSubstitutionRegister::getTypeName.");
    }
  };

/**
 * @brief Indexes only intra amino-acid substitutions. Every group
 * represents a substitutions for the same amino acid. Met and Trp are
 * not taken into account due their non-degenerescence.
 *
 * @author Juraj Michalik
 */

  class AAInteriorSubstitutionRegister :
    public GeneralSubstitutionRegister
  {
    std::map<std::string, size_t> categoryCorrespondance_;

  public:
    AAInteriorSubstitutionRegister(const CodonSubstitutionModel* model) :
      GeneralSubstitutionRegister(model),
      categoryCorrespondance_()
    {
      size_t categoryIndex = 1;
      for (size_t i = 1; i <= model->getAlphabet()->getSize(); ++i)
      {
        int state1 = model->getAlphabet()->getStateAt(i).getNum();
        for (size_t j = i + 1; j <= model->getAlphabet()->getSize(); ++j)
        {
          int state2 = model->getAlphabet()->getStateAt(j).getNum();
          if (!(model->getGeneticCode()->isStop(state1)) && !(model->getGeneticCode()->isStop(state2)))
          {
            if (model->getGeneticCode()->translate(state1) == model->getGeneticCode()->translate(state2))
            {
              std::string aminoAcid = model->getGeneticCode()->getTargetAlphabet()->intToChar(model->getGeneticCode()->translate(state1));
              if (categoryCorrespondance_.find(aminoAcid) == categoryCorrespondance_.end())
              {
                categoryCorrespondance_[aminoAcid] = categoryIndex;
                categoryIndex++;
              }
              matrix_(i, j) = categoryCorrespondance_[aminoAcid];
              matrix_(j, i) = categoryCorrespondance_[aminoAcid];
            }
          }
        }
      }
      updateTypes_();
    }

    AAInteriorSubstitutionRegister* clone() const { return new AAInteriorSubstitutionRegister(*this); }

    ~AAInteriorSubstitutionRegister() {}


    std::string getTypeName(size_t type) const
    {
      if (types_.find(type) != types_.end())
      {
        for (std::map<std::string, size_t>::const_iterator it = categoryCorrespondance_.begin(); it != categoryCorrespondance_.end(); it++)
        {
          if (it->second == type)
            return TextTools::toString(it->first);
        }
      }
      throw Exception("Bad type number " + TextTools::toString(type) + " in GeneralSubstitutionRegister::getTypeName.");
    }
  };

/**
 * @brief Indexes only substitutions between amino-acids. Groups are
 * distinguished by their direction.
 *
 * @author Juraj Michalik
 */

  class AAExteriorSubstitutionRegister :
    public GeneralSubstitutionRegister
  {
    std::map<std::string, size_t> categoryCorrespondance_;

  public:
    AAExteriorSubstitutionRegister (const CodonSubstitutionModel* model) :
      GeneralSubstitutionRegister(model),
      categoryCorrespondance_()
    {
      size_t categoryIndex = 1;
      for (size_t i = 1; i <= model->getAlphabet()->getSize(); ++i)
      {
        int state1 = model->getAlphabet()->getStateAt(i).getNum();
        for (size_t j = i + 1; j <= model->getAlphabet()->getSize(); ++j)
        {
          int state2 = model->getAlphabet()->getStateAt(j).getNum();
          if (!(model->getGeneticCode()->isStop(state1)) && !(model->getGeneticCode()->isStop(state2)))
          {
            if (model->getGeneticCode()->translate(state1) != model->getGeneticCode()->translate(state2))
            {
              std::string aminoAcid1 = model->getGeneticCode()->getTargetAlphabet()->intToChar(model->getGeneticCode()->translate(state1));
              std::string aminoAcid2 = model->getGeneticCode()->getTargetAlphabet()->intToChar(model->getGeneticCode()->translate(state2));
              bool AA1IsNotInGroup = ((categoryCorrespondance_.find(aminoAcid1 + "->" + aminoAcid2) == categoryCorrespondance_.end()));
              bool AA2IsNotInGroup = ((categoryCorrespondance_.find(aminoAcid2 + "->" + aminoAcid1) == categoryCorrespondance_.end()));
              if (AA1IsNotInGroup)
              {
                categoryCorrespondance_[aminoAcid1 + "->" + aminoAcid2] = categoryIndex;
                categoryIndex++;
              }
              this->matrix_(i, j) = categoryCorrespondance_[aminoAcid1 + "->" + aminoAcid2];
              if (AA2IsNotInGroup)
              {
                categoryCorrespondance_[aminoAcid2 + "->" + aminoAcid1] = categoryIndex;
                categoryIndex++;
              }
              this->matrix_(j, i) = categoryCorrespondance_[aminoAcid2 + "->" + aminoAcid1];
            }
          }
        }
      }
      updateTypes_();
    }

    AAExteriorSubstitutionRegister* clone() const { return new AAExteriorSubstitutionRegister(*this); }

    ~AAExteriorSubstitutionRegister() {}

    std::string getTypeName(size_t type) const
    {
      if (types_.find(type) != types_.end())
      {
        for (std::map<std::string, size_t>::const_iterator it = categoryCorrespondance_.begin(); it != categoryCorrespondance_.end(); it++)
        {
          if (it->second == type)
            return TextTools::toString(it->first);
        }
      }
      throw Exception("Bad type number " + TextTools::toString(type) + " in GeneralSubstitutionRegister::getTypeName.");
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

  class TsTvSubstitutionRegister :
    public AbstractSubstitutionRegister
  {
  private:
    /**
     *  @brief useful for codon alphabet
     *
     */
    
    const GeneticCode* code_;

  public:
    TsTvSubstitutionRegister(const NucleotideSubstitutionModel* model) :
      AbstractSubstitutionRegister(model),
      code_()
    {}

    TsTvSubstitutionRegister(const CodonSubstitutionModel* model) :
      AbstractSubstitutionRegister(model),
      code_(model->getGeneticCode())
    {}

    TsTvSubstitutionRegister(const TsTvSubstitutionRegister& reg) :
      AbstractSubstitutionRegister(reg),
      code_(reg.code_)
    {}

    TsTvSubstitutionRegister& operator=(const TsTvSubstitutionRegister& reg)
    {
      AbstractSubstitutionRegister::operator=(reg);
      code_ = reg.code_;
      return *this;
    }

    TsTvSubstitutionRegister* clone() const { return new TsTvSubstitutionRegister(*this); }

  public:
    size_t getNumberOfSubstitutionTypes() const { return 2; }

    size_t getType(size_t fromState, size_t toState) const
    {
      int x = model_->getAlphabetStateAsInt(fromState);
      int y = model_->getAlphabetStateAsInt(toState);
      if (x == y)
        return 0;                     // nothing happens

      const CodonAlphabet* cAlpha = dynamic_cast<const CodonAlphabet*>(model_->getAlphabet());
      if (cAlpha)
      {
        if (code_->getSourceAlphabet()->isGap(x)
            || code_->getSourceAlphabet()->isGap(y)
            || code_->isStop(x)
            || code_->isStop(y))
          return 0;
        
        size_t countPos = 0;
        if (cAlpha->getFirstPosition(x) != cAlpha->getFirstPosition(y))
          countPos++;
        if (cAlpha->getSecondPosition(x) != cAlpha->getSecondPosition(y))
          countPos++;
        if (cAlpha->getThirdPosition(x) != cAlpha->getThirdPosition(y))
          countPos++;
        if (countPos > 1)
          return 0;
      }

      int d=abs(x-y);

      while (d!=0)
      {
        if (d==2)
          return 1;                     // This is a transition
        else
          d/=4;
      }
      
      return 2; // This is a transversion
    }

    std::string getTypeName(size_t type) const
    {
      if (type == 0)
      {
        return "no substitution";
      }
      else if (type == 1)
      {
        return "transition";
      }
      else if (type == 2)
      {
        return "transversion";
      }
      else
      {
        throw Exception("TsTvSubstitutionRegister::getTypeName. Bad substitution type.");
      }
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
  class DnDsSubstitutionRegister :
    public AbstractSubstitutionRegister
  {
  private:
    const GeneticCode* code_;
    bool countMultiple_;

  public:
    DnDsSubstitutionRegister(const CodonSubstitutionModel* model, bool countMultiple = false) :
      AbstractSubstitutionRegister(model),
      code_(model->getGeneticCode()),
      countMultiple_(countMultiple)
    {}

    DnDsSubstitutionRegister(const DnDsSubstitutionRegister& reg) :
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
    size_t getNumberOfSubstitutionTypes() const { return 2; }

    size_t getType(size_t fromState, size_t toState) const
    {
      int x = model_->getAlphabetStateAsInt(fromState);
      int y = model_->getAlphabetStateAsInt(toState);
      const CodonAlphabet* cAlpha = dynamic_cast<const CodonAlphabet*>(model_->getAlphabet());
      if (code_->getSourceAlphabet()->isGap(x)
          || code_->getSourceAlphabet()->isGap(y)
          || code_->isStop(x)
          || code_->isStop(y))
        return 0;
      if (x == y)
        return 0;                     // nothing happens
      if (!countMultiple_)
      {
        size_t countPos = 0;
        if (cAlpha->getFirstPosition(x) != cAlpha->getFirstPosition(y))
          countPos++;
        if (cAlpha->getSecondPosition(x) != cAlpha->getSecondPosition(y))
          countPos++;
        if (cAlpha->getThirdPosition(x) != cAlpha->getThirdPosition(y))
          countPos++;
        if (countPos > 1)
          return 0;
      }
      return code_->areSynonymous(x, y) ? 1 : 2;
    }

    std::string getTypeName (size_t type) const
    {
      if (type == 0)
      {
        return "no substitution";
      }
      else if (type == 1)
      {
        return "synonymous";
      }
      else if (type == 2)
      {
        return "non synonymous";
      }
      else
      {
        throw Exception("DnDsSubstitutionRegister::getTypeName. Bad substitution type.");
      }
    }
  };

/**
 * @brief Count conservative and radical amino-acid substitutions.
 *
 * Categories are defined as in Sainudiin et al (2005).
 * - 1 a conservative (C) substitution
 * - 2 a radical (R) substitution
 *
 * @author Jonathan Romiguier, Emeric Figuet, Julien Dutheil
 */
  class KrKcSubstitutionRegister :
    public AbstractSubstitutionRegister
  {
  private:
    std::vector< std::vector<bool> > types_;

  public:
    KrKcSubstitutionRegister(const ProteinSubstitutionModel* model) :
      AbstractSubstitutionRegister(model),
      types_(20)
    {
      for (size_t i = 0; i < 20; ++i) {
        types_[i].resize(20);
        for (size_t j = 0; j < 20; ++j) {
          types_[i][j] = false;
        }
      }

      types_[0][7] = true;
      types_[0][14] = true;
      types_[0][19] = true;
      types_[7][0] = true;
      types_[7][14] = true;
      types_[7][19] = true;
      types_[14][0] = true;
      types_[14][7] = true;
      types_[14][19] = true;
      types_[19][0] = true;
      types_[19][7] = true;
      types_[19][14] = true;
    
      types_[1][5] = true;
      types_[1][6] = true;
      types_[1][8] = true;
      types_[1][11] = true;
      types_[1][17] = true;
      types_[1][18] = true;
      types_[5][6] = true;
      types_[5][8] = true;
      types_[5][11] = true;
      types_[5][1] = true;
      types_[5][17] = true;
      types_[5][18] = true;
      types_[6][8] = true;
      types_[6][11] = true;
      types_[6][5] = true;
      types_[6][1] = true;
      types_[6][17] = true;
      types_[6][18] = true;
      types_[8][6] = true;
      types_[8][11] = true;
      types_[8][5] = true;
      types_[8][1] = true;
      types_[8][17] = true;
      types_[8][18] = true;
      types_[11][1] = true;
      types_[11][5] = true;
      types_[11][6] = true;
      types_[11][8] = true;
      types_[11][17] = true;
      types_[11][18] = true;
      types_[17][1] = true;
      types_[17][5] = true;
      types_[17][6] = true;
      types_[17][8] = true;
      types_[17][11] = true;
      types_[17][18] = true;
      types_[18][1] = true;
      types_[18][6] = true;
      types_[18][8] = true;
      types_[18][11] = true;
      types_[18][5] = true;
      types_[18][17] = true;
    
      types_[2][3] = true;
      types_[2][4] = true;
      types_[2][15] = true;
      types_[2][16] = true;
      types_[3][2] = true;
      types_[3][4] = true;
      types_[3][15] = true;
      types_[3][16] = true;
      types_[4][2] = true;
      types_[4][3] = true;
      types_[4][15] = true;
      types_[4][16] = true;
      types_[15][2] = true;
      types_[15][3] = true;
      types_[15][4] = true;
      types_[15][16] = true;
      types_[16][2] = true;
      types_[16][3] = true;
      types_[16][4] = true;
      types_[16][15] = true;
    
      types_[9][10] = true;
      types_[9][12] = true;
      types_[9][13] = true;
      types_[10][9] = true;
      types_[10][12] = true;
      types_[10][13] = true;
      types_[12][9] = true;
      types_[12][10] = true;
      types_[12][13] = true;
      types_[13][9] = true;
      types_[13][10] = true;
      types_[13][12] = true;
    }

    KrKcSubstitutionRegister* clone() const { return new KrKcSubstitutionRegister(*this); }

  public:
    size_t getNumberOfSubstitutionTypes() const { return 2; }

    size_t getType(size_t fromState, size_t toState) const
    {
      if (fromState == toState)
        return 0;  // nothing happens
      int x = model_->getAlphabetStateAsInt(fromState);
      int y = model_->getAlphabetStateAsInt(toState);
      return types_[static_cast<size_t>(x)][static_cast<size_t>(y)] ? 1 : 2;
    }

    std::string getTypeName(size_t type) const
    {
      if (type == 0)
      {
        return "no substitution";
      }
      else if (type == 1)
      {
        return "conservative";
      }
      else if (type == 2)
      {
        return "radical";
      }
      else
      {
        throw Exception("KrKcSubstitutionRegister::getTypeName. Bad substitution type.");
      }
    }
  };



} // end of namespace bpp.

#endif // _SUBSTITUTIONREGISTER_H_
