// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_NUCLEOTIDE_YPR_H
#define BPP_PHYL_MODEL_NUCLEOTIDE_YPR_H

#include <Bpp/Seq/Alphabet/RNY.h>

#include "../AbstractSubstitutionModel.h"
#include "NucleotideSubstitutionModel.h"

// From Utils:
#include <Bpp/Exceptions.h>

using namespace std;

namespace bpp
{
/**
 * @brief YpR  model.
 * @author Laurent Guéguen
 *
 * Model YpR, on RNY triplets, with independent positions.
 *
 * This model is made from neighbourhood parameters and
 * a nucleotidic probabilistic model such that the R/Y
 * condition is respected : with letters A, T, C, G
 * @f[
 * M=\begin{pmatrix}
 * . & \beta_T & \beta_C & \alpha_G \\
 * \beta_A & . & \alpha_C & \beta_G \\
 * \beta_A & \alpha_T & . & \beta_G \\
 * \alpha_A & \beta_T & \beta_C & .
 * \end{pmatrix}
 * @f]
 *
 * From this model, the rates are a multiplication of the rates for the
 * letters of the triplet.
 * For the first letter, on alphabet R, C, T:
 * @f[
 * M_1=\begin{pmatrix}
 * . & \beta_C & \beta_T \\
 * \beta_A + \beta_G & . & \alpha_T \\
 * \beta_A + \beta_G & \alpha_C & .
 * \end{pmatrix}
 * @f]
 * For the second letter, on alphabet A, G, C, T:
 * @f[
 * M_2=\begin{pmatrix}
 * . & \alpha_G & \beta_C & \beta_T \\
 * \alpha_A & . & \beta_C & \beta_T \\
 * \beta_A & \beta_G & . & \alpha_T \\
 * \beta_A & \beta_G & \alpha_C & .
 * \end{pmatrix}
 * @f]
 * For the third letter, on alphabet A, G, Y:
 * @f[
 * M_3 = \begin{pmatrix}
 * . & \alpha_G & \beta_C + \beta_T\\
 * \alpha_A & . & \beta_C + \beta_T\\
 * \beta_A & \beta_G & .
 * \end{pmatrix}.
 * @f]
 * And the model is the "union" of the 3 matrices, since in the
 * model  only the mutations concerning only one position  are
 * not null.
 * In addition to this, neighbour dependency parameters
 * (in inherited models) give the extra mutation rates for transitions
 * inside YpR dinucleotides.
 * For example from CpG to CpA and TpG, in relative proportion
 * to the single transition rate:
 *  CGx -> CAx   += rcGA*G->A
 *  CGx -> TGx   += rCgT*C->T
 *  xCG -> xCA   += rcGA*G->A
 *  xCG -> xTG   += rCgT*C->T
 *
 * Other YpR neighbour dependency rates are:
 *   TG -> CG, TG -> TA, CA -> TA, CA -> CG, TA -> CA, TA -> TC
 *
 * The generator in normalized such that there is on average
 *  a substitution per site per unit of time ON THE CENTRAL POSITION
 *  of the triplet.
 * @see AbstractSubstitutionModel
 */
class YpR :
  public AbstractSubstitutionModel
{
protected:
  std::unique_ptr<NucleotideSubstitutionModelInterface> pmodel_;

  std::string nestedPrefix_;

protected:
  /**
   * @brief Build a new YpR substitution model, with no dependency
   *   parameters
   */
  YpR(
      std::shared_ptr<const RNY>,
      std::unique_ptr<NucleotideSubstitutionModelInterface> const,
      const std::string& prefix);

  YpR(const YpR&, const std::string& prefix);

  YpR(const YpR& ypr);

  YpR& operator=(const YpR& ypr)
  {
    AbstractParameterAliasable::operator=(ypr);
    AbstractSubstitutionModel::operator=(ypr);
    nestedPrefix_ = ypr.nestedPrefix_;
    pmodel_.reset(ypr.pmodel_->clone());
    return *this;
  }

public:
  virtual ~YpR() {}

protected:
  void updateMatrices_(double, double, double, double,
      double, double, double, double);

  virtual void updateMatrices_() override;

  string getNestedPrefix() const
  {
    return nestedPrefix_;
  }

public:
  //  virtual std::string getName() const;

  const NucleotideSubstitutionModelInterface& nestedModel() const { return *pmodel_; }

  size_t getNumberOfStates() const override { return 36; }

  virtual void setNamespace(const std::string&) override;

  void fireParameterChanged(const ParameterList& parameters) override
  {
    AbstractSubstitutionModel::fireParameterChanged(parameters);
    pmodel_->matchParametersValues(parameters);
    updateMatrices_();
  }

  // Check that the model is good for YpR
  void checkModel(const SubstitutionModelInterface& model) const;
};
}


// //////////////////////////////////////
// //////// YpR_symetrical

namespace bpp
{
/**
 * @brief symetrical YpR  model.
 *
 * Model YpR, on RNY triplets with symetrical dependency parameters
 *
 * The neighbour dependency parameters are noted as: XyZ for
 *  XpY -> ZpY substitution. Since the process are symetrical,
 *  only these are necessary.
 * @see YpR
 *
 */

class YpR_Sym :
  public YpR
{
public:
  /**
   * @brief Build a new YpR_Sym substitution model.
   * @param CgT, TgC, CaT, TaC neighbour dependency parameters
   * @param alph RNY alphabet
   * @param pm Substitution model.
   */
  YpR_Sym(
      std::shared_ptr<const RNY> alph,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pm,
      double CgT = 0., double TgC = 0.,
      double CaT = 0., double TaC = 0.);

  YpR_Sym(const YpR_Sym&);

  virtual ~YpR_Sym() {}

  YpR_Sym* clone() const override { return new YpR_Sym(*this); }

  std::string getName() const override;

protected:
  void updateMatrices_() override;
};

// //////////////////////////////////////
// //////// YpR_general

/**
 * @brief General YpR  model.
 *
 * Model YpR, on RNY triplets with general dependency parameters
 *
 * The neighbour dependency parameters are noted as: XyZ for
 *  XpY -> ZpY substitution and xYZ for XpY -> XpZ substitution.
 *
 * @see YpR
 */
class YpR_Gen :
  public YpR
{
public:
  /**
   * @brief Build a new YpR_Gen substitution model.
   * @param CgT, cGA, TgC, tGA, CaT, cAG, TaC, tAG neighbour
   * dependency parameters
   * @param alph RNY alphabet
   * @param pm Substitution model.
   */
  YpR_Gen(
      std::shared_ptr<const RNY> alph,
      std::unique_ptr<NucleotideSubstitutionModelInterface> pm,
      double CgT = 0., double cGA = 0.,
      double TgC = 0., double tGA = 0.,
      double CaT = 0., double cAG = 0.,
      double TaC = 0., double tAG = 0.);

  YpR_Gen(const YpR_Gen&);

  virtual ~YpR_Gen() {}

  YpR_Gen* clone() const override { return new YpR_Gen(*this); }

  std::string getName() const override;

protected:
  void updateMatrices_() override;
};
}
#endif // BPP_PHYL_MODEL_NUCLEOTIDE_YPR_H
