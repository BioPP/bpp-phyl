//
// File: AbstractSubstitutionModel.cpp
// Created by: Julien Dutheil
// Created on: Tue May 27 10:31:49 2003
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

#include "AbstractSubstitutionModel.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Matrix/EigenValue.h>
#include <Bpp/Numeric/NumConstants.h>

// From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

AbstractSubstitutionModel::AbstractSubstitutionModel(const Alphabet* alpha, StateMap* stateMap, const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  alphabet_(alpha),
  stateMap_(stateMap),
  size_(alpha->getSize()),
  rate_(1),
  generator_(size_, size_),
  freq_(size_),
  exchangeability_(size_, size_),
  pijt_(size_, size_),
  dpijt_(size_, size_),
  d2pijt_(size_, size_),
  eigenDecompose_(true),
  eigenValues_(size_),
  iEigenValues_(size_),
  isDiagonalizable_(false),
  rightEigenVectors_(size_, size_),
  isNonSingular_(false),
  leftEigenVectors_(size_, size_),
  vPowGen_(),
  tmpMat_(size_, size_)
{
  for (size_t i = 0; i < size_; i++)
  {
    freq_[i] = 1.0 / static_cast<double>(size_);
  }
}

/******************************************************************************/

void AbstractSubstitutionModel::updateMatrices()
{
  // if the object is not an AbstractReversibleSubstitutionModel,
  // computes the exchangeability_ Matrix (otherwise the generator_
  // has been computed from the exchangeability_)

  if (!dynamic_cast<AbstractReversibleSubstitutionModel*>(this)) {
    for (size_t i = 0; i < size_; i++)
    {
      for (size_t j = 0; j < size_; j++)
      {
        exchangeability_(i, j) = generator_(i, j) / freq_[j];
      }
    }
  }

  // Compute eigen values and vectors:
  if (enableEigenDecomposition())
  {
    EigenValue<double> ev(generator_);
    rightEigenVectors_ = ev.getV();
    eigenValues_ = ev.getRealEigenValues();
    iEigenValues_ = ev.getImagEigenValues();
    try
    {
      MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
      isNonSingular_ = true;
      isDiagonalizable_ = true;
      for (size_t i = 0; i < size_ && isDiagonalizable_; i++)
      {
        if (abs(iEigenValues_[i]) > NumConstants::TINY())
          isDiagonalizable_ = false;
      }
    }
    catch (ZeroDivisionException& e)
    {
      ApplicationTools::displayMessage("Singularity during diagonalization. Taylor series used instead.");

      isNonSingular_ = false;
      isDiagonalizable_ = false;
      MatrixTools::Taylor(generator_, 30, vPowGen_);
    }
  }
}


/******************************************************************************/

const Matrix<double>& AbstractSubstitutionModel::getPij_t(double t) const
{
  if (t ==0)
  {
    MatrixTools::getId(size_, pijt_);
  }
  else if (isNonSingular_)
  {
    if (isDiagonalizable_)
    {
      MatrixTools::mult<double>(rightEigenVectors_, VectorTools::exp(eigenValues_ * (rate_ * t)), leftEigenVectors_, pijt_);
    }
    else
    {
      std::vector<double> vdia(size_);
      std::vector<double> vup(size_ - 1);
      std::vector<double> vlo(size_ - 1);
      double c = 0, s = 0;
      double l = rate_ * t;
      for (size_t i = 0; i < size_; i++)
      {
        vdia[i] = std::exp(eigenValues_[i] * l);
        if (iEigenValues_[i] != 0)
        {
          s = std::sin(iEigenValues_[i] * l);
          c = std::cos(iEigenValues_[i] * l);
          vup[i] = vdia[i] * s;
          vlo[i] = -vup[i];
          vdia[i] *= c;
          vdia[i + 1] = vdia[i]; // trick to avoid computation
          i++;
        }
        else
        {
          if (i < size_ - 1)
          {
            vup[i] = 0;
            vlo[i] = 0;
          }
        }
      }
      MatrixTools::mult<double>(rightEigenVectors_, vdia, vup, vlo, leftEigenVectors_, pijt_);
    }
  }
  else
  {
    MatrixTools::getId(size_, pijt_);
    double s = 1.0;
    double v = rate_ * t;
    size_t m = 0;
    while (v > 0.5)    // exp(r*t*A)=(exp(r*t/(2^m) A))^(2^m)
    {
      m += 1;
      v /= 2;
    }
    for (size_t i = 1; i < vPowGen_.size(); i++)
    {
      s *= v / static_cast<double>(i);
      MatrixTools::add(pijt_, s, vPowGen_[i]);
    }
    while (m > 0)  // recover the 2^m
    {
      MatrixTools::mult(pijt_, pijt_, tmpMat_);
      MatrixTools::copy(tmpMat_, pijt_);
      m--;
    }
  }
//  MatrixTools::print(pijt_);

  // Check to avoid numerical issues
  if (t<= NumConstants::SMALL())
    for (size_t i = 0; i < size_; i++)
      for (size_t j = 0; j < size_; j++)
        if (pijt_(i,j)<0.)
          pijt_(i,j)=0.;

  return pijt_;
}

const Matrix<double>& AbstractSubstitutionModel::getdPij_dt(double t) const
{
  if (isNonSingular_)
  {
    if (isDiagonalizable_)
    {
      MatrixTools::mult(rightEigenVectors_, rate_ * eigenValues_ * VectorTools::exp(eigenValues_ * (rate_ * t)), leftEigenVectors_, dpijt_);
    }
    else
    {
      std::vector<double> vdia(size_);
      std::vector<double> vup(size_ - 1);
      std::vector<double> vlo(size_ - 1);
      double c, s, e;
      double l = rate_ * t;
      for (size_t i = 0; i < size_; i++)
      {
        e = std::exp(eigenValues_[i] * l);
        if (iEigenValues_[i] != 0)
        {
          s = std::sin(iEigenValues_[i] * l);
          c = std::cos(iEigenValues_[i] * l);
          vdia[i] = rate_ * (eigenValues_[i] * c - iEigenValues_[i] * s) * e;
          vup[i] = rate_ * (eigenValues_[i] * s + iEigenValues_[i] * c) * e;
          vlo[i] = -vup[i];
          vdia[i + 1] = vdia[i]; // trick to avoid computation
          i++;
        }
        else
        {
          if (i < size_ - 1)
          {
            vup[i] = 0;
            vlo[i] = 0;
          }
        }
      }
      MatrixTools::mult<double>(rightEigenVectors_, vdia, vup, vlo, leftEigenVectors_, dpijt_);
    }
  }
  else
  {
    MatrixTools::getId(size_, dpijt_);
    double s = 1.0;
    double v = rate_ * t;
    size_t m = 0;
    while (v > 0.5)    // r*A*exp(t*r*A)=r*A*(exp(r*t/(2^m) A))^(2^m)
    {
      m += 1;
      v /= 2;
    }
    for (size_t i = 1; i < vPowGen_.size(); i++)
    {
      s *= v / static_cast<double>(i);
      MatrixTools::add(dpijt_, s, vPowGen_[i]);
    }
    while (m > 0)  // recover the 2^m
    {
      MatrixTools::mult(dpijt_, dpijt_, tmpMat_);
      MatrixTools::copy(tmpMat_, dpijt_);
      m--;
    }
    MatrixTools::scale(dpijt_, rate_);
    MatrixTools::mult(vPowGen_[1], dpijt_, tmpMat_);
    MatrixTools::copy(tmpMat_, dpijt_);
  }
  return dpijt_;
}

const Matrix<double>& AbstractSubstitutionModel::getd2Pij_dt2(double t) const
{
  if (isNonSingular_)
  {
    if (isDiagonalizable_)
    {
      MatrixTools::mult(rightEigenVectors_, VectorTools::sqr(rate_ * eigenValues_) * VectorTools::exp(eigenValues_ * (rate_ * t)), leftEigenVectors_, d2pijt_);
    }
    else
    {
      std::vector<double> vdia(size_);
      std::vector<double> vup(size_ - 1);
      std::vector<double> vlo(size_ - 1);
      double c, s, e;
      double l = rate_ * t;
      for (size_t i = 0; i < size_; i++)
      {
        e = std::exp(eigenValues_[i] * l);
        if (iEigenValues_[i] != 0)
        {
          s = std::sin(iEigenValues_[i] * l);
          c = std::cos(iEigenValues_[i] * l);
          vdia[i] = NumTools::sqr(rate_)
                    * ((NumTools::sqr(eigenValues_[i]) - NumTools::sqr(iEigenValues_[i])) * c
                       - 2 * eigenValues_[i] * iEigenValues_[i] * s) * e;
          vup[i] = NumTools::sqr(rate_)
                   * ((NumTools::sqr(eigenValues_[i]) - NumTools::sqr(iEigenValues_[i])) * s
                      - 2 * eigenValues_[i] * iEigenValues_[i] * c) * e;
          vlo[i] = -vup[i];
          vdia[i + 1] = vdia[i]; // trick to avoid computation
          i++;
        }
        else
        {
          if (i < size_ - 1)
          {
            vup[i] = 0;
            vlo[i] = 0;
          }
        }
      }
      MatrixTools::mult<double>(rightEigenVectors_, vdia, vup, vlo, leftEigenVectors_, d2pijt_);
    }
  }
  else
  {
    MatrixTools::getId(size_, d2pijt_);
    double s = 1.0;
    double v = rate_ * t;
    size_t m = 0;
    while (v > 0.5)    // r^2*A^2*exp(t*r*A)=r^2*A^2*(exp(r*t/(2^m) A))^(2^m)
    {
      m += 1;
      v /= 2;
    }
    for (size_t i = 1; i < vPowGen_.size(); i++)
    {
      s *= v / static_cast<double>(i);
      MatrixTools::add(d2pijt_, s, vPowGen_[i]);
    }
    while (m > 0)  // recover the 2^m
    {
      MatrixTools::mult(d2pijt_, d2pijt_, tmpMat_);
      MatrixTools::copy(tmpMat_, d2pijt_);
      m--;
    }
    MatrixTools::scale(d2pijt_, rate_ * rate_);
    MatrixTools::mult(vPowGen_[2], d2pijt_, tmpMat_);
    MatrixTools::copy(tmpMat_, d2pijt_);
  }
  return d2pijt_;
}

/******************************************************************************/

double AbstractSubstitutionModel::getInitValue(size_t i, int state) const throw (IndexOutOfBoundsException, BadIntException)
{
  if (i >= size_)
    throw IndexOutOfBoundsException("AbstractSubstitutionModel::getInitValue", i, 0, size_ - 1);
  if (state < 0 || !alphabet_->isIntInAlphabet(state))
    throw BadIntException(state, "AbstractSubstitutionModel::getInitValue. Character " + alphabet_->intToChar(state) + " is not allowed in model.");
  vector<int> states = alphabet_->getAlias(state);

  for (size_t j = 0; j < states.size(); j++)
  {
     if (getAlphabetStateAsInt(i) == states[j])
      return 1.;
  }
  return 0.;
}

/******************************************************************************/

void AbstractSubstitutionModel::setFreqFromData(const SequenceContainer& data, double pseudoCount)
{
  map<int, int> counts;
  SequenceContainerTools::getCounts(data, counts);
  double t = 0;
  map<int, double> freqs;

  for (int i = 0; i < static_cast<int>(size_); i++)
  {
    t += (counts[i] + pseudoCount);
  }
  for (int i = 0; i < static_cast<int>(size_); i++)
  {
    freqs[i] = (static_cast<double>(counts[i]) + pseudoCount) / t;
  }

  // Re-compute generator and eigen values:
  setFreq(freqs);
}

/******************************************************************************/

void AbstractSubstitutionModel::setFreq(map<int, double>& freqs)
{
  for (size_t i = 0; i < size_; ++i)
  {
    freq_[i] = freqs[static_cast<int>(i)];
  }
  // Re-compute generator and eigen values:
  updateMatrices();
}

/******************************************************************************/

double AbstractSubstitutionModel::getScale() const
{
  vector<double> v;
  MatrixTools::diag(generator_, v);
  return -VectorTools::scalar<double, double>(v, freq_);
}

/******************************************************************************/

void AbstractSubstitutionModel::setScale(double scale) {
  MatrixTools::scale(generator_, scale);
}

/******************************************************************************/

double AbstractSubstitutionModel::getRate() const
{
  return rate_;
}

/******************************************************************************/

void AbstractSubstitutionModel::setRate(double rate)
{
  if (rate <= 0)
    throw Exception("Bad value for rate: " + TextTools::toString(rate));
  
  if (hasParameter("rate"))
    setParameterValue("rate", rate);
  else
    rate_ = rate;
}

void AbstractSubstitutionModel::addRateParameter()
{
  addParameter_(new Parameter(getNamespace() + "rate", rate_, &Parameter::R_PLUS_STAR));
}

/******************************************************************************/

void AbstractReversibleSubstitutionModel::updateMatrices()
{
  RowMatrix<double> Pi;
  MatrixTools::diag(freq_, Pi);
  MatrixTools::mult(exchangeability_, Pi, generator_); // Diagonal elements of the exchangability matrix will be ignored.
  // Compute diagonal elements of the generator:
  for (size_t i = 0; i < size_; i++)
  {
    double lambda = 0;
    for (size_t j = 0; j < size_; j++)
    {
      if (j != i)
        lambda += generator_(i, j);
    }
    generator_(i, i) = -lambda;
  }
  // Normalization:
  double scale = getScale();
  MatrixTools::scale(generator_, 1. / scale);

  // Normalize exchangeability matrix too:
  MatrixTools::scale(exchangeability_, 1. / scale);
  // Compute diagonal elements of the exchangeability matrix:
  for (size_t i = 0; i < size_; i++)
  {
    exchangeability_(i, i) = generator_(i, i) / freq_[i];
  }
  AbstractSubstitutionModel::updateMatrices();
}

/******************************************************************************/

