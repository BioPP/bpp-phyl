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

AbstractTransitionModel::AbstractTransitionModel(const Alphabet* alpha, std::shared_ptr<const StateMap> stateMap, const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  alphabet_(alpha),
  stateMap_(stateMap),
  size_(alpha->getSize()),
  rate_(1),
  freq_(size_),
  pijt_(size_, size_),
  dpijt_(size_, size_),
  d2pijt_(size_, size_)
{
  if (computeFrequencies())
    for (auto& fr : freq_)
      fr = 1.0 / static_cast<double>(size_);
}

/******************************************************************************/

  double AbstractTransitionModel::getRate() const
  {
    return rate_;
  }

/******************************************************************************/

void AbstractTransitionModel::setRate(double rate)
{
  if (rate <= 0)
    throw Exception("Bad value for rate: " + TextTools::toString(rate));
  
  if (hasParameter("rate"))
    setParameterValue("rate", rate);
  else
    rate_ = rate;
}

void AbstractTransitionModel::addRateParameter()
{
  addParameter_(new Parameter(getNamespace() + "rate", rate_, Parameter::R_PLUS_STAR));
}

/******************************************************************************/
double AbstractTransitionModel::getInitValue(size_t i, int state) const
{
  if (i >= size_)
    throw IndexOutOfBoundsException("AbstractTransitionModel::getInitValue", i, 0, size_ - 1);
  if (state < 0 || !alphabet_->isIntInAlphabet(state))
    throw BadIntException(state, "AbstractTransitionModel::getInitValue. Character " + alphabet_->intToChar(state) + " is not allowed in model.");
  vector<int> states = alphabet_->getAlias(state);

  for (size_t j = 0; j < states.size(); j++)
  {
     if (getAlphabetStateAsInt(i) == states[j])
      return 1.;
  }
  return 0.;
}

/******************************************************************************/

void AbstractTransitionModel::setFreqFromData(const SequencedValuesContainer& data, double pseudoCount)
{
  map<int, double> freqs;
  SequenceContainerTools::getFrequencies(data, freqs, pseudoCount);
  // Re-compute generator and eigen values:
  setFreq(freqs);
}

/******************************************************************************/

void AbstractTransitionModel::setFreq(map<int, double>& freqs)
{
  for (auto i : freqs)
    freq_[(size_t)i.first]=i.second;

  // Re-compute generator and eigen values:
  updateMatrices();
}



/******************************************************************************/

AbstractSubstitutionModel::AbstractSubstitutionModel(const Alphabet* alpha, std::shared_ptr<const StateMap> stateMap, const std::string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractTransitionModel(alpha, stateMap, prefix),
  isScalable_(true),
  generator_(size_, size_),
  computeFreq_(true),
  exchangeability_(size_, size_),
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
}


/******************************************************************************/

void AbstractSubstitutionModel::updateMatrices()
{

  // Compute eigen values and vectors:
  if (enableEigenDecomposition())
  {
    // Look for null lines (such as stop lines)
    // ie null diagonal elements

    size_t nbStop=0;
    size_t salph = getNumberOfStates();
    vector<bool> vnull(salph); // vector of the indices of lines with
                               // only zeros

    for (size_t i = 0; i < salph; i++)
    {
      bool flag=(abs(generator_(i, i)) < NumConstants::TINY());

      if (flag)
        for (size_t j = 0; j < salph; j++)
          if (abs(generator_(j, i)) >= NumConstants::TINY())
          {
            flag=false;
            break;
          }

      if (flag)
      {
        nbStop++;
        vnull[i]=true;
      }
      else
        vnull[i]=false;
    }
        
    if (nbStop != 0)
    {
      size_t salphok=salph - nbStop;
      
      RowMatrix<double> gk(salphok, salphok);
      size_t gi = 0, gj = 0;

      for (size_t i = 0; i < salph; i++)
      {
        if (!vnull[i])
        {
          gj = 0;
          for (size_t j = 0; j < salph; j++)
          {
            if (!vnull[j])
            {
              gk(i - gi, j - gj) = generator_(i, j);
            }
            else
              gj++;
          }
        }
        else
          gi++;
      }

      EigenValue<double> ev(gk);
      eigenValues_ = ev.getRealEigenValues();
      iEigenValues_ = ev.getImagEigenValues();

      for (size_t i = 0; i < nbStop; i++)
      {
        eigenValues_.push_back(0);
        iEigenValues_.push_back(0);
      }

      RowMatrix<double> rev = ev.getV();
      rightEigenVectors_.resize(salph, salph);
      gi = 0;
      for (size_t i = 0; i < salph; i++)
      {
        if (vnull[i])
        {
          gi++;
          for (size_t j = 0; j < salph; j++)
          {
            rightEigenVectors_(i, j) = 0;
          }

          rightEigenVectors_(i, salphok + gi - 1) = 1;
        }
        else
        {
          for (size_t j = 0; j < salphok; j++)
          {
            rightEigenVectors_(i, j) = rev(i - gi, j);
          }

          for (size_t j = salphok; j < salph; j++)
          {
            rightEigenVectors_(i, j) = 0;
          }
        }
      }
    }
    else
    {
      EigenValue<double> ev(generator_);
      rightEigenVectors_ = ev.getV();
      eigenValues_ = ev.getRealEigenValues();
      iEigenValues_ = ev.getImagEigenValues();
      nbStop = 0;
    }

    /// Now check inversion and diagonalization
    try
    {
      MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);

      // is it diagonalizable ?
      isDiagonalizable_ = true;

      if (!dynamic_cast<ReversibleSubstitutionModel*>(this))
      {
        for (auto& vi : iEigenValues_)
        {
          if (abs(vi) > NumConstants::TINY())
          {
            isDiagonalizable_ = false;
            break;
          }
        }
      }
      
      // looking for the vector of 0 eigenvalues

      vector<size_t> vNullEv;
      double fact=0.1;
      while (vNullEv.size()==0 && fact<1000)
      {
        fact*=10;
        
        for (size_t i = 0; i< salph - nbStop; i++)
          if ((abs(eigenValues_[i]) < fact*NumConstants::SMALL()) && (abs(iEigenValues_[i]) < NumConstants::SMALL()))
            vNullEv.push_back(i);
      }
      

      // pb to find unique null eigenvalue      
      isNonSingular_=(vNullEv.size()==1);
      
      size_t nulleigen;
      
      double val;
      if (!isNonSingular_)
      {
        //look or check which non-stop right eigen vector elements are
        //equal.
        for (auto cnull : vNullEv)
        {
          size_t i = 0;
          while (vnull[i])
            i++;
          
          val = rightEigenVectors_(i, cnull);
          i++;
          
          while (i < salph)
          {
            if (!vnull[i])
            {
              if (abs((rightEigenVectors_(i, cnull) - val)/val) > NumConstants::SMALL())
                break;
            }
            i++;
          }
          
          if (i >= salph)
          {
            isNonSingular_ = true;
            nulleigen=cnull;
            break;
          }
        }
      }
      else
        nulleigen=vNullEv[0];
      
      if (isNonSingular_)
      {
        eigenValues_[nulleigen] = 0; // to avoid approximation errors on long long branches
        iEigenValues_[nulleigen] = 0; // to avoid approximation errors on long long branches

        if (computeFrequencies())
        {
          for (size_t i = 0; i < salph; i++)
            freq_[i] = leftEigenVectors_(nulleigen, i);
        
          double x = VectorTools::sum(freq_);        
          freq_ /= x;
        }
      }
      else
      {
        ApplicationTools::displayMessage("AbstractSubstitutionModel::updateMatrices : Unable to find eigenvector for eigenvalue 0. Taylor series used instead.");
        isDiagonalizable_ = false;
      }
    }
    // if rightEigenVectors_ is singular
    catch (ZeroDivisionException& e)
    {
      ApplicationTools::displayMessage("AbstractSubstitutionModel::updateMatrices : Singularity during diagonalization. Taylor series used instead.");
      isNonSingular_ = false;
      isDiagonalizable_ = false;
    }

    if (!isNonSingular_)
    {
      double min = generator_(0, 0);
      for (size_t i = 1; i < salph; i++)
      {
        if (min > generator_(i, i))
          min = generator_(i, i);
      }

      setScale(-1 / min);

      if (vPowGen_.size() == 0)
        vPowGen_.resize(30);

      
      if (computeFrequencies())
      {
        MatrixTools::getId(salph, tmpMat_);    // to compute the equilibrium frequency  (Q+Id)^256
        MatrixTools::add(tmpMat_, generator_);
        MatrixTools::pow(tmpMat_, 256, vPowGen_[0]);

        for (size_t i = 0; i < salph; i++)
          freq_[i] = vPowGen_[0](0, i);
      }

      MatrixTools::getId(salph, vPowGen_[0]);
    }

    // normalization
    normalize();
    
    if (!isNonSingular_)
      MatrixTools::Taylor(generator_, 30, vPowGen_);
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
  // if (t<= NumConstants::SMALL())
    for (size_t i = 0; i < size_; i++)
      for (size_t j = 0; j < size_; j++)
        if (pijt_(i,j)<0.)
        {
          if (std::abs(pijt_(i,j))>NumConstants::SMALL())
          {
            throw Exception("There is an issue in the computation of transition matrix of " + getName() + " : pijt_(" + to_string(i) + "," + to_string(j) + ", " + to_string(t) + ")=" + to_string(pijt_(i,j)));
          }
          pijt_(i,j)=0.;
          
        }
  return pijt_;
}

/******************************************************************************/

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

/******************************************************************************/

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

double AbstractSubstitutionModel::getScale() const
{
  vector<double> v;
  MatrixTools::diag(generator_, v);
  return -VectorTools::scalar<double, double>(v, freq_);
}

/******************************************************************************/

void AbstractSubstitutionModel::setScale(double scale)
{
  if (isScalable_)
  {
    MatrixTools::scale(generator_, scale);
    eigenValues_ *= scale;
    iEigenValues_ *= scale;
  }
}


/******************************************************************************/

void AbstractSubstitutionModel::setDiagonal()
{
  for (size_t i = 0; i < size_; i++)
  {
    double lambda=0;
    Vdouble& row=generator_.getRow(i);
    
    for (size_t j = 0; j < size_; j++)
    {
      if (j != i)
        lambda += row[j];
    }
    row[i] = -lambda;
  }
}


/******************************************************************************/

void AbstractSubstitutionModel::normalize() 
{
  if (isScalable_)
    setScale(1/getScale());
}

/******************************************************************************/

void AbstractReversibleSubstitutionModel::updateMatrices()
{
  MatrixTools::hadamardMult(exchangeability_, freq_, generator_, false); // Diagonal elements of the exchangeability matrix will be ignored.

  // Normalization:
  setDiagonal();
  normalize();
  
  AbstractSubstitutionModel::updateMatrices();
}

/******************************************************************************/

