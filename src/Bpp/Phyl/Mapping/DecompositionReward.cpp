//
// File: DecompositionReward.h
// Created by: Laurent Guéguen
// Created on: mercredi 27 mars 2013, à 12h 36
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

#include "DecompositionReward.h"

#include "Bpp/Numeric/Matrix/MatrixTools.h"
#include <vector>
#include <typeinfo>

using namespace bpp;
using namespace std;

/******************************************************************************/

DecompositionReward::DecompositionReward(const SubstitutionModel* model, AlphabetIndex1* alphIndex) :
  AbstractReward(alphIndex),
  DecompositionMethods(model),
  rewards_(nbStates_, nbStates_),
  currentLength_(-1.)
{
  //Check compatiblity between model and alphabet Index:
  if (typeid(model->getAlphabet()) != typeid(alphIndex_->getAlphabet()))
    throw Exception("DecompositionReward (constructor): alphabets do not match between alphabet index and model.");

  //Initialize the B matrice. This is done once for all,
  //unless the number of states changes:

  initBMatrices_();
  initRewards_();
  
  fillBMatrice_();
  computeProducts_();
}				

DecompositionReward::DecompositionReward(const StateMap& statemap, AlphabetIndex1* alphIndex) :
  AbstractReward(alphIndex),
  DecompositionMethods(statemap),
  rewards_(nbStates_, nbStates_),
  currentLength_(-1.)
{
}				

/******************************************************************************/

void DecompositionReward::initRewards_()
{
  rewards_.resize(nbStates_, nbStates_);
}

/******************************************************************************/

void DecompositionReward::fillBMatrice_()
{
  vector<int> supportedStates = model_->getAlphabetStates();
  for (size_t j = 0; j < nbStates_; ++j) 
    bMatrices_[0](j, j) = getAlphabetIndex()->getIndex(supportedStates[j]);
}

/******************************************************************************/

void DecompositionReward::computeRewards_(double length) const
{
  computeExpectations(rewards_, length);

  // Now we must divide by pijt:
  RowMatrix<double> P = model_->getPij_t(length);
  for (size_t j = 0; j < nbStates_; j++) {
    for (size_t k = 0; k < nbStates_; k++) {
      rewards_(j, k) /= P(j, k);
      if (std::isnan(rewards_(j, k)) || std::isnan(-rewards_(j, k)) || std::isinf(rewards_(j, k)))
        rewards_(j, k) = 0.;
    }
  }
}

/******************************************************************************/

Matrix<double>* DecompositionReward::getAllRewards(double length) const
{
  if (!model_)
    throw Exception("DecompositionReward::getAllRewards: model not defined.");
  
  if (length < 0)
    throw Exception("DecompositionReward::getAllRewards. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeRewards_(length);
    currentLength_ = length;
  }
  return new RowMatrix<double>(rewards_);
}

/******************************************************************************/

void DecompositionReward::storeAllRewards(double length, Eigen::MatrixXd& mat) const
{
  if (!model_)
    throw Exception("DecompositionReward::storeAllRewards: model not defined.");

  if (length < 0)
    throw Exception("DecompositionReward::storeAllRewards. Negative branch length: " + TextTools::toString(length) + ".");
  
  if (length != currentLength_)
  {
    computeRewards_(length);
    currentLength_ = length;
  }

  mat.resize(Eigen::Index(nbStates_), Eigen::Index(nbStates_));

  for (size_t j = 0; j < nbStates_; j++) 
    for (size_t k = 0; k < nbStates_; k++) 
      mat(Eigen::Index(j), Eigen::Index(k)) = rewards_(j, k);
}

/******************************************************************************/

double DecompositionReward::getReward(size_t initialState, size_t finalState, double length) const
{
  if (length < 0)
    throw Exception("DecompositionReward::getRewards. Negative branch length: " + TextTools::toString(length) + ".");
  if (length != currentLength_)
  {
    computeRewards_(length);
    currentLength_ = length;
  }
  return rewards_(initialState, finalState);
}

/******************************************************************************/

void DecompositionReward::setSubstitutionModel(const SubstitutionModel* model)
{
  DecompositionMethods::setSubstitutionModel(model);

  if (!model)
    return;
  
  //Check compatiblity between model and substitution register:
  if (typeid(model->getAlphabet()) != typeid(alphIndex_->getAlphabet()))
    throw Exception("DecompositionReward::setSubstitutionModel: alphabets do not match between alphabet index and model.");


  initRewards_();
  
  fillBMatrice_();  
  computeProducts_();

  //Recompute rewards:

  computeRewards_(currentLength_);
}

/******************************************************************************/

void DecompositionReward::alphabetIndexHasChanged()
{
  if (!model_)
    return;
  
  //Check compatiblity between model and substitution register:
  if (typeid(model_->getAlphabet()) != typeid(alphIndex_->getAlphabet()))
    throw Exception("DecompositionReward::AlphabetIndexHasChanged: alphabets do not match between alphbaet index and model.");

  fillBMatrice_();
  computeProducts_();

  //Recompute rewards:
  if (currentLength_ > 0)
   computeRewards_(currentLength_);
}


