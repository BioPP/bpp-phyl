//
// File: MutationProcess.cpp
// Created by: Julien Dutheil
// Created on: Wed Mar 12 16:11:44 2003
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

#include "MutationProcess.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

using namespace bpp;

/******************************************************************************/

size_t AbstractMutationProcess::mutate(size_t state) const
{
  double alea = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
  for (size_t j = 0; j < size_; j++)
  {
    if (alea < repartition_[state][j]) return j;
  }
  throw Exception("AbstractMutationProcess::mutate. Repartition function is incomplete for state " + TextTools::toString(state));
}

/******************************************************************************/

size_t AbstractMutationProcess::mutate(size_t state, unsigned int n) const
{
  size_t s = state;
  for (unsigned int k = 0; k < n; k++)
  {
    double alea = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
    for (size_t j = 1; j < size_ + 1; j++)
    {
      if (alea < repartition_[s][j])
      {
        s = j;
        break;
      }
    }
  }
  return s;
}

/******************************************************************************/

double AbstractMutationProcess::getTimeBeforeNextMutationEvent(size_t state) const
{
  return RandomTools::randExponential(-1. / model_->Qij(state, state));
}

/******************************************************************************/

size_t AbstractMutationProcess::evolve(size_t initialState, double time) const
{
  double t = 0;
  size_t currentState = initialState;
  t += getTimeBeforeNextMutationEvent(currentState);
  while (t < time)
  {
    currentState = mutate(currentState);
    t += getTimeBeforeNextMutationEvent(currentState);
  }
  return currentState;
}

/******************************************************************************/

MutationPath AbstractMutationProcess::detailedEvolve(size_t initialState, double time) const
{
  MutationPath mp(model_->getAlphabet(), initialState, time);
  double t = 0;
  size_t currentState = initialState;
  t += getTimeBeforeNextMutationEvent(currentState);
  while (t < time)
  {
    currentState = mutate(currentState);
    mp.addEvent(currentState, t);
    t += getTimeBeforeNextMutationEvent(currentState);
  }
  return mp;
}

/******************************************************************************/

SimpleMutationProcess::SimpleMutationProcess(const SubstitutionModel* model) :
  AbstractMutationProcess(model)
{
  size_ = model->getNumberOfStates();
  repartition_ = VVdouble(size_);
  // Each element contains the probabilities concerning each character in the alphabet.

  // We will now initiate each of these probability vector.
  RowMatrix<double> Q = model->getGenerator();
  for (size_t i = 0; i < size_; i++)
  {
    repartition_[i] = Vdouble(size_);
    double cum = 0;
    double sum_Q = 0;
    for (size_t j = 0; j < size_; j++)
    {
      if (j != i) sum_Q += Q(i, j);
    }
    for (size_t j = 0; j < size_; j++)
    {
      if (j != i)
      {
        cum += model->Qij(i, j) / sum_Q;
        repartition_[i][j] = cum;
      }
      else repartition_[i][j] = -1;
      // Forbiden value: does not correspond to a change.
    }
  }
  // Note that I use cumulative probabilities in repartition_ (hence the name).
  // These cumulative probabilities are useful for the 'mutate(...)' function.
}

SimpleMutationProcess::~SimpleMutationProcess() {}

/******************************************************************************/

size_t SimpleMutationProcess::evolve(size_t initialState, double time) const
{
  // Compute all cumulative pijt:
  Vdouble pijt(size_);
  pijt[0] = model_->Pij_t(initialState, 0, time);
  for (size_t i = 1; i < size_; i++)
  {
    pijt[i] = pijt[i - 1] + model_->Pij_t(initialState, i, time);
  }
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
  for (size_t i = 0; i < size_; i++)
  {
    if (rand < pijt[i]) return i;
  }
  throw Exception("SimpleSimulationProcess::evolve(intialState, time): error all pijt do not sum to one (total sum = " + TextTools::toString(pijt[size_ - 1]) + ").");
}

/******************************************************************************/

SelfMutationProcess::SelfMutationProcess(size_t alphabetSize) :
  AbstractMutationProcess(0)
{
  size_ = alphabetSize;
  repartition_ = VVdouble(size_);
  // Each element contains the probabilities concerning each character in the alphabet.

  // We will now initiate each of these probability vector.
  for (size_t i = 0; i < size_; i++)
  {
    repartition_[i] = Vdouble(size_);
    for (size_t j = 0; j < size_; j++)
    {
      repartition_[i][j] = static_cast<double>(j + 1) / static_cast<double>(size_);
    }
  }
  // Note that I use cumulative probabilities in repartition_ (hence the name).
  // These cumulative probabilities are useful for the 'mutate(...)' function.
}

SelfMutationProcess::~SelfMutationProcess() {}

/******************************************************************************/

