//
// File: TwoParameterBinarySubstitutionModel.cpp
// Created by: Keren Halabi
// Created on: 2018
//

/*
   Copyright or � or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "TwoParameterBinarySubstitutionModel.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace std;

/******************************************************************************/

TwoParameterBinarySubstitutionModel::TwoParameterBinarySubstitutionModel(const BinaryAlphabet* alpha, double mu, double pi0) :
  AbstractParameterAliasable("TwoParameterBinary."),
  AbstractReversibleSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "Binary."),
  mu_(mu),
  pi0_(pi0)
{
  addParameter_(new Parameter(getNamespace() + "mu", mu_, &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter(getNamespace() + "pi0", pi0_, &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  updateMatrices();
}

/******************************************************************************/

void TwoParameterBinarySubstitutionModel::updateMatrices()
{
  mu_ = getParameterValue("mu");
  pi0_ = getParameterValue("pi0");

  // Frequences:
  freq_[0] = pi0_;
  freq_[1] = 1 - pi0_;

  // Exchangeability matrix:
  exchangeability_(0,0) = -1 * mu_ * (1-pi0_) / pi0_;
  exchangeability_(0,1) = 1;
  exchangeability_(1,0) = 1;
  exchangeability_(1,1) = -1 * mu_ * pi0_ / (1-pi0_);

  AbstractReversibleSubstitutionModel::updateMatrices();
}

/******************************************************************************/

void TwoParameterBinarySubstitutionModel::setMuBounds(double lb, double ub)
{
  IntervalConstraint* bounds = new IntervalConstraint(lb, ub, true, true); 
  getParameter_("mu").setConstraint(bounds, true);
}