//
// File: L95.cpp
// Created by: Julien Dutheil
// Created on: Tue Nov 4 11:46 2008
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

#include "L95.h"
#include "../FrequenciesSet/NucleotideFrequenciesSet.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/
 
L95::L95(
	const NucleicAlphabet* alphabet,
	double alpha, double beta, double gamma, double kappa, double theta):
    AbstractParameterAliasable("L95."),
    AbstractNucleotideSubstitutionModel(alphabet, new CanonicalStateMap(alphabet, false), "L95."), alpha_(alpha), beta_(beta), gamma_(gamma), kappa_(kappa), theta_(theta)
{

  addParameter_(new Parameter("L95.alpha", alpha, &Parameter::PROP_CONSTRAINT_IN));
  addParameter_(new Parameter("L95.beta", beta, &Parameter::PROP_CONSTRAINT_IN));
  addParameter_(new Parameter("L95.gamma", gamma, &Parameter::PROP_CONSTRAINT_IN));
  addParameter_(new Parameter("L95.kappa", kappa, new IntervalConstraint(0, 1000, false, false, NumConstants::MILLI()), true));
  addParameter_(new Parameter("L95.theta", theta, new IntervalConstraint(0, 1, false, false, NumConstants::MILLI()), true));

  updateMatrices();
}

/******************************************************************************/
	
void L95::updateMatrices()
{
  alpha_  = getParameterValue("alpha");
  beta_  = getParameterValue("beta");
  gamma_  = getParameterValue("gamma");
  kappa_  = getParameterValue("kappa");
  theta_  = getParameterValue("theta");
  
  freq_[0] = (1-theta_)/2;
  freq_[1] = theta_/2;
  freq_[2] = theta_/2;
  freq_[3] = (1-theta_)/2;
	
  // Generator matrix:
  generator_(0,0) = -kappa_ * theta_ - gamma_;
  generator_(0,1) = kappa_* beta_ * theta_;
  generator_(0,2) = kappa_ * (1-beta_) * theta_;
  generator_(0,3) = gamma_;
  generator_(1,0) = kappa_ * alpha_ * ( 1- theta_);
  generator_(1,1) = -kappa_ * (1 - theta_) + gamma_ - 1;
  generator_(1,2) = 1 - gamma_;
  generator_(1,3) = kappa_ * (1 - theta_) * (1 - alpha_);
  generator_(2,0) = kappa_ * (1 - theta_) * (1 - alpha_);
  generator_(2,1) = 1 - gamma_;
  generator_(2,2) = -kappa_ * (1 - theta_) + gamma_ - 1;
  generator_(2,3) = kappa_ * alpha_ * (1 - theta_);
  generator_(3,0) = gamma_;
  generator_(3,1) = kappa_ * (1-beta_) * theta_;
  generator_(3,2) = kappa_* beta_ * theta_;
  generator_(3,3) = -kappa_ * theta_ - gamma_;

  MatrixTools::scale(generator_, 1. / (2*kappa_*theta_*(1-theta_)+gamma_+theta_-2*theta_*gamma_));
  AbstractSubstitutionModel::updateMatrices();
}

/******************************************************************************/

void L95::setFreq(map<int, double>& freqs)
{
  setParameterValue("theta",freqs[1]+freqs[2]);
  updateMatrices();
}

/******************************************************************************/

