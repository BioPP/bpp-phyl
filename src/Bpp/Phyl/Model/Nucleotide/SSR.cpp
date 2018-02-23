//
// File: SSR.cpp
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

#include "SSR.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/
 
SSR::SSR(
  const NucleicAlphabet* alpha,
  double beta,
  double gamma,
  double delta,
  double theta):
  AbstractParameterAliasable("SSR."),
  AbstractReversibleNucleotideSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "SSR."),
  beta_(beta), gamma_(gamma), delta_(delta), theta_(theta),
  piA_((1. - theta) / 2.), piC_(theta / 2.), piG_(theta / 2.), piT_((1. - theta) / 2.)
{
  addParameter_(new Parameter("SSR.beta" , beta , &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("SSR.gamma", gamma, &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("SSR.delta", delta, &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("SSR.theta" , theta , &Parameter::PROP_CONSTRAINT_EX));
  updateMatrices();
}

/******************************************************************************/
	
void SSR::updateMatrices()
{
  beta_  = getParameterValue("beta");
  gamma_ = getParameterValue("gamma");
  delta_ = getParameterValue("delta");
  theta_ = getParameterValue("theta");
  
  freq_[0] = piA_ = (1. - theta_)/2.;
  freq_[1] = piC_ = theta_/2.;
  freq_[2] = piG_ = theta_/2;
  freq_[3] = piT_ = (1. - theta_)/2.;
	
  // Exchangeability matrix:
  exchangeability_(0,0) = -gamma_*piT_-piG_-beta_*piC_;
  exchangeability_(1,0) = beta_;
  exchangeability_(0,1) = beta_;
  exchangeability_(2,0) = 1.;
  exchangeability_(0,2) = 1.;
  exchangeability_(3,0) = gamma_;
  exchangeability_(0,3) = gamma_;
  exchangeability_(1,1) = -piT_-delta_*piG_-beta_*piA_;
  exchangeability_(1,2) = delta_;
  exchangeability_(2,1) = delta_;
  exchangeability_(1,3) = 1.;
  exchangeability_(3,1) = 1.;
  exchangeability_(2,2) = -beta_*piT_-delta_*piC_-piA_;
  exchangeability_(2,3) = beta_;
  exchangeability_(3,2) = beta_;
  exchangeability_(3,3) = -beta_*piG_-piC_-gamma_*piA_;
  
  AbstractReversibleSubstitutionModel::updateMatrices();
}

/******************************************************************************/

void SSR::setFreq(map<int, double>& freqs)
{
  piC_ = freqs[1];
  piG_ = freqs[2];
  setParameterValue("theta",piC_ + piG_);
  updateMatrices();
}

/******************************************************************************/

