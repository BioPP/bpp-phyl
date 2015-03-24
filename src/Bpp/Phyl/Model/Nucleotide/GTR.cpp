//
// File: GTR.cpp
// Created by: Julien Dutheil
// Created on: Tue Oct 25 10:21 2005
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

#include "GTR.h"
#include "../FrequenciesSet/NucleotideFrequenciesSet.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/

GTR::GTR(
    const NucleicAlphabet* alpha,
    double a,
    double b,
    double c,
    double d,
    double e,
    double piA,
    double piC,
    double piG,
    double piT) :
  AbstractParameterAliasable("GTR."),
  AbstractReversibleNucleotideSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "GTR."),
  a_(a), b_(b), c_(c), d_(d), e_(e), piA_(piA), piC_(piC), piG_(piG), piT_(piT), theta_(piG + piC), theta1_(piA / (1. - theta_)), theta2_(piG / theta_), p_()
{
  addParameter_(new Parameter("GTR.a", a, &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("GTR.b", b, &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("GTR.c", c, &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("GTR.d", d, &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("GTR.e", e, &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("GTR.theta", theta_, &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter("GTR.theta1", theta1_, &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter("GTR.theta2", theta2_, &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  updateMatrices();
}

/******************************************************************************/
  
void GTR::updateMatrices()
{
  a_ = getParameterValue("a");
  b_ = getParameterValue("b");
  c_ = getParameterValue("c");
  d_ = getParameterValue("d");
  e_ = getParameterValue("e");
  theta_  = getParameterValue("theta");
  theta1_ = getParameterValue("theta1");
  theta2_ = getParameterValue("theta2");
  piA_ = theta1_ * (1. - theta_);
  piC_ = (1. - theta2_) * theta_;
  piG_ = theta2_ * theta_;
  piT_ = (1. - theta1_) * (1. - theta_);
  p_ = 2*(a_*piC_*piT_+b_*piA_*piT_+c_*piG_*piT_+d_*piA_*piC_+e_*piC_*piG_+piA_*piG_);

  freq_[0] = piA_;
  freq_[1] = piC_;
  freq_[2] = piG_;
  freq_[3] = piT_;
  
  // Exchangeability matrix:
  exchangeability_(0,0) = (-b_*piT_-piG_-d_*piC_)/(piA_ * p_);
  exchangeability_(1,0) = d_/p_;
  exchangeability_(0,1) = d_/p_;
  exchangeability_(2,0) = 1/p_;
  exchangeability_(0,2) = 1/p_;
  exchangeability_(3,0) = b_/p_;
  exchangeability_(0,3) = b_/p_;
  exchangeability_(1,1) = (-a_*piT_-e_*piG_-d_*piA_)/(piC_ * p_);
  exchangeability_(1,2) = e_/p_;
  exchangeability_(2,1) = e_/p_;
  exchangeability_(1,3) = a_/p_;
  exchangeability_(3,1) = a_/p_;
  exchangeability_(2,2) = (-c_*piT_-e_*piC_-piA_)/(piG_ * p_);
  exchangeability_(2,3) = c_/p_;
  exchangeability_(3,2) = c_/p_;
  exchangeability_(3,3) = (-c_*piG_-a_*piC_-b_*piA_)/(piT_ * p_);

  AbstractReversibleSubstitutionModel::updateMatrices();
}

/******************************************************************************/

void GTR::setFreq(map<int, double>& freqs)
{
  piA_ = freqs[0];
  piC_ = freqs[1];
  piG_ = freqs[2];
  piT_ = freqs[3];
  vector<string> thetas(3);
  thetas[0] = getNamespace() + "theta";
  thetas[1] = getNamespace() + "theta1";
  thetas[2] = getNamespace() + "theta2";
  ParameterList pl = getParameters().subList(thetas);
  pl[0].setValue(piC_ + piG_);
  pl[1].setValue(piA_ / (piA_ + piT_));
  pl[2].setValue(piG_ / (piC_ + piG_));
  setParametersValues(pl);
}

/******************************************************************************/

