//
// File: F84.cpp
// Created by: Julien Dutheil
// Created on: Tue May 23 11:13 2005
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

#include "F84.h"
#include "../FrequenciesSet/NucleotideFrequenciesSet.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/

F84::F84(
    const NucleicAlphabet* alpha,
    double kappa,
    double piA,
    double piC,
    double piG,
    double piT) :
  AbstractParameterAliasable("F84."),
  AbstractReversibleNucleotideSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "F84."),
  kappa_(kappa), piA_(piA), piC_(piC), piG_(piG), piT_(piT), piY_(), piR_(),
  r_(), k1_(), k2_(), theta_(piG + piC), theta1_(piA / (1. - theta_)), theta2_(piG / theta_),
  l_(), exp1_(), exp2_(), p_(size_, size_)
{
  addParameter_(new Parameter("F84.kappa", kappa, &Parameter::R_PLUS));
  addParameter_(new Parameter("F84.theta" , theta_ , &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter("F84.theta1", theta1_, &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter("F84.theta2", theta2_, &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  updateMatrices();
}

/******************************************************************************/

void F84::updateMatrices()
{
  kappa_ = getParameterValue("kappa");
  theta_  = getParameterValue("theta");
  theta1_ = getParameterValue("theta1");
  theta2_ = getParameterValue("theta2");
  piA_ = theta1_ * (1. - theta_);
  piC_ = (1. - theta2_) * theta_;
  piG_ = theta2_ * theta_;
  piT_ = (1. - theta1_) * (1. - theta_);
  piR_ = piA_ + piG_;
  piY_ = piT_ + piC_;
  r_ = 1. / (1 - piA_ * piA_ - piC_ * piC_ - piG_ * piG_ - piT_*piT_ + 2. * kappa_ * (piC_ * piT_ / piY_ + piA_ * piG_ / piR_));
  k1_ = 1;
  k2_ = kappa_ + 1;
  
  // Frequences:
  freq_[0] = piA_;
  freq_[1] = piC_;
  freq_[2] = piG_;
  freq_[3] = piT_;

  //Generator:
  generator_(0, 0) = -(                  piC_ + (1 + kappa_ / piR_) * piG_ +       piT_);
  generator_(1, 1) = -(      piA_ +                   piG_ + (1 + kappa_ / piY_) * piT_); 
  generator_(2, 2) = -((1 + kappa_ / piR_) * piA_ +       piC_             +       piT_);
  generator_(3, 3) = -(      piA_ + (1 + kappa_ / piY_) * piC_ +       piG_            );

  generator_(1, 0) = piA_;
  generator_(3, 0) = piA_;
  generator_(0, 1) = piC_;
  generator_(2, 1) = piC_;
  generator_(1, 2) = piG_;
  generator_(3, 2) = piG_;
  generator_(0, 3) = piT_;
  generator_(2, 3) = piT_;
  
  generator_(2, 0) = (1 + kappa_ / piR_) * piA_;
  generator_(3, 1) = (1 + kappa_ / piY_) * piC_;
  generator_(0, 2) = (1 + kappa_ / piR_) * piG_;
  generator_(1, 3) = (1 + kappa_ / piY_) * piT_;
  
  // Normalization:
  MatrixTools::scale(generator_, r_);
  
  // Exchangeability:
  exchangeability_(0,0) = generator_(0,0) / piA_;
  exchangeability_(0,1) = generator_(0,1) / piC_; 
  exchangeability_(0,2) = generator_(0,2) / piG_; 
  exchangeability_(0,3) = generator_(0,3) / piT_;

  exchangeability_(1,0) = generator_(1,0) / piA_; 
  exchangeability_(1,1) = generator_(1,1) / piC_; 
  exchangeability_(1,2) = generator_(1,2) / piG_; 
  exchangeability_(1,3) = generator_(1,3) / piT_; 
  
  exchangeability_(2,0) = generator_(2,0) / piA_; 
  exchangeability_(2,1) = generator_(2,1) / piC_; 
  exchangeability_(2,2) = generator_(2,2) / piG_; 
  exchangeability_(2,3) = generator_(2,3) / piT_; 
  
  exchangeability_(3,0) = generator_(3,0) / piA_;
  exchangeability_(3,1) = generator_(3,1) / piC_; 
  exchangeability_(3,2) = generator_(3,2) / piG_; 
  exchangeability_(3,3) = generator_(3,3) / piT_;

  // Eigen values:
  eigenValues_[0] = 0;
  eigenValues_[1] = -r_ * (1+kappa_);
  eigenValues_[2] = -r_ * (1+kappa_); 
  eigenValues_[3] = -r_;
  
  // Eigen vectors:
  leftEigenVectors_(0,0) = piA_;
  leftEigenVectors_(0,1) = piC_;
  leftEigenVectors_(0,2) = piG_;
  leftEigenVectors_(0,3) = piT_;

  leftEigenVectors_(1,0) = 0.;
  leftEigenVectors_(1,1) = piT_ / piY_;
  leftEigenVectors_(1,2) = 0.;
  leftEigenVectors_(1,3) = -piT_ / piY_;

  leftEigenVectors_(2,0) = piG_ / piR_;
  leftEigenVectors_(2,1) = 0.;
  leftEigenVectors_(2,2) = -piG_ / piR_;
  leftEigenVectors_(2,3) = 0.;

  leftEigenVectors_(3,0) = piA_*piY_ / piR_;
  leftEigenVectors_(3,1) = -piC_;
  leftEigenVectors_(3,2) = piG_*piY_ / piR_;
  leftEigenVectors_(3,3) = -piT_;

  rightEigenVectors_(0,0) = 1.;
  rightEigenVectors_(0,1) = 0.;
  rightEigenVectors_(0,2) = 1.;
  rightEigenVectors_(0,3) = 1.;
  
  rightEigenVectors_(1,0) = 1.;
  rightEigenVectors_(1,1) = 1.;
  rightEigenVectors_(1,2) = 0.;;
  rightEigenVectors_(1,3) = -piR_ / piY_;

  rightEigenVectors_(2,0) = 1.;
  rightEigenVectors_(2,1) = 0.;
  rightEigenVectors_(2,2) = -piA_ / piG_;
  rightEigenVectors_(2,3) = 1.;

  rightEigenVectors_(3,0) = 1.;
  rightEigenVectors_(3,1) = -piC_ / piT_;
  rightEigenVectors_(3,2) = 0.;
  rightEigenVectors_(3,3) = -piR_ / piY_;
}
  
/******************************************************************************/

double F84::Pij_t(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-k1_*l_);
  exp2_ = exp(-k2_*l_);
  
  switch(i) {
    //A
    case 0 : {
      switch(j) {
        case 0 : return piA_ * (1. + (piY_/piR_) * exp1_) + (piG_/piR_) * exp2_; //A
        case 1 : return piC_ * (1. -               exp1_);                       //C
        case 2 : return piG_ * (1. + (piY_/piR_) * exp1_) - (piG_/piR_) * exp2_; //G
        case 3 : return piT_ * (1. -               exp1_);                       //T, U
      }
    } 
    //C
    case 1 : {
      switch(j) {
        case 0 : return piA_ * (1. -               exp1_);                       //A
        case 1 : return piC_ * (1. + (piR_/piY_) * exp1_) + (piT_/piY_) * exp2_; //C
        case 2 : return piG_ * (1. -               exp1_);                       //G
        case 3 : return piT_ * (1. + (piR_/piY_) * exp1_) - (piT_/piY_) * exp2_; //T, U
      }
    }
    //G
    case 2 : {
      switch(j) {
        case 0 : return piA_ * (1. + (piY_/piR_) * exp1_) - (piA_/piR_) * exp2_; //A
        case 1 : return piC_ * (1. -               exp1_);                       //C
        case 2 : return piG_ * (1. + (piY_/piR_) * exp1_) + (piA_/piR_) * exp2_; //G
        case 3 : return piT_ * (1. -               exp1_);                       //T, U
      }
    }
    //T, U
    case 3 : {
      switch(j) {
        case 0 : return piA_ * (1. -               exp1_);                       //A
        case 1 : return piC_ * (1. + (piR_/piY_) * exp1_) - (piC_/piY_) * exp2_; //C
        case 2 : return piG_ * (1. -               exp1_);                       //G
        case 3 : return piT_ * (1. + (piR_/piY_) * exp1_) + (piC_/piY_) * exp2_; //T, U
      }
    }
  }
  return 0;
}

/******************************************************************************/

double F84::dPij_dt(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-k1_*l_);
  exp2_ = exp(-k2_*l_);
  
  switch(i) {
    //A
    case 0 : {
      switch(j) {
        case 0 : return rate_ * r_ * (piA_ * -(piY_/piR_) * exp1_ - (piG_/piR_) * k2_ * exp2_); //A
        case 1 : return rate_ * r_ * (piC_ *                exp1_);                             //C
        case 2 : return rate_ * r_ * (piG_ * -(piY_/piR_) * exp1_ + (piG_/piR_) * k2_ * exp2_); //G
        case 3 : return rate_ * r_ * (piT_ *                exp1_);                             //T, U
      }
    } 
    //C
    case 1 : {
      switch(j) {
        case 0 : return rate_ * r_ * (piA_ *                exp1_);                             //A
        case 1 : return rate_ * r_ * (piC_ * -(piR_/piY_) * exp1_ - (piT_/piY_) * k2_ * exp2_); //C
        case 2 : return rate_ * r_ * (piG_ *                exp1_);                             //G
        case 3 : return rate_ * r_ * (piT_ * -(piR_/piY_) * exp1_ + (piT_/piY_) * k2_ * exp2_); //T, U
      }
    }
    //G
    case 2 : {
      switch(j) {
        case 0 : return rate_ * r_ * (piA_ * -(piY_/piR_) * exp1_ + (piA_/piR_) * k2_ * exp2_); //A
        case 1 : return rate_ * r_ * (piC_ *                exp1_);                             //C
        case 2 : return rate_ * r_ * (piG_ * -(piY_/piR_) * exp1_ - (piA_/piR_) * k2_ * exp2_); //G
        case 3 : return rate_ * r_ * (piT_ *                exp1_);                             //T, U
      }
    }
    //T, U
    case 3 : {
      switch(j) {
        case 0 : return rate_ * r_ * (piA_ *                exp1_);                             //A
        case 1 : return rate_ * r_ * (piC_ * -(piR_/piY_) * exp1_ + (piC_/piY_) * k2_ * exp2_); //C
        case 2 : return rate_ * r_ * (piG_ *                exp1_);                             //G
        case 3 : return rate_ * r_ * (piT_ * -(piR_/piY_) * exp1_ - (piC_/piY_) * k2_ * exp2_); //T, U
      }
    }
  }
  return 0;
}

/******************************************************************************/

double F84::d2Pij_dt2(size_t i, size_t j, double d) const
{
  double r_2 = rate_ * rate_ * r_ * r_;
  l_ = rate_ * r_ * d;
  double k2_2 = k2_ * k2_;
  exp1_ = exp(-k1_*l_);
  exp2_ = exp(-k2_*l_);
  
  switch(i) {
    //A
    case 0 : {
      switch(j) {
        case 0 : return r_2 * (piA_ * (piY_/piR_) * exp1_ + (piG_/piR_) * k2_2 * exp2_); //A
        case 1 : return r_2 * (piC_ *             - exp1_);                              //C
        case 2 : return r_2 * (piG_ * (piY_/piR_) * exp1_ - (piG_/piR_) * k2_2 * exp2_); //G
        case 3 : return r_2 * (piT_ *             - exp1_);                              //T, U
      }
    } 
    //C
    case 1 : {
      switch(j) {
        case 0 : return r_2 * (piA_ *             - exp1_);                              //A
        case 1 : return r_2 * (piC_ * (piR_/piY_) * exp1_ + (piT_/piY_) * k2_2 * exp2_); //C
        case 2 : return r_2 * (piG_ *             - exp1_);                              //G
        case 3 : return r_2 * (piT_ * (piR_/piY_) * exp1_ - (piT_/piY_) * k2_2 * exp2_); //T, U
      }
    }
    //G
    case 2 : {
      switch(j) {
        case 0 : return r_2 * (piA_ * (piY_/piR_) * exp1_ - (piA_/piR_) * k2_2 * exp2_); //A
        case 1 : return r_2 * (piC_ *             - exp1_);                              //C
        case 2 : return r_2 * (piG_ * (piY_/piR_) * exp1_ + (piA_/piR_) * k2_2 * exp2_); //G
        case 3 : return r_2 * (piT_ *             - exp1_);                              //T, U
      }
    }
    //T, U
    case 3 : {
      switch(j) {
        case 0 : return r_2 * (piA_ *             - exp1_);                              //A
        case 1 : return r_2 * (piC_ * (piR_/piY_) * exp1_ - (piC_/piY_) * k2_2 * exp2_); //C
        case 2 : return r_2 * (piG_ *             - exp1_);                              //G
        case 3 : return r_2 * (piT_ * (piR_/piY_) * exp1_ + (piC_/piY_) * k2_2 * exp2_); //T, U
      }
    }
  }
  return 0;
}

/******************************************************************************/

const Matrix<double> & F84::getPij_t(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-k1_*l_);
  exp2_ = exp(-k2_*l_);

  //A
  p_(0, 0) = piA_ * (1. + (piY_/piR_) * exp1_) + (piG_/piR_) * exp2_; //A
  p_(0, 1) = piC_ * (1. -               exp1_);                       //C
  p_(0, 2) = piG_ * (1. + (piY_/piR_) * exp1_) - (piG_/piR_) * exp2_; //G
  p_(0, 3) = piT_ * (1. -               exp1_);                       //T, U

  //C
  p_(1, 0) = piA_ * (1. -               exp1_);                       //A
  p_(1, 1) = piC_ * (1. + (piR_/piY_) * exp1_) + (piT_/piY_) * exp2_; //C
  p_(1, 2) = piG_ * (1. -               exp1_);                       //G
  p_(1, 3) = piT_ * (1. + (piR_/piY_) * exp1_) - (piT_/piY_) * exp2_; //T, U

  //G
  p_(2, 0) = piA_ * (1. + (piY_/piR_) * exp1_) - (piA_/piR_) * exp2_; //A
  p_(2, 1) = piC_ * (1. -               exp1_);                       //C
  p_(2, 2) = piG_ * (1. + (piY_/piR_) * exp1_) + (piA_/piR_) * exp2_; //G
  p_(2, 3) = piT_ * (1. -               exp1_);                       //T, U

  //T, U
  p_(3, 0) = piA_ * (1. -               exp1_);                       //A
  p_(3, 1) = piC_ * (1. + (piR_/piY_) * exp1_) - (piC_/piY_) * exp2_; //C
  p_(3, 2) = piG_ * (1. -               exp1_);                       //G
  p_(3, 3) = piT_ * (1. + (piR_/piY_) * exp1_) + (piC_/piY_) * exp2_; //T, U

  return p_;
}

const Matrix<double> & F84::getdPij_dt(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-k1_*l_);
  exp2_ = exp(-k2_*l_);

  //A
  p_(0, 0) = rate_ * r_ * (piA_ * -(piY_/piR_) * exp1_ - (piG_/piR_) * k2_ * exp2_); //A
  p_(0, 1) = rate_ * r_ * (piC_ *                exp1_);                             //C
  p_(0, 2) = rate_ * r_ * (piG_ * -(piY_/piR_) * exp1_ + (piG_/piR_) * k2_ * exp2_); //G
  p_(0, 3) = rate_ * r_ * (piT_ *                exp1_);                             //T, U

  //C
  p_(1, 0) = rate_ * r_ * (piA_ *                exp1_);                             //A
  p_(1, 1) = rate_ * r_ * (piC_ * -(piR_/piY_) * exp1_ - (piT_/piY_) * k2_ * exp2_); //C
  p_(1, 2) = rate_ * r_ * (piG_ *                exp1_);                             //G
  p_(1, 3) = rate_ * r_ * (piT_ * -(piR_/piY_) * exp1_ + (piT_/piY_) * k2_ * exp2_); //T, U

  //G
  p_(2, 0) = rate_ * r_ * (piA_ * -(piY_/piR_) * exp1_ + (piA_/piR_) * k2_ * exp2_); //A
  p_(2, 1) = rate_ * r_ * (piC_ *                exp1_);                             //C
  p_(2, 2) = rate_ * r_ * (piG_ * -(piY_/piR_) * exp1_ - (piA_/piR_) * k2_ * exp2_); //G
  p_(2, 3) = rate_ * r_ * (piT_ *                exp1_);                             //T, U

  //T, U
  p_(3, 0) = rate_ * r_ * (piA_ *                exp1_);                             //A
  p_(3, 1) = rate_ * r_ * (piC_ * -(piR_/piY_) * exp1_ + (piC_/piY_) * k2_ * exp2_); //C
  p_(3, 2) = rate_ * r_ * (piG_ *                exp1_);                             //G
  p_(3, 3) = rate_ * r_ * (piT_ * -(piR_/piY_) * exp1_ - (piC_/piY_) * k2_ * exp2_); //T, U

  return p_;
}

const Matrix<double> & F84::getd2Pij_dt2(double d) const
{
  double r_2 = rate_ * rate_ * r_ * r_;
  l_ = rate_ * r_ * d;
  double k2_2 = k2_ * k2_;
  exp1_ = exp(-k1_*l_);
  exp2_ = exp(-k2_*l_);

  //A
  p_(0, 0) = r_2 * (piA_ * (piY_/piR_) * exp1_ + (piG_/piR_) * k2_2 * exp2_); //A
  p_(0, 1) = r_2 * (piC_ *             - exp1_);                              //C
  p_(0, 2) = r_2 * (piG_ * (piY_/piR_) * exp1_ - (piG_/piR_) * k2_2 * exp2_); //G
  p_(0, 3) = r_2 * (piT_ *             - exp1_);                              //T, U

  //C
  p_(1, 0) = r_2 * (piA_ *             - exp1_);                              //A
  p_(1, 1) = r_2 * (piC_ * (piR_/piY_) * exp1_ + (piT_/piY_) * k2_2 * exp2_); //C
  p_(1, 2) = r_2 * (piG_ *             - exp1_);                              //G
  p_(1, 3) = r_2 * (piT_ * (piR_/piY_) * exp1_ - (piT_/piY_) * k2_2 * exp2_); //T, U

  //G
  p_(2, 0) = r_2 * (piA_ * (piY_/piR_) * exp1_ - (piA_/piR_) * k2_2 * exp2_); //A
  p_(2, 1) = r_2 * (piC_ *             - exp1_);                              //C
  p_(2, 2) = r_2 * (piG_ * (piY_/piR_) * exp1_ + (piA_/piR_) * k2_2 * exp2_); //G
  p_(2, 3) = r_2 * (piT_ *             - exp1_);                              //T, U
 
  //T, U
  p_(3, 0) = r_2 * (piA_ *             - exp1_);                              //A
  p_(3, 1) = r_2 * (piC_ * (piR_/piY_) * exp1_ - (piC_/piY_) * k2_2 * exp2_); //C
  p_(3, 2) = r_2 * (piG_ *             - exp1_);                              //G
  p_(3, 3) = r_2 * (piT_ * (piR_/piY_) * exp1_ + (piC_/piY_) * k2_2 * exp2_); //T, U

  return p_;
}

/******************************************************************************/

void F84::setFreq(map<int, double>& freqs)
{
  piA_ = freqs[0];
  piC_ = freqs[1];
  piG_ = freqs[2];
  piT_ = freqs[3];
  vector<string> thetas(3);
  thetas[0] = getNamespace()+ "theta";
  thetas[1] = getNamespace() + "theta1";
  thetas[2] = getNamespace() + "theta2";
  ParameterList pl = getParameters().subList(thetas);
  pl[0].setValue(piC_ + piG_);
  pl[1].setValue(piA_ / (piA_ + piT_));
  pl[2].setValue(piG_ / (piC_ + piG_));
  setParametersValues(pl);
}

/******************************************************************************/

