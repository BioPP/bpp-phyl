//
// File: TN93.cpp
// Created by: Julien Dutheil
// Created on: Thu Jan 22 10:26:51 2004
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

#include "TN93.h"
#include "../FrequenciesSet/NucleotideFrequenciesSet.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Matrix/EigenValue.h>

// From bpp-seq:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

// From the STL:
#include <cmath>

using namespace std;

using namespace bpp;

/******************************************************************************/

TN93::TN93(
  const NucleicAlphabet* alpha,
  double kappa1,
  double kappa2,
  double piA,
  double piC,
  double piG,
  double piT):
  AbstractParameterAliasable("TN93."),
  AbstractReversibleNucleotideSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "TN93."),
  kappa1_(kappa1), kappa2_(kappa2), 
  piA_(piA), piC_(piC), piG_(piG), piT_(piT), piY_(), piR_(),
  r_(), k1_(), k2_(),
  theta_(piG + piC), theta1_(piA / (1. - theta_)), theta2_(piG / theta_),
  exp1_(), exp21_(), exp22_(), l_(), p_(size_, size_)
{
  addParameter_(new Parameter("TN93.kappa1", kappa1, &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("TN93.kappa2", kappa2, &Parameter::R_PLUS_STAR));
  addParameter_(new Parameter("TN93.theta" , theta_ , &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter("TN93.theta1", theta1_, &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  addParameter_(new Parameter("TN93.theta2", theta2_, &FrequenciesSet::FREQUENCE_CONSTRAINT_SMALL));
  p_.resize(size_, size_);
  updateMatrices();
}

/******************************************************************************/

void TN93::updateMatrices()
{
  kappa1_ = getParameterValue("kappa1");
  kappa2_ = getParameterValue("kappa2");
  theta_  = getParameterValue("theta" );
  theta1_ = getParameterValue("theta1");
  theta2_ = getParameterValue("theta2");
  piA_ = theta1_ * (1. - theta_);
  piC_ = (1. - theta2_) * theta_;
  piG_ = theta2_ * theta_;
  piT_ = (1. - theta1_) * (1. - theta_);
  piR_ = piA_ + piG_;
  piY_ = piT_ + piC_;
  r_ = 1. / (2. * (piA_ * piC_ + piC_ * piG_ + piA_ * piT_ + piG_ * piT_ + kappa2_ * piC_ * piT_ + kappa1_ * piA_ * piG_));
  k1_ = kappa2_ * piY_ + piR_;
  k2_ = kappa1_ * piR_ + piY_;

  freq_[0] = piA_;
  freq_[1] = piC_;
  freq_[2] = piG_;
  freq_[3] = piT_;
    
  generator_(0, 0) = -(                       piC_ + kappa1_*piG_ +         piT_);
  generator_(1, 1) = -(        piA_ +                        piG_ + kappa2_*piT_); 
  generator_(2, 2) = -(kappa1_*piA_ +         piC_                +         piT_);
  generator_(3, 3) = -(        piA_ + kappa2_*piC_ +         piG_               );

  generator_(1, 0) = piA_;
  generator_(3, 0) = piA_;
  generator_(0, 1) = piC_;
  generator_(2, 1) = piC_;
  generator_(1, 2) = piG_;
  generator_(3, 2) = piG_;
  generator_(0, 3) = piT_;
  generator_(2, 3) = piT_;
  
  generator_(2, 0) = kappa1_ * piA_;
  generator_(3, 1) = kappa2_ * piC_;
  generator_(0, 2) = kappa1_ * piG_;
  generator_(1, 3) = kappa2_ * piT_;
  
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

  // We are not sure that the values are computed in this order :p
  // Eigen values:
  //eigenValues_[0] = 0;
  //eigenValues_[1] = -r * (kappa2 * piY + piR);
  //eigenValues_[2] = -r * (kappa1 * piR + piY); 
  //eigenValues_[3] = -r;
  
  // Eigen vectors and values:
  EigenValue<double> ev(generator_);
  rightEigenVectors_ = ev.getV();
  MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
  eigenValues_ = ev.getRealEigenValues();
}
  
/******************************************************************************/

double TN93::Pij_t(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);
  exp22_ = exp(-k2_ * l_);
  exp21_ = exp(-k1_ * l_);
  
  switch(i)
  {
    //A
    case 0 : {
      switch(j) {
        case 0 : return piA_ * (1. + (piY_/piR_) * exp1_) + (piG_/piR_) * exp22_; //A
        case 1 : return piC_ * (1. -               exp1_);                        //C
        case 2 : return piG_ * (1. + (piY_/piR_) * exp1_) - (piG_/piR_) * exp22_; //G
        case 3 : return piT_ * (1. -               exp1_);                        //T, U
      }
    } 
    //C
    case 1 : {
      switch(j) {
        case 0 : return piA_ * (1. -               exp1_);                         //A
        case 1 : return piC_ * (1. + (piR_/piY_) * exp1_) + (piT_/piY_) * exp21_; //C
        case 2 : return piG_ * (1. -               exp1_);                        //G
        case 3 : return piT_ * (1. + (piR_/piY_) * exp1_) - (piT_/piY_) * exp21_; //T, U
      }
    }
    //G
    case 2 : {
      switch(j) {
        case 0 : return piA_ * (1. + (piY_/piR_) * exp1_) - (piA_/piR_) * exp22_; //A
        case 1 : return piC_ * (1. -               exp1_);                        //C
        case 2 : return piG_ * (1. + (piY_/piR_) * exp1_) + (piA_/piR_) * exp22_; //G
        case 3 : return piT_ * (1. -               exp1_);                        //T, U
      }
    }
    //T, U
    case 3 : {
      switch(j) {
        case 0 : return piA_ * (1. -               exp1_);                        //A
        case 1 : return piC_ * (1. + (piR_/piY_) * exp1_) - (piC_/piY_) * exp21_; //C
        case 2 : return piG_ * (1. -               exp1_);                        //G
        case 3 : return piT_ * (1. + (piR_/piY_) * exp1_) + (piC_/piY_) * exp21_; //T, U
      }
    }
  }
  return 0;
}

/******************************************************************************/

double TN93::dPij_dt(size_t i, size_t j, double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);
  exp22_ = exp(-k2_ * l_);
  exp21_ = exp(-k1_ * l_);
  
  switch(i)
  {
    //A
    case 0 : {
      switch(j) {
        case 0 : return rate_ * r_ * (piA_ * -(piY_/piR_) * exp1_ - (piG_/piR_) * k2_ * exp22_); //A
        case 1 : return rate_ * r_ * (piC_ *                exp1_);                              //C
        case 2 : return rate_ * r_ * (piG_ * -(piY_/piR_) * exp1_ + (piG_/piR_) * k2_ * exp22_); //G
        case 3 : return rate_ * r_ * (piT_ *                exp1_);                              //T, U
      }
    } 
    //C
    case 1 : {
      switch(j) {
        case 0 : return rate_ * r_ * (piA_ *                exp1_);                              //A
        case 1 : return rate_ * r_ * (piC_ * -(piR_/piY_) * exp1_ - (piT_/piY_) * k1_ * exp21_); //C
        case 2 : return rate_ * r_ * (piG_ *                exp1_);                              //G
        case 3 : return rate_ * r_ * (piT_ * -(piR_/piY_) * exp1_ + (piT_/piY_) * k1_ * exp21_); //T, U
      }
    }
    //G
    case 2 : {
      switch(j) {
        case 0 : return rate_ * r_ * (piA_ * -(piY_/piR_) * exp1_ + (piA_/piR_) * k2_ * exp22_); //A
        case 1 : return rate_ * r_ * (piC_ *                exp1_);                              //C
        case 2 : return rate_ * r_ * (piG_ * -(piY_/piR_) * exp1_ - (piA_/piR_) * k2_ * exp22_); //G
        case 3 : return rate_ * r_ * (piT_ *                exp1_);                              //T, U
      }
    }
    //T, U
    case 3 : {
      switch(j) {
        case 0 : return rate_ * r_ * (piA_ *                exp1_);                              //A
        case 1 : return rate_ * r_ * (piC_ * -(piR_/piY_) * exp1_ + (piC_/piY_) * k1_ * exp21_); //C
        case 2 : return rate_ * r_ * (piG_ *                exp1_);                              //G
        case 3 : return rate_ * r_ * (piT_ * -(piR_/piY_) * exp1_ - (piC_/piY_) * k1_ * exp21_); //T, U
      }
    }
  }
  return 0;
}

/******************************************************************************/

double TN93::d2Pij_dt2(size_t i, size_t j, double d) const
{
  double r_2 = rate_ * rate_ * r_ * r_;
  l_ = rate_ * r_ * d;
  double k1_2 = k1_ * k1_;
  double k2_2 = k2_ * k2_;
  exp1_ = exp(-l_);
  exp22_ = exp(-k2_ * l_);
  exp21_ = exp(-k1_ * l_);

  switch(i)
  {
    //A
    case 0 : {
      switch(j) {
        case 0 : return r_2 * (piA_ * (piY_/piR_) * exp1_ + (piG_/piR_) * k2_2 * exp22_); //A
        case 1 : return r_2 * (piC_ *             - exp1_);                               //C
        case 2 : return r_2 * (piG_ * (piY_/piR_) * exp1_ - (piG_/piR_) * k2_2 * exp22_); //G
        case 3 : return r_2 * (piT_ *             - exp1_);                               //T, U
      }
    }
    //C
    case 1 : {
      switch(j) {
        case 0 : return r_2 * (piA_ *             - exp1_);                               //A
        case 1 : return r_2 * (piC_ * (piR_/piY_) * exp1_ + (piT_/piY_) * k1_2 * exp21_); //C
        case 2 : return r_2 * (piG_ *             - exp1_);                               //G
        case 3 : return r_2 * (piT_ * (piR_/piY_) * exp1_ - (piT_/piY_) * k1_2 * exp21_); //T, U
      }
    }
    //G
    case 2 : {
      switch(j) {
        case 0 : return r_2 * (piA_ * (piY_/piR_) * exp1_ - (piA_/piR_) * k2_2 * exp22_); //A
        case 1 : return r_2 * (piC_ *             - exp1_);                               //C
        case 2 : return r_2 * (piG_ * (piY_/piR_) * exp1_ + (piA_/piR_) * k2_2 * exp22_); //G
        case 3 : return r_2 * (piT_ *             - exp1_);                               //T, U
      }
    }
    //T, U
    case 3 : {
      switch(j) {
        case 0 : return r_2 * (piA_ *             - exp1_);                               //A
        case 1 : return r_2 * (piC_ * (piR_/piY_) * exp1_ - (piC_/piY_) * k1_2 * exp21_); //C
        case 2 : return r_2 * (piG_ *             - exp1_);                               //G
        case 3 : return r_2 * (piT_ * (piR_/piY_) * exp1_ + (piC_/piY_) * k1_2 * exp21_); //T, U
      }
    }
  }
  return 0;
}

/******************************************************************************/

const Matrix<double> & TN93::getPij_t(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);
  exp22_ = exp(-k2_ * l_);
  exp21_ = exp(-k1_ * l_);

  //A
  p_(0, 0) = piA_ * (1. + (piY_/piR_) * exp1_) + (piG_/piR_) * exp22_; //A
  p_(0, 1) = piC_ * (1. -               exp1_);                        //C
  p_(0, 2) = piG_ * (1. + (piY_/piR_) * exp1_) - (piG_/piR_) * exp22_; //G
  p_(0, 3) = piT_ * (1. -               exp1_);                        //T, U

  //C
  p_(1, 0) = piA_ * (1. -               exp1_);                        //A
  p_(1, 1) = piC_ * (1. + (piR_/piY_) * exp1_) + (piT_/piY_) * exp21_; //C
  p_(1, 2) = piG_ * (1. -               exp1_);                        //G
  p_(1, 3) = piT_ * (1. + (piR_/piY_) * exp1_) - (piT_/piY_) * exp21_; //T, U

  //G
  p_(2, 0) = piA_ * (1. + (piY_/piR_) * exp1_) - (piA_/piR_) * exp22_; //A
  p_(2, 1) = piC_ * (1. -               exp1_);                        //C
  p_(2, 2) = piG_ * (1. + (piY_/piR_) * exp1_) + (piA_/piR_) * exp22_; //G
  p_(2, 3) = piT_ * (1. -               exp1_);                        //T, U

  //T, U
  p_(3, 0) = piA_ * (1. -               exp1_);                        //A
  p_(3, 1) = piC_ * (1. + (piR_/piY_) * exp1_) - (piC_/piY_) * exp21_; //C
  p_(3, 2) = piG_ * (1. -               exp1_);                        //G
  p_(3, 3) = piT_ * (1. + (piR_/piY_) * exp1_) + (piC_/piY_) * exp21_; //T, U

  return p_;
}

const Matrix<double> & TN93::getdPij_dt(double d) const
{
  l_ = rate_ * r_ * d;
  exp1_ = exp(-l_);
  exp22_ = exp(-k2_ * l_);
  exp21_ = exp(-k1_ * l_);

  //A
  p_(0, 0) = rate_ * r_ * (piA_ * -(piY_/piR_) * exp1_ - (piG_/piR_) * k2_ * exp22_); //A
  p_(0, 1) = rate_ * r_ * (piC_ *                exp1_);                              //C
  p_(0, 2) = rate_ * r_ * (piG_ * -(piY_/piR_) * exp1_ + (piG_/piR_) * k2_ * exp22_); //G
  p_(0, 3) = rate_ * r_ * (piT_ *                exp1_);                              //T, U

  //C
  p_(1, 0) = rate_ * r_ * (piA_ *                exp1_);                              //A
  p_(1, 1) = rate_ * r_ * (piC_ * -(piR_/piY_) * exp1_ - (piT_/piY_) * k1_ * exp21_); //C
  p_(1, 2) = rate_ * r_ * (piG_ *                exp1_);                              //G
  p_(1, 3) = rate_ * r_ * (piT_ * -(piR_/piY_) * exp1_ + (piT_/piY_) * k1_ * exp21_); //T, U

  //G
  p_(2, 0) = rate_ * r_ * (piA_ * -(piY_/piR_) * exp1_ + (piA_/piR_) * k2_ * exp22_); //A
  p_(2, 1) = rate_ * r_ * (piC_ *                exp1_);                              //C
  p_(2, 2) = rate_ * r_ * (piG_ * -(piY_/piR_) * exp1_ - (piA_/piR_) * k2_ * exp22_); //G
  p_(2, 3) = rate_ * r_ * (piT_ *                exp1_);                              //T, U

  //T, U
  p_(3, 0) = rate_ * r_ * (piA_ *                exp1_);                              //A
  p_(3, 1) = rate_ * r_ * (piC_ * -(piR_/piY_) * exp1_ + (piC_/piY_) * k1_ * exp21_); //C
  p_(3, 2) = rate_ * r_ * (piG_ *                exp1_);                              //G
  p_(3, 3) = rate_ * r_ * (piT_ * -(piR_/piY_) * exp1_ - (piC_/piY_) * k1_ * exp21_); //T, U

  return p_;
}

const Matrix<double> & TN93::getd2Pij_dt2(double d) const
{
  double r_2 = rate_ * rate_ * r_ * r_;
  l_ = rate_ * r_ * d;
  double k1_2 = k1_ * k1_;
  double k2_2 = k2_ * k2_;
  exp1_ = exp(-l_);
  exp22_ = exp(-k2_ * l_);
  exp21_ = exp(-k1_ * l_);

  //A
  p_(0, 0) = r_2 * (piA_ * (piY_/piR_) * exp1_ + (piG_/piR_) * k2_2 * exp22_); //A
  p_(0, 1) = r_2 * (piC_ *             - exp1_);                               //C
  p_(0, 2) = r_2 * (piG_ * (piY_/piR_) * exp1_ - (piG_/piR_) * k2_2 * exp22_); //G
  p_(0, 3) = r_2 * (piT_ *             - exp1_);                               //T, U

  //C
  p_(1, 0) = r_2 * (piA_ *             - exp1_);                               //A
  p_(1, 1) = r_2 * (piC_ * (piR_/piY_) * exp1_ + (piT_/piY_) * k1_2 * exp21_); //C
  p_(1, 2) = r_2 * (piG_ *             - exp1_);                               //G
  p_(1, 3) = r_2 * (piT_ * (piR_/piY_) * exp1_ - (piT_/piY_) * k1_2 * exp21_); //T, U

  //G
  p_(2, 0) = r_2 * (piA_ * (piY_/piR_) * exp1_ - (piA_/piR_) * k2_2 * exp22_); //A
  p_(2, 1) = r_2 * (piC_ *             - exp1_);                               //C
  p_(2, 2) = r_2 * (piG_ * (piY_/piR_) * exp1_ + (piA_/piR_) * k2_2 * exp22_); //G
  p_(2, 3) = r_2 * (piT_ *             - exp1_);                               //T, U

  //T, U
  p_(3, 0) = r_2 * (piA_ *             - exp1_);                               //A
  p_(3, 1) = r_2 * (piC_ * (piR_/piY_) * exp1_ - (piC_/piY_) * k1_2 * exp21_); //C
  p_(3, 2) = r_2 * (piG_ *             - exp1_);                               //G
  p_(3, 3) = r_2 * (piT_ * (piR_/piY_) * exp1_ + (piC_/piY_) * k1_2 * exp21_); //T, U

  return p_;
}

/******************************************************************************/

void TN93::setFreq(std::map<int, double>& freqs)
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

