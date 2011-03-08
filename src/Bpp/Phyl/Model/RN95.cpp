//
// File: RN95.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 24 février 2011, à 20h 42
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "RN95.h"

#include <Bpp/Numeric/Matrix/MatrixTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

// From the STL:
#include <cmath>

using namespace bpp;
using namespace std;

/******************************************************************************/
 
RN95::RN95(
	const NucleicAlphabet* alphabet,
        double alpha1,
        double alpha2,
        double alpha3,
        double alpha4,
        double beta1,
        double beta2,
        double beta3): 
  AbstractSubstitutionModel(alphabet, "RN95."),
  alpha1_(alpha1),
  alpha2_(alpha2),
  alpha3_(alpha3),
  alpha4_(alpha4),
  beta1_(beta1),
  beta2_(beta2),
  beta3_(beta3)
{
  addParameter_(Parameter("RN95.alpha1" , alpha1 , &Parameter::R_PLUS_STAR));
  addParameter_(Parameter("RN95.alpha2" , alpha2 , &Parameter::R_PLUS_STAR));
  addParameter_(Parameter("RN95.alpha3" , alpha3 , &Parameter::R_PLUS_STAR));
  addParameter_(Parameter("RN95.alpha4" , alpha4 , &Parameter::R_PLUS_STAR));
  addParameter_(Parameter("RN95.beta1" , beta1 , &Parameter::R_PLUS_STAR));
  addParameter_(Parameter("RN95.beta2" , beta2 , &Parameter::R_PLUS_STAR));
  addParameter_(Parameter("RN95.beta3" , beta3 , &Parameter::R_PLUS_STAR));

  updateMatrices();
}

/******************************************************************************/
	
void RN95::updateMatrices()
{
  alpha1_  = getParameterValue("alpha1");
  alpha2_  = getParameterValue("alpha2");
  alpha3_  = getParameterValue("alpha3");
  alpha4_  = getParameterValue("alpha4");
  beta1_  = getParameterValue("beta1");
  beta2_  = getParameterValue("beta2");
  beta3_  = getParameterValue("beta3");

  
  // Generator matrix:

  generator_(0,1) = beta3_;
  generator_(0,2) = alpha4_;
  generator_(0,3) = beta2_;

  generator_(0,0)= -(beta3_+ alpha4_+ beta2_);

  generator_(1,0) = beta1_;
  generator_(1,2) = 1;
  generator_(1,3) = alpha2_;

  generator_(1,1)= -(beta1_+ alpha2_+ 1);

  generator_(2,0) = alpha1_;
  generator_(2,1) = beta3_;
  generator_(2,3) = beta2_;

  generator_(2,2)= -(beta3_+ alpha1_+ beta2_);

  generator_(3,0) = beta1_;
  generator_(3,1) = alpha3_;
  generator_(3,2) = 1;

  generator_(1,1)= -(beta1_+ alpha3_+ 1);

  
  EigenValue<double> ev(generator_);
  rightEigenVectors_ = ev.getV();
  MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);
  Vdouble vi = ev.getImagEigenValues();
  eigenValues_ = ev.getRealEigenValues();
  unsigned int  i;
  double x;

  // looking for the 0 eigenvector for which the eigen vector
  // elements are of the same sign.
  // There is a threshold value in case of approximations

  double seuil=0.00000001;
  bool flag=true;
  unsigned int nulleigen;
  unsigned int salph=4;
  
  while(flag){
    seuil*=10;
    nulleigen=0;
    int signe=0;

    while (flag && (nulleigen < salph)) {
      if (abs(eigenValues_[nulleigen]) < 0.0000001 && abs(vi[nulleigen]) < 0.0000001){
        signe=0;
        i=0;
        while (signe==0 && i< salph){
          x=leftEigenVectors_(nulleigen, i);
          signe=x>seuil?1:x< -seuil?-1:0;
          i++;
        }
        if (signe==0)
          nulleigen++;
        else {
          while (i<salph){
            x=leftEigenVectors_(nulleigen, i);
            if ((signe==-1 && x>seuil) || (signe==1 && x<-seuil))
              break;
                
            i++;
          }
          if (i<salph)
            nulleigen++;
          else
            flag=false;
        }
      }
      else
        nulleigen++;
    }
  }
    
  for (i = 0; i < salph; i++)
    {
      freq_[i] = leftEigenVectors_(nulleigen, i);
    }

  x = 0;
  for (i = 0; i < salph; i++)
    {
      x += freq_[i];
    }

  for (i = 0; i < salph; i++)
    {
      freq_[i] /= x;
    }

  // normalization

  x = 0;
  for (i = 0; i < salph; i++)
    {
      x += freq_[i] * generator_(i, i);
    }

  MatrixTools::scale(generator_, -1. / x);

  for (i = 0; i < salph; i++)
    {
      eigenValues_[i] /= -x;
    }
  
}

/******************************************************************************/

void RN95::setFreq(map<int, double>& freqs)
{
}

/******************************************************************************/

