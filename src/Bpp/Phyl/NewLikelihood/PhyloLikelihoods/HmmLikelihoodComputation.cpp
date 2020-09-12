//
// File: HmmLikelihood_DF.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 26 11:57 2007
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

#include "HmmLikelihood_DF.h"

// from the STL:
#include <iostream>
#include <algorithm>
using namespace bpp;
using namespace std;

/***********************************/
/* Forward Likelihood */
/***********************************/
   
void ForwardHmmLikelihood_DF::compute()
{
  const auto& hmmEq = accessValueConstCast<Eigen::VectorXd>(*this->dependency(0));

  const auto& hmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(1));
  const auto& hmmEmis = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(2));

  auto & scales = this->accessValueMutable ();

  auto& condLik = dynamic_pointer_cast<CondLikelihood>(condLik_)->accessValueMutable();

  size_t nbSites = (size_t)hmmEmis.cols();

  Eigen::VectorXd tmp;
  
  //Initialisation:

  Eigen::Ref<Eigen::VectorXd>  col0 = parCondLik_.col(0);
  
  col0 = hmmTrans * hmmEq;

  col0.array() *= hmmEmis.col(0).array();

  scales(0) = col0.sum();
  
  condLik.col(0) = col0 / scales(0);

  //Iteration
  for (size_t i = 1; i < nbSites; i++)
  {
    Eigen::Ref<Eigen::VectorXd> coli = parCondLik_.col(i);
    
    coli =  hmmTrans * condLik.col(i-1);

    coli.array() *= hmmEmis.col(i).array();

    scales(i) = coli.sum();

    condLik.col(i) = coli / scales(i);
  }
  
}

NodeRef ForwardHmmLikelihood_DF::derive (Context & c, const Node_DF & node)
{
  if (&node == this) {
    return ConstantOne<Eigen::RowVectorXd>::create (c, targetDimension_);
  }
    
  /*
   * 1st order derivatives of Forward Likelihood Arrays
   *
   * Dependencies are:
   *  Value<VectorXd> : Starting vector of states probabililies
   *  Value<MatrixXd> : TransitionMatrix
   *  Value<MatrixXd> : Matrix of Emission likelihoods states X sites
   *
   *  ForwardHmmLikelihood_DF : Forward Computations
   *
   *  Value<VectorXd> : Derivatives of starting vector of states probabililies
   *  Value<MatrixXd> : Derivatives of TransitionMatrix
   *  Value<MatrixXd> : Derivatives Matrix of Emission likelihoods states X sites
   */
      
  return ForwardHmmDLikelihood_DF::create (c, {
      this->dependency(0),
        this->dependency(1),
        this->dependency(2),
        this->shared_from_this(),
        this->dependency(0)->derive (c, node),
        this->dependency(1)->derive (c, node),
        this->dependency(2)->derive (c, node)}, targetDimension_);
}

/***********************************/
/* 1st order derivative of Forward Likelihood */
/***********************************/
   

void ForwardHmmDLikelihood_DF::compute()
{
  const auto& hmmEq = accessValueConstCast<Eigen::VectorXd>(*this->dependency(0));

  const auto& hmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(1));
  const auto& hmmEmis = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(2));

  
  auto forwardNode = dynamic_pointer_cast<ForwardHmmLikelihood_DF>(this->dependency(3));
  
  const auto& condLik = forwardNode -> getForwardCondLikelihood()->getTargetValue();

  const auto& parCondLik = forwardNode -> getParCondLik();
  
  const auto& scales = forwardNode -> getTargetValue();
  

  
  const auto& dHmmEq = accessValueConstCast<Eigen::VectorXd>(*this->dependency(4));

  const auto& dHmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(5));
  
  const auto& dHmmEmis = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(6));

  auto & dScales = this->accessValueMutable ();

  
  auto& dCondLik = dynamic_pointer_cast<CondLikelihood>(dCondLik_)->accessValueMutable();

  size_t nbSites = (size_t)hmmEmis.cols();

  Eigen::VectorXd tmp(hmmEmis.rows());
  
  //Initialisation:
  Eigen::Ref<Eigen::VectorXd>  col0 = dParCondLik_.col(0);

  col0 = dHmmTrans * hmmEq + hmmTrans * dHmmEq;

  tmp.array() = col0.array() * hmmEmis.col(0).array()
    +  dHmmEmis.col(0).array() *  parCondLik.col(0).array();

  dScales(0) = tmp.sum();

  dCondLik.col(0) = (tmp - condLik.col(0) * dScales(0))/scales(0);

  //Iteration
  for (size_t i = 1; i < nbSites; i++)
  {
    Eigen::Ref<Eigen::VectorXd>  coli = dParCondLik_.col(i);
    
    coli = dHmmTrans * condLik.col(i-1) + hmmTrans * dCondLik.col(i-1);

    tmp.array() = coli.array() * hmmEmis.col(i).array() + parCondLik.col(i).array() * dHmmEmis.col(i).array();
  
    dScales(i) = tmp.sum();
  
    dCondLik.col(i) = (tmp - condLik.col(i) * dScales(i))/scales(i);
  }
  
}

NodeRef ForwardHmmDLikelihood_DF::derive (Context & c, const Node_DF & node)
{
  if (&node == this) {
    return ConstantOne<Eigen::RowVectorXd>::create (c, targetDimension_);
  }

  /*
   *
   * 2nd order derivatives of Forward Likelihood Arrays
   *
   * Dependencies are:
   *  Value<VectorXd> : Starting vector of states probabililies
   *  Value<MatrixXd> : TransitionMatrix
   *  Value<MatrixXd> : Matrix of Emission likelihoods states X sites
   *
   *  ForwardHmmLikelihood_DF : Forward Computations
   *
   *  Value<VectorXd> : 1st Derivatives of starting vector of states probabililies
   *  Value<MatrixXd> : 1st Derivatives of TransitionMatrix
   *  Value<MatrixXd> : 1st Derivatives Matrix of Emission likelihoods states X sites
   *
   *  ForwardHmmDLikelihood_DF : 1st order derivatives Forward Computations
   *
   *  Value<VectorXd> : 2nd Derivatives of starting vector of states probabililies
   *  Value<MatrixXd> : 2nd Derivatives of TransitionMatrix
   *  Value<MatrixXd> : 2nd Derivatives Matrix of Emission likelihoods states X sites
   */
      
  return ForwardHmmD2Likelihood_DF::create (c, {
      this->dependency(0),
        this->dependency(1),
        this->dependency(2),

        this->dependency(3),

        this->dependency(4),
        this->dependency(5),
        this->dependency(6),
        
        this->shared_from_this(),
        
        this->dependency(4)->derive (c, node),
        this->dependency(5)->derive (c, node),
        this->dependency(6)->derive (c, node)},

    targetDimension_);

}



/***********************************/
/* 2nd order derivative of Forward Likelihood */
/***********************************/
   

void ForwardHmmD2Likelihood_DF::compute()
{
  const auto& hmmEq = accessValueConstCast<Eigen::VectorXd>(*this->dependency(0));

  const auto& hmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(1));
  const auto& hmmEmis = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(2));

  
  auto forwardNode = dynamic_pointer_cast<ForwardHmmLikelihood_DF>(this->dependency(3));
  
  const auto& condLik = forwardNode -> getForwardCondLikelihood()->getTargetValue();

  const auto& parCondLik = forwardNode -> getParCondLik();
  
  const auto& scales = forwardNode -> getTargetValue();
  

  
  const auto& dHmmEq = accessValueConstCast<Eigen::VectorXd>(*this->dependency(4));

  const auto& dHmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(5));
  
  const auto& dHmmEmis = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(6));


  auto forwardDNode = dynamic_pointer_cast<ForwardHmmDLikelihood_DF>(this->dependency(7));
  
  const auto& dCondLik = forwardDNode -> getForwardDCondLikelihood()->getTargetValue();

  const auto& parDCondLik = forwardDNode -> getParDCondLik();
  
  const auto& dScales = forwardDNode -> getTargetValue();
  


  const auto& d2HmmEq = accessValueConstCast<Eigen::VectorXd>(*this->dependency(8));

  const auto& d2HmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(9));
  
  const auto& d2HmmEmis = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(10));


  //////////////////////////////////////
  
  auto & d2Scales = this->accessValueMutable ();
  
  Eigen::VectorXd d2CondLik(hmmEmis.rows());

  size_t nbSites = (size_t)hmmEmis.cols();

  Eigen::VectorXd d2ParCondLik, tmp(hmmEmis.rows());
  
  //Initialisation:
  d2ParCondLik = d2HmmTrans * hmmEq + 2 * dHmmTrans * dHmmEq + hmmTrans * d2HmmEq;

  tmp.array() = d2ParCondLik.array() * hmmEmis.col(0).array()
    + 2 * parDCondLik.col(0).array() * dHmmEmis.col(0).array()
    + parCondLik.col(0).array() * d2HmmEmis.col(0).array();
  
  d2Scales(0) = tmp.sum();
  
  d2CondLik = (tmp - 2 * dCondLik.col(0) * dScales(0) - condLik.col(0) * d2Scales(0))/scales(0);

  //Iteration
  for (size_t i = 1; i < nbSites; i++)
  {
    d2ParCondLik = d2HmmTrans * condLik.col(i)
      + 2 * dHmmTrans * dCondLik.col(i) + hmmTrans * d2CondLik;

    tmp.array() = d2ParCondLik.array() * hmmEmis.col(i).array()
      + 2 * parDCondLik.col(i).array() * dHmmEmis.col(i).array()
      + parCondLik.col(i).array() * d2HmmEmis.col(i).array();
  
    d2Scales(i) = tmp.sum();
  
    d2CondLik = (tmp - 2 * dCondLik.col(i) * dScales(i) - condLik.col(i) * d2Scales(i))/scales(i);
  }
}

/*****************************************
 ** Backward 
 *****************************************/


void BackwardHmmLikelihood_DF::compute()
{
  const auto& scales = accessValueConstCast<Eigen::RowVectorXd>(*this->dependency(0));

  const auto& hmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(1));
  const auto& hmmEmis = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(2));

  auto & backward = this->accessValueMutable ();

  size_t nbSites = (size_t)hmmEmis.cols();
  size_t nbStates = (size_t)hmmEmis.rows();

  Eigen::VectorXd tmp(nbStates);
  
  //Initialisation:
  backward.col(nbSites-1) = Eigen::VectorXd::Constant(nbStates, 1);

  //Iteration
  for (size_t i = nbSites-1; i > 0; i--)
  {
    tmp.array() =  hmmEmis.col(i).array() * backward.col(i).array();

    backward.col(i-1) = hmmTrans * tmp;

    backward.col(i-1) /= scales(i);
  }
  
}


