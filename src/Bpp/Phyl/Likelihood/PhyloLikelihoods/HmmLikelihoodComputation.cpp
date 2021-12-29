//
// File: HmmLikelihoodComputation.cpp
// Authors:
//   Julien Dutheil
// Created: 2007-10-26 11:57:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
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
using namespace numeric;

/***********************************/
/* Forward Likelihood */
/***********************************/

void ForwardHmmLikelihood_DF::compute()
{
  const auto& hmmEq = accessValueConstCast<Eigen::VectorXd>(*this->dependency(0));

  const auto& hmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(1));
  const auto& hmmEmis = accessValueConstCast<MatrixLik>(*this->dependency(2));

  auto& condLik = dynamic_pointer_cast<CondLikelihood>(condLik_)->accessValueMutable();

  auto nbSites = hmmEmis.cols();
  auto nbStates = hmmEmis.rows();

  VDataLik tscales((size_t)nbSites);

  VectorLik tmp((int)nbStates);

  // Initialisation:
  parCondLik_[0] = hmmTrans * hmmEq;

  cwise(tmp) = cwise(parCondLik_[0]) * cwise(hmmEmis.col(0));
  tscales[0] = tmp.sum();

  for (auto s=0;s<nbStates;s++)
    condLik(s,0) = convert(tmp(s) / tscales[0]);
  // tmp = condLik * scales

  // Iteration
  for (auto i = 1; i < nbSites; i++)
  {
    parCondLik_[(size_t)i] =  hmmTrans * condLik.col(i - 1);

    cwise(tmp) = cwise(parCondLik_[(size_t)i]) * cwise(hmmEmis.col(i));    
    tscales[(size_t)i] = tmp.sum();

    // tmp = condLik * scales
    for (auto s=0;s<nbStates;s++)
      condLik(s,i) = convert(tmp(s) / tscales[(size_t)i]);
  }

  copyBppToEigen(tscales, this->accessValueMutable ());
}

NodeRef ForwardHmmLikelihood_DF::derive (Context& c, const Node_DF& node)
{
  if (&node == this)
  {
    return ConstantOne<RowLik>::create (c, targetDimension_);
  }

  /*
   * 1st order derivatives of Forward Likelihood Arrays
   *
   * Dependencies are:
   *  Value<VectorXd> : Starting vector of states probabililies
   *  Value<MatrixXd> : TransitionMatrix
   *  Value<MatrixLik> : Matrix of Emission likelihoods states X sites
   *
   *  ForwardHmmLikelihood_DF : Forward Computations
   *
   *  Value<VectorXd> : Derivatives of starting vector of states probabililies
   *  Value<MatrixXd> : Derivatives of TransitionMatrix
   *  Value<MatrixLik> : Derivatives Matrix of Emission likelihoods states X sites
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
  const auto& hmmEmis = accessValueConstCast<MatrixLik>(*this->dependency(2));


  auto forwardNode = dynamic_pointer_cast<ForwardHmmLikelihood_DF>(this->dependency(3));

  const auto& condLik = forwardNode->getForwardCondLikelihood()->getTargetValue();

  const auto& parCondLik = forwardNode->getParCondLik();

  const auto& scales = forwardNode->getTargetValue();

  const auto& dHmmEq = accessValueConstCast<Eigen::VectorXd>(*this->dependency(4));

  const auto& dHmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(5));

  const auto& dHmmEmis = accessValueConstCast<MatrixLik>(*this->dependency(6));

  auto nbSites = hmmEmis.cols();
  const int nbStates = (int)hmmEmis.rows();
  
  VDataLik tdScales((size_t)nbSites);

  auto& dCondLik = dynamic_pointer_cast<CondLikelihood>(dCondLik_)->accessValueMutable();

  VectorLik dtmp(nbStates);

  // Initialisation
  dParCondLik_[0] = dHmmTrans * hmmEq + hmmTrans * dHmmEq;

  cwise(dtmp) = (cwise(dParCondLik_[0]) * cwise(hmmEmis.col(0))
                +  cwise(dHmmEmis.col(0)) *  cwise(parCondLik[0]));
  tdScales[0] = dtmp.sum();

  // dtmp = dCondLik * scales + CondLik * dScales

  for (auto s=0;s<nbStates;s++)
    dCondLik(s,0)=convert((dtmp(s) - condLik(s,0) * tdScales[0])/scales(0));
  
  // Iteration
  for (auto i = 1; i < nbSites; i++)
  {
    dParCondLik_[(size_t)i] = dHmmTrans * condLik.col(i - 1) + hmmTrans * dCondLik.col(i - 1);

    cwise(dtmp) = cwise(dParCondLik_[(size_t)i]) * cwise(hmmEmis.col(i)) + cwise(parCondLik[(size_t)i]) * cwise(dHmmEmis.col(i));
    tdScales[(size_t)i] = dtmp.sum();

    for (auto s=0;s<nbStates;s++)
      dCondLik(s,i)=convert((dtmp(s) - condLik(s,i) * tdScales[(size_t)i])/scales(i));
  }

  copyBppToEigen(tdScales, this->accessValueMutable ());
}

NodeRef ForwardHmmDLikelihood_DF::derive (Context& c, const Node_DF& node)
{
  if (&node == this)
  {
    return ConstantOne<RowLik>::create (c, targetDimension_);
  }

  /*
   *
   * 2nd order derivatives of Forward Likelihood Arrays
   *
   * Dependencies are:
   *  Value<VectorXd> : Starting vector of states probabililies
   *  Value<MatrixXd> : TransitionMatrix
   *  Value<MatrixLik> : Matrix of Emission likelihoods states X sites
   *
   *  ForwardHmmLikelihood_DF : Forward Computations
   *
   *  Value<VectorXd> : 1st Derivatives of starting vector of states probabililies
   *  Value<MatrixXd> : 1st Derivatives of TransitionMatrix
   *  Value<MatrixLik> : 1st Derivatives Matrix of Emission likelihoods states X sites
   *
   *  ForwardHmmDLikelihood_DF : 1st order derivatives Forward Computations
   *
   *  Value<VectorXd> : 2nd Derivatives of starting vector of states probabililies
   *  Value<MatrixXd> : 2nd Derivatives of TransitionMatrix
   *  Value<MatrixLik> : 2nd Derivatives Matrix of Emission likelihoods states X sites
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
  const auto& hmmEmis = accessValueConstCast<MatrixLik>(*this->dependency(2));


  auto forwardNode = dynamic_pointer_cast<ForwardHmmLikelihood_DF>(this->dependency(3));

  const auto& condLik = forwardNode->getForwardCondLikelihood()->getTargetValue();

  const auto& parCondLik = forwardNode->getParCondLik();

  const auto& scales = forwardNode->getTargetValue();


  const auto& dHmmEq = accessValueConstCast<Eigen::VectorXd>(*this->dependency(4));

  const auto& dHmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(5));

  const auto& dHmmEmis = accessValueConstCast<MatrixLik>(*this->dependency(6));


  auto forwardDNode = dynamic_pointer_cast<ForwardHmmDLikelihood_DF>(this->dependency(7));

  const auto& dCondLik = forwardDNode->getForwardDCondLikelihood()->getTargetValue();

  const auto& parDCondLik = forwardDNode->getParDCondLik();

  const auto& dScales = forwardDNode->getTargetValue();


  const auto& d2HmmEq = accessValueConstCast<Eigen::VectorXd>(*this->dependency(8));

  const auto& d2HmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(9));

  const auto& d2HmmEmis = accessValueConstCast<MatrixLik>(*this->dependency(10));


  //////////////////////////////////////

  auto nbSites = hmmEmis.cols();
  const int nbStates = (int)hmmEmis.rows();

  VDataLik td2Scales((size_t)nbSites);

  Eigen::VectorXd d2CondLik((int)hmmEmis.rows());
 
  Eigen::VectorXd d2ParCondLik;
  VectorLik d2tmp(nbStates);

  // Initialisation:
  d2ParCondLik = d2HmmTrans * hmmEq + 2 * dHmmTrans * dHmmEq + hmmTrans * d2HmmEq;

  cwise(d2tmp) = cwise(d2ParCondLik) * cwise(hmmEmis.col(0))
               + 2 * cwise(parDCondLik[0]) * cwise(dHmmEmis.col(0))
               + cwise(parCondLik[0]) * cwise(d2HmmEmis.col(0));

  td2Scales[0] = d2tmp.sum();

  // d2tmp = d2CondLik * scales + 2 * dCondLik * dScales + CondLik * d2Scales
  
  for (auto s=0;s<nbStates;s++)
    d2CondLik(s,0) = convert((d2tmp(s) - 2 * dCondLik(s,0) * dScales(0) - condLik(s,0) * td2Scales[0]) / scales(0));
 
  // Iteration
  for (auto i = 1; i < nbSites; i++)
  {
    d2ParCondLik = d2HmmTrans * condLik.col(i)
                   + 2 * dHmmTrans * dCondLik.col(i) + hmmTrans * d2CondLik;

    cwise(d2tmp) = cwise(d2ParCondLik) * cwise(hmmEmis.col(i))
                 + 2 * cwise(parDCondLik[(size_t)i]) * cwise(dHmmEmis.col(i))
                 + cwise(parCondLik[(size_t)i]) * cwise(d2HmmEmis.col(i));

    td2Scales[(size_t)i] = d2tmp.sum();

    for (auto s=0;s<nbStates;s++)
      d2CondLik(s,i) = convert((d2tmp(s) - 2 * dCondLik(s,i) * dScales(i) - condLik(s,i) * td2Scales[(size_t)i]) / scales(i));

  }

  copyBppToEigen(td2Scales, this->accessValueMutable ());
}

/*****************************************
** Backward
*****************************************/
void BackwardHmmLikelihood_DF::compute()
{
  const auto& scales = accessValueConstCast<RowLik>(*this->dependency(0));

  const auto& hmmTrans = accessValueConstCast<Eigen::MatrixXd>(*this->dependency(1));
  const auto& hmmEmis = accessValueConstCast<MatrixLik>(*this->dependency(2));


  auto nbSites = hmmEmis.cols();
  auto nbStates = hmmEmis.rows();

  std::vector<Eigen::VectorXd> tScales((size_t)nbSites);

  VectorLik tmp((int)nbStates);

  // Initialisation:
  tScales[(size_t)(nbSites - 1)] = Eigen::VectorXd::Constant(nbStates, 1);

  // Iteration
  for (auto i = nbSites - 1; i > 0; i--)
  {
    tScales[(size_t)(i - 1)].resize(nbStates);
    
    cwise(tmp) = cwise(hmmEmis.col(i)) * cwise(tScales[(size_t)i]);

    auto tmp2 = hmmTrans * tmp;

    for (auto s=0;s<nbStates;s++)
      tScales[(size_t)(i - 1)](s) = convert(tmp2(s)/scales(i));
  }

  copyBppToEigen(tScales, this->accessValueMutable ());
}
