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

#include "HmmLikelihoodComputation.h"

// from the STL:
#include <iostream>
#include <algorithm>
using namespace bpp;
using namespace std;

HmmLikelihood_DF::HmmLikelihood_DF(
  Context& context,
  std::shared_ptr<HmmStateAlphabet> hiddenAlphabet,
  std::shared_ptr<HmmTransitionMatrix> transitionMatrix,
  std::shared_ptr<HmmPhyloEmissionProbabilities> emissionProbabilities,
  const std::string& prefix) :
  AlignedLikelihoodCalculation(context),
  context_(context),
  hiddenAlphabet_(hiddenAlphabet),
  emissionProbabilities_(emissionProbabilities),
  matrix_(),
  hmmEq_(),
  hmmTrans_(),
  hmmEmis_(),
  forwardLik_(),
  backwardLik_(),
  hiddenPostProb_(),
  nbStates_(),
  nbSites_()
{
  if (!hiddenAlphabet)
    throw Exception("HmmLikelihood_DF: null pointer passed for HmmStateAlphabet.");
  if (!transitionMatrix)
    throw Exception("HmmLikelihood_DF: null pointer passed for HmmTransitionMatrix.");
  if (!emissionProbabilities)
    throw Exception("HmmLikelihood_DF: null pointer passed for HmmEmissionProbabilities.");

  nbStates_ = Eigen::Index(hiddenAlphabet_->getNumberOfStates());
  nbSites_ = Eigen::Index(emissionProbabilities_->getNumberOfPositions());

  if (emissionProbabilities_->getNumberOfStates() != size_t(nbStates_))
    throw BadSizeException("HmmLikelihood_DF: HmmStateAlphabet and HmmEmissionProbabilities do not have the same number of states.", emissionProbabilities_->getNumberOfStates(), size_t(nbStates_));

  if (transitionMatrix->getNumberOfStates() != size_t(nbStates_))
    throw BadSizeException("HmmLikelihood_DF: HmmStateAlphabet and HmmTransitionMatrix do not have the same number of states.", transitionMatrix->getNumberOfStates(), size_t(nbStates_));


  ////////////////////////////
  // Adding transition matrix

  // parameters of the transition matrix
  const auto& param = transitionMatrix->getParameters();
  ParameterList paramList;

  for (size_t i = 0; i < param.size(); i++)
  {
    paramList.shareParameter(ConfiguredParameter::create(context_, param[i]));
  }

  shareParameters_(paramList);

  // make TransitionMatrix DF
  matrix_ = ConfiguredParametrizable::createConfigured<HmmTransitionMatrix, ConfiguredTransitionMatrix>(context_, *transitionMatrix, paramList, "");

  // for derivates
  auto deltaNode = NumericMutable<double>::create(context_, 0.001);
  auto config = NumericalDerivativeType::ThreePoints;

  matrix_->config.delta = deltaNode;
  matrix_->config.type = config;

  // equilibrium

  auto eqdim = VectorDimension(Eigen::Index(nbStates_));

  hmmEq_ = ConfiguredParametrizable::createVector<ConfiguredTransitionMatrix, EquilibriumFrequenciesFromTransitionMatrix, Eigen::VectorXd>(context_, {matrix_}, eqdim);

  // transition

  auto transdim = MatrixDimension(nbStates_, nbStates_);
  hmmTrans_ = ConfiguredParametrizable::createMatrix<ConfiguredTransitionMatrix, TransitionMatrixFromTransitionMatrix, Eigen::MatrixXd>(context_, {matrix_}, transdim);

  // emission

  hmmEmis_ = emissionProbabilities_->getEmissionProbabilities();

  // Manage parameters:
  addParameters_(hiddenAlphabet_->getParameters());
  addParameters_(emissionProbabilities_->getParameters());

  // forward computation
  forwardLik_ = ForwardHmmLikelihood_DF::create(context_, {hmmEq_, hmmTrans_, hmmEmis_}, MatrixDimension(nbStates_, nbSites_));

  setLikelihoodNode(SumOfLogarithms<RowLik>::create (getContext_(), {forwardLik_}, RowVectorDimension (nbSites_)));

  // backward computation
  backwardLik_ = BackwardHmmLikelihood_DF::create(context_, {forwardLik_, hmmTrans_, hmmEmis_}, MatrixDimension(nbStates_, nbSites_));


  // Then Hidden Posterior Probabilities

  auto forwardNode = dynamic_pointer_cast<ForwardHmmLikelihood_DF>(forwardLik_);

  hiddenPostProb_ = CWiseMul<MatrixLik, std::tuple<MatrixLik, MatrixLik> >::create(context_, {forwardNode->getForwardCondLikelihood(), backwardLik_}, MatrixDimension(nbStates_, nbSites_));

  // and site likelihoods

  auto statesLog = CWiseMul<MatrixLik, std::tuple<MatrixLik, MatrixLik> >::create(context_, {hiddenPostProb_, hmmEmis_}, MatrixDimension(nbStates_, nbSites_));

  setSiteLikelihoods(CWiseAdd<RowLik, MatrixLik>::create(context_, {statesLog}, RowVectorDimension(nbSites_)));
}

void HmmLikelihood_DF::setNamespace(const std::string& nameSpace)
{
  AbstractParametrizable::setNamespace(nameSpace);

  hiddenAlphabet_->setNamespace(nameSpace);
  matrix_->setNamespace(nameSpace);
  emissionProbabilities_->setNamespace(nameSpace);
}
