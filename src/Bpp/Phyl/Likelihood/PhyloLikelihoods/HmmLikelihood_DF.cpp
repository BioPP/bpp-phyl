// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "HmmLikelihoodComputation.h"
#include "HmmLikelihood_DF.h"

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

  hiddenPostProb_ = CWiseMul<Eigen::MatrixXd, std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>>::create(context_, {forwardNode->getForwardCondLikelihood(), backwardLik_}, MatrixDimension(nbStates_, nbSites_));

  // and site likelihoods

  auto statesLog = CWiseMul<MatrixLik, std::tuple<Eigen::MatrixXd, MatrixLik>>::create(context_, {hiddenPostProb_, hmmEmis_}, MatrixDimension(nbStates_, nbSites_));

  setSiteLikelihoods(CWiseAdd<RowLik, MatrixLik>::create(context_, {statesLog}, RowVectorDimension(nbSites_)));
}

void HmmLikelihood_DF::setNamespace(const std::string& nameSpace)
{
  AbstractParametrizable::setNamespace(nameSpace);

  hiddenAlphabet_->setNamespace(nameSpace);
  matrix_->setNamespace(nameSpace);
  emissionProbabilities_->setNamespace(nameSpace);
}
