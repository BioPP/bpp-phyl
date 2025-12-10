// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "MvaFrequencySet.h"

using namespace bpp;

#include <cmath>
using namespace std;

MvaFrequencySet::MvaFrequencySet(shared_ptr<const ProteicAlphabet> alpha) :
  AbstractFrequencySet(
      make_shared<CanonicalStateMap>(alpha, false),
      "MVA.",
      "MVAprotein"),
  model_(),
  tPpalAxes_(),
  rowCoords_(),
  nbrOfAxes_(0),
  columnWeights_(),
  paramValues_()
{}

void MvaFrequencySet::initSet(std::shared_ptr<const Coala> coala)
{
  model_ = coala;
  setNbrOfAxes(coala->getNbrOfAxes());
  setTransposeMatrixOfPpalAxes(coala->getTppalAxesMatrix());
  setMatrixOfRowCoords(coala->getRowCoordinates());
  setVectorOfColumnWeights(coala->getColumnWeights());
  defineParameters();
  updateFrequencies();
}

void MvaFrequencySet::defineParameters()
{
  for (unsigned int i = 0; i < nbrOfAxes_; i++)
  {
    const vector<double> rCoords = rowCoords_.col(i);
    double maxCoord = VectorTools::max(rCoords);
    double minCoord = VectorTools::min(rCoords);
    double sd = VectorTools::sd<double, double>(rCoords);
    auto constraint = std::make_shared<IntervalConstraint>(minCoord - sd, maxCoord + sd, true, true);
    if (paramValues_.find("RootAxPos" + TextTools::toString(i)) != paramValues_.end())
      addParameter_(new Parameter(getNamespace() + "RootAxPos" + TextTools::toString(i), TextTools::toDouble(paramValues_["RootAxPos" + TextTools::toString(i)].substr(0, 8)), constraint));
    else
      addParameter_(new Parameter(getNamespace() + "RootAxPos" + TextTools::toString(i), 0., constraint));
  }
}

void MvaFrequencySet::fireParameterChanged(const ParameterList& parameters)
{
  updateFrequencies();
}

void MvaFrequencySet::updateFrequencies()
{
  if (nbrOfAxes_ == 0)
    throw Exception("The number of axes kept by the MVA analysis was not set. You should initialize it with the setNbrOfAxes function");
  vector<double> positions;

  for (unsigned int i = 0; i < nbrOfAxes_; i++)
  {
    positions.push_back(getParameterValue("RootAxPos" + TextTools::toString(i)));
  }

  vector<double> tmpFreqs(20, 0.0);
  vector<double> freqs(20, 0.0);

  computeReverseCOA(positions, tmpFreqs);
  computeCoordsFirstSpaceCOA(tmpFreqs, freqs);

  setFrequencies_(freqs);

  bool norm = false;
  for (unsigned int i = 0; i < 20; i++)
  {
    if (getFreq_(i) < 0.001)
    {
      norm = true;
      getFreq_(i) = 0.001;
    }
    if (getFreq_(i) > 0.5)
    {
      norm = true;
      getFreq_(i) = 0.5;
    }
  }
  if (norm == true)
  {
    double s = VectorTools::sum(getFrequencies());
    for (size_t i = 0; i < 20; ++i)
    {
      getFreq_(i) = getFreq_(i) / s;
    }
  }
}

void MvaFrequencySet::setFrequencies(const vector<double>& frequencies)
{}

void MvaFrequencySet::computeReverseCOA(const std::vector<double>& positions, std::vector<double>& tmpFreqs)
{
  for (unsigned int i = 0; i < 20; i++)
  {
    for (unsigned int j = 0; j < nbrOfAxes_; j++)
    {
      tmpFreqs[i] = tmpFreqs[i] + tPpalAxes_(j, i) * positions[j];
    }
  }
}

void MvaFrequencySet::computeCoordsFirstSpaceCOA(std::vector<double>& tmpFreqs, std::vector<double>& freqs)
{
  if (freqs.size() != tmpFreqs.size())
    throw Exception("MvaFrequencySet::computeCoordsFirstSpaceCOA : error in the size of the vectors");
  // The vector of amino acid frequencies is calculated from the original column weights
  for (unsigned int i = 0; i < tmpFreqs.size(); i++)
  {
    freqs[i] = (tmpFreqs[i] + 1) * columnWeights_[i];
  }
}
