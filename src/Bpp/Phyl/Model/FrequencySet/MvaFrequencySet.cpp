//
// File: MvaFrequencySet.cpp
// Authors:
//   Mathieu Groussin
// Created: 2013-01-12 00:00:00
//

/*
  Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)
  
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


#include "MvaFrequencySet.h"

using namespace bpp;

#include <cmath>
using namespace std;

MvaFrequencySet::MvaFrequencySet(shared_ptr<const ProteicAlphabet> alpha) :
  AbstractFrequencySet(
      make_shared<CanonicalStateMap>(alpha, false), 
      "MVA.",
      "MVAprotein"),
  tPpalAxes_(),
  rowCoords_(),
  nbrOfAxes_(0),
  model_(),
  columnWeights_(),
  paramValues_()
{}

void MvaFrequencySet::initSet(const CoalaCore& coala)
{
  setNbrOfAxes(coala.getNbrOfAxes());
  setTransposeMatrixOfPpalAxes(coala.getTppalAxesMatrix());
  setMatrixOfRowCoords(coala.getRowCoordinates());
  setVectorOfColumnWeights(coala.getColumnWeights());
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
    std::shared_ptr<Constraint> constraint = std::make_shared<IntervalConstraint>(minCoord - sd, maxCoord + sd, true, true);
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
