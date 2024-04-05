// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Stat/Mva/CorrespondenceAnalysis.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>

#include "Coala.h"

using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>
#include <string>

using namespace std;


/******************************************************************************/

Coala::Coala(
  shared_ptr<const ProteicAlphabet> alpha,
  const ProteinSubstitutionModelInterface& model,
  unsigned int nbAxes,
  bool param) :
  AbstractParameterAliasable("Coala."),
  AbstractReversibleProteinSubstitutionModel(alpha, model.getStateMap(), "Coala."),
  CoalaCore(nbAxes),
  init_(true),
  nbrOfAxes_(nbAxes),
  file_(),
  param_(param)
{
  setNamespace(getName() + ".");

  // Setting the exchangeability matrix
  exchangeability_ = model.exchangeabilityMatrix();
  updateMatrices_();
}

/******************************************************************************/

void Coala::readFromFile(string& file)
{
  ifstream in(file.c_str(), ios::in);
  // Read exchangeability matrix:
  for (unsigned int i = 1; i < 20; ++i)
  {
    string line = FileTools::getNextLine(in);
    StringTokenizer st(line);
    for (unsigned int j = 0; j < i; ++j)
    {
      double s = TextTools::toDouble(st.nextToken());
      exchangeability_(i, j) = exchangeability_(j, i) = s;
    }
  }

  // Now build diagonal of the exchangeability matrix:
  for (unsigned int i = 0; i < 20; ++i)
  {
    double sum = 0;
    for (unsigned int j = 0; j < 20; ++j)
    {
      if (j != i)
        sum += exchangeability_(i, j);
    }
    exchangeability_(i, i) = -sum;
  }

  // Closing stream:
  in.close();
}


/******************************************************************************/

void Coala::computeEquilibriumFrequencies()
{
  // Computes the equilibrium frequencies from a set of coordinates along the principal axes of the COA.
  if (init_)
    init_ = false;
  else
  {
    // We get the coordinates:
    vector<double> coord;
    for (unsigned int i = 0; i < nbrOfAxes_; ++i)
    {
      coord.push_back(parameter("AxPos" + TextTools::toString(i)).getValue());
    }

    // Now, frequencies are computed from the vector of coordinates and the transpose of the principal axes matrix (P_):
    vector<double> tmpFreqs;
    tmpFreqs = prodMatrixVector(P_, coord);
    for (unsigned int i = 0; i < tmpFreqs.size(); ++i)
    {
      tmpFreqs[i] = (tmpFreqs[i] + 1) * colWeights_[i];
    }
    freq_ = tmpFreqs;

    // Frequencies are not allowed to be lower than 10^-3 or higher than 0.5:
    bool norm = false;
    for (unsigned int i = 0; i < 20; ++i)
    {
      if (freq_[i] < 0.001)
      {
        norm = true;
        freq_[i] = 0.001;
      }
      if (freq_[i] > 0.2)
      {
        norm = true;
        freq_[i] = 0.2;
      }
    }
    if (norm == true)
    {
      double s = VectorTools::sum(freq_);
      for (size_t i = 0; i < 20; i++)
      {
        freq_[i] = freq_[i] / s;
      }
    }
  }
}

/******************************************************************************/

void Coala::updateMatrices_()
{
  computeEquilibriumFrequencies();
  AbstractReversibleSubstitutionModel::updateMatrices_();
}

/******************************************************************************/

void Coala::setFreqFromData(const SequenceDataInterface& data, double pseudoCount)
{
  // Compute the COA from the observed frequencies, add the axis position parameters and update the Markov matrix
  auto plist = computeCOA(data, pseudoCount, param_);
  addParameters_(plist);
  updateMatrices_();
}

/******************************************************************************/
