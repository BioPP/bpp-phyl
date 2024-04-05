// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Stat/Mva/CorrespondenceAnalysis.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Seq/Container/SequenceContainer.h>
#include <Bpp/Seq/ProbabilisticSequence.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>

#include "CoalaCore.h"

using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/******************************************************************************/

CoalaCore::CoalaCore(size_t nbAxes) :
  init_(true),
  nbrOfAxes_(nbAxes),
  P_(),
  R_(),
  colWeights_(),
  paramValues_()
{}

/******************************************************************************/

ParameterList CoalaCore::computeCOA(const SequenceDataInterface& data, double pseudoCount, bool param)
{
  ParameterList pList;
  // Now we perform the Correspondence Analysis on from the matrix of observed frequencies computed on the alignment, to obtain the matrix of principal axes.
  // First, the matrix of amino acid frequencies is calculated from the alignment:
  vector<string> seqKeys = data.getSequenceKeys();
  vector< map<int, double>> freqs(seqKeys.size()); // One map per sequence
  // Each map is filled with the corresponding frequencies, which are then normalized.
  for (size_t i = 0; i < seqKeys.size(); ++i)
  {
    const SequenceContainerInterface* sc = dynamic_cast<const SequenceContainerInterface*>(&data);
    const ProbabilisticSequenceContainerInterface* psc = dynamic_cast<const ProbabilisticSequenceContainerInterface*>(&data);

    shared_ptr<CruxSymbolListInterface> seq(sc ?
        dynamic_cast<CruxSymbolListInterface*>(sc->sequence(seqKeys[i]).clone()) :
        dynamic_cast<CruxSymbolListInterface*>(psc->sequence(seqKeys[i]).clone()));

    SymbolListTools::changeGapsToUnknownCharacters(*seq);
    std::map<int, double> counts;
    SequenceTools::getCounts(*seq, counts);

    // add pseudoCounts
    for (int k = 0; k < 20; ++k)
    {
      counts[k] += pseudoCount;
    }
    
    // Unknown characters are now ignored:
    double t = 0;
    for (int k = 0; k < 20; ++k)
    {
      t += counts[k];
    }
    
    for (int k = 0; k < 20; ++k)
    {
      freqs.at(i)[k] = counts[k]/t;
    }
  }

  // The matrix of observed frequencies is filled. If an amino acid is
  // completely absent from the alignment, its frequency is set to
  // 10^-6.
  RowMatrix<double> freqMatrix(seqKeys.size(), 20);
  for (size_t i = 0; i < freqs.size(); ++i)
  {
    bool normalize = false;
    for (size_t j = 0; j < 20; ++j)
    {
      map<int, double>::iterator it = freqs[i].find(static_cast<int>(j));
      if (it != freqs[i].end())
      {
        freqMatrix(i, j) = (*it).second;
      }
      else
      {
        freqMatrix(i, j) = 0.000001;
        normalize = true;
      }
    }
    if (normalize)
    {
      double sum = 0;
      for (size_t k = 0; k < 20; ++k)
      {
        sum += freqMatrix(i, k);
      }
      for (size_t l = 0; l < 20; ++l)
      {
        freqMatrix(i, l) = freqMatrix(i, l) / sum;
      }
    }
  }

  // The COA analysis:
  CorrespondenceAnalysis coa(freqMatrix, 19);
  // Matrix of principal axes:
  RowMatrix<double> ppalAxes = coa.getPrincipalAxes();
  // The transpose of the matrix of principal axes is computed:
  MatrixTools::transpose(ppalAxes, P_);
  // The matrix of row coordinates is stored:
  R_ = coa.getRowCoordinates();
  // The column weights are retrieved:
  colWeights_ = coa.getColumnWeights();

  if (param)
  {
    // Parameters are defined:
    size_t nbAxesConserved = coa.getNbOfKeptAxes();
    if (nbrOfAxes_ > nbAxesConserved)
    {
      ApplicationTools::displayWarning("The specified number of parameters per branch (" + TextTools::toString(nbrOfAxes_) +
          ") is higher than the number of axes (" + TextTools::toString(nbAxesConserved) +
          ")... The number of parameters per branch is now equal to the number of axes kept by the COA analysis (" + TextTools::toString(nbAxesConserved) + ")");
      nbrOfAxes_ = nbAxesConserved;
    }
    for (unsigned int i = 0; i < nbrOfAxes_; ++i)
    {
      const vector<double> rCoords = R_.col(i);
      double maxCoord = VectorTools::max(rCoords);
      double minCoord = VectorTools::min(rCoords);
      double sd = VectorTools::sd<double, double>(rCoords);
      auto constraint = make_shared<IntervalConstraint>(minCoord - sd, maxCoord + sd, true, true);
      if (paramValues_.hasParameter("AxPos" + TextTools::toString(i)))
        pList.addParameter(new Parameter("Coala.AxPos" + TextTools::toString(i), paramValues_.getParameterValue("AxPos" + TextTools::toString(i).substr(0, 8)), constraint));
      else
        pList.addParameter(new Parameter("Coala.AxPos" + TextTools::toString(i), 0., constraint));
    }
  }
  return pList;
}

/******************************************************************************/
/* Function that computes the product of a matrix P of size nbrOfAxes_x20 with a vector V of size nbrOfAxes_, and returns a vector of size 20.*/

vector<double> CoalaCore::prodMatrixVector(RowMatrix<double>& P, vector<double>& V)
{
  vector<double> E(20, 0.0);

  for (unsigned int i = 0; i < 20; i++)
  {
    for (unsigned int j = 0; j < V.size(); j++)
    {
      E[i] = E[i] + P(j, i) * V[j];
    }
  }
  return E;
}

/******************************************************************************/
