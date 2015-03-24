//
// File: CoalaCore.cpp
// Created by: Mathieu Groussin
// Created on: Sun Mar 13 12:00:00 2011
//

/*
   Copyright or (C) or Copr. Bio++ Development Team, (November 16, 2004)

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


#include "CoalaCore.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Stat/Mva/CorrespondenceAnalysis.h>

#include <Bpp/Seq/SequenceTools.h>


using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/******************************************************************************/

CoalaCore::CoalaCore(size_t nbAxes, const string& exch) :
  init_(true),
  nbrOfAxes_(nbAxes),
  exch_(exch),
  P_(),
  R_(),
  colWeights_(),
  paramValues_()
{}

/******************************************************************************/

ParameterList CoalaCore::computeCOA(const SequenceContainer& data, bool param)
{
  ParameterList pList;
  // Now we perform the Correspondence Analysis on from the matrix of observed frequencies computed on the alignment, to obtain the matrix of principal axes.
  // First, the matrix of amino acid frequencies is calculated from the alignment:
  vector<string> names = data.getSequencesNames();
  vector< map<int, double> > freqs(names.size()); // One map per sequence
  // Each map is filled with the corresponding frequencies, which are then normalized.
  for (size_t i = 0; i < names.size(); ++i)
  {
    Sequence* seq = new BasicSequence(data.getSequence(names[i]));
    SymbolListTools::changeGapsToUnknownCharacters(*seq);
    SequenceTools::getFrequencies(*seq, freqs.at(i));
    // Unknown characters are now ignored:
    double t = 0;
    for (int k = 0; k < 20; ++k)
    {
      t += freqs.at(i)[k];
    }
    for (int k = 0; k < 20; k++)
    {
      freqs.at(i)[k] /= t;
    }
    delete seq;
  }

  // The matrix of observed frequencies is filled. If an amino acid is completely absent from the alignment, its frequency is set to 10^-6.
  RowMatrix<double> freqMatrix(names.size(), 20);
  for (size_t i = 0; i < freqs.size(); i++)
  {
    bool normalize = false;
    for (size_t j = 0; j < 20; j++)
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
      for (size_t k = 0; k < 20; k++)
      {
        sum += freqMatrix(i, k);
      }
      for (size_t l = 0; l < 20; l++)
      {
        freqMatrix(i, l) = freqMatrix(i, l) / sum;
      }
    }
  }

  // The COA analysis:
  CorrespondenceAnalysis* coa = new CorrespondenceAnalysis(freqMatrix, 19);
  // Matrix of principal axes:
  RowMatrix<double> ppalAxes = coa->getPrincipalAxes();
  // The transpose of the matrix of principal axes is computed:
  MatrixTools::transpose(ppalAxes, P_);
  // The matrix of row coordinates is stored:
  R_ = coa->getRowCoordinates();
  // The column weights are retrieved:
  colWeights_ = coa->getColumnWeights();

  if (param)
  {
    // Parameters are defined:
    size_t nbAxesConserved = coa->getNbOfKeptAxes();
    if (nbrOfAxes_ > nbAxesConserved)
    {
      ApplicationTools::displayWarning("The specified number of parameters per branch (" + TextTools::toString(nbrOfAxes_) +
                                       ") is higher than the number of axes (" + TextTools::toString(nbAxesConserved) +
                                       ")... The number of parameters per branch is now equal to the number of axes kept by the COA analysis (" + TextTools::toString(nbAxesConserved) + ")");
      nbrOfAxes_ = nbAxesConserved;
    }
    for (unsigned int i = 0; i < nbrOfAxes_; i++)
    {
      const vector<double> rCoords = R_.col(i);
      double maxCoord = VectorTools::max(rCoords);
      double minCoord = VectorTools::min(rCoords);
      double sd = VectorTools::sd<double, double>(rCoords);
      IntervalConstraint* constraint = new IntervalConstraint(minCoord - sd, maxCoord + sd, true, true);
      if (paramValues_.find("AxPos" + TextTools::toString(i)) != paramValues_.end())
        pList.addParameter(Parameter("Coala.AxPos" + TextTools::toString(i), TextTools::toDouble(paramValues_["AxPos" + TextTools::toString(i)].substr(0, 8)), constraint));
      else
        pList.addParameter(Parameter("Coala.AxPos" + TextTools::toString(i), 0., constraint));
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
