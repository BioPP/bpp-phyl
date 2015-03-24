//
// File: Coala.cpp
// Created by: Mathieu Groussin
// Created on: Sun Mar 13 12:00:00 2011
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


#include "Coala.h"

#include <Bpp/Io/FileTools.h>
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

Coala::Coala(
  const ProteicAlphabet* alpha,
  const ProteinSubstitutionModel& model,  
  unsigned int nbAxes,
  bool param) :
  AbstractParameterAliasable("Coala."),
  AbstractReversibleProteinSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), "Coala."),
  CoalaCore(nbAxes, model.getName()),
  init_(true),
  nbrOfAxes_(nbAxes),
  exch_(model.getName()),
  file_(),
  param_(param)
{
  setNamespace(getName() + ".");

  // Setting the exchangeability matrix
  exchangeability_ = model.getExchangeabilityMatrix();
  updateMatrices();
}

/******************************************************************************/

void Coala::readFromFile(string& file)
{
  ifstream in(file.c_str(), ios::in);
  // Read exchangeability matrix:
  for (unsigned int i = 1; i < 20; i++)
  {
    string line = FileTools::getNextLine(in);
    StringTokenizer st(line);
    for (unsigned int j = 0; j < i; j++)
    {
      double s = TextTools::toDouble(st.nextToken());
      exchangeability_(i, j) = exchangeability_(j, i) = s;
    }
  }

  // Now build diagonal of the exchangeability matrix:
  for (unsigned int i = 0; i < 20; i++)
  {
    double sum = 0;
    for (unsigned int j = 0; j < 20; j++)
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
    for (unsigned int i = 0; i < nbrOfAxes_; i++)
    {
      coord.push_back(getParameter("AxPos" + TextTools::toString(i)).getValue());
    }

    // Now, frequencies are computed from the vector of coordinates and the transpose of the principal axes matrix (P_):
    vector<double> tmpFreqs;
    tmpFreqs = prodMatrixVector(P_, coord);
    for (unsigned int i = 0; i < tmpFreqs.size(); i++)
    {
      tmpFreqs[i] = (tmpFreqs[i] + 1) * colWeights_[i];
    }
    freq_ = tmpFreqs;

    // Frequencies are not allowed to be lower than 10^-3 or higher than 0.5:
    bool norm = false;
    for (unsigned int i = 0; i < 20; i++)
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

void Coala::updateMatrices()
{
  computeEquilibriumFrequencies();
  AbstractReversibleSubstitutionModel::updateMatrices();
}

/******************************************************************************/

void Coala::setFreqFromData(const SequenceContainer& data, double pseudoCount)
{
  // Compute the COA from the observed frequencies, add the axis position parameters and update the Markov matrix
  ParameterList pList = computeCOA(data, param_);
  addParameters_(pList);
  updateMatrices();
}

/******************************************************************************/
