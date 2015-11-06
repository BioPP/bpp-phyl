//
// File: MixtureSequenceEvolution.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 30 avril 2015, à 17h 20
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

#include "MixtureSequenceEvolution.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

MixtureSequenceEvolution::MixtureSequenceEvolution(
  SubstitutionProcessCollection* processColl,
  std::vector<size_t>& nProc) :
  MultiProcessSequenceEvolution(processColl, nProc, "Mixture."),
  simplex_(nProc.size(), 1, false, "Mixture.")
{
  // initialize parameters:
  addParameters_(simplex_.getParameters());
}

void MixtureSequenceEvolution::setSubProcessProb(const Simplex& si)
{
  simplex_.setFrequencies(si.getFrequencies());
  matchParametersValues(simplex_.getParameters());
}

void MixtureSequenceEvolution::setNamespace(const std::string& nameSpace)
{
  deleteParameters_(simplex_.getParameters().getParameterNames());

  simplex_.setNamespace(nameSpace);

  addParameters_(simplex_.getParameters());
}


void MixtureSequenceEvolution::fireParameterChanged(const ParameterList& parameters)
{
  MultiProcessSequenceEvolution::fireParameterChanged(parameters);
  simplex_.matchParametersValues(parameters);
}

ParameterList MixtureSequenceEvolution::getNonDerivableParameters() const
{
  // patch, to be fixed properly later
  return getIndependentParameters();
  
  ParameterList pl = MultiProcessSequenceEvolution::getNonDerivableParameters();
  pl.addParameters(simplex_.getParameters());
  
  return pl;
}

