//
// File: TripletSubstitutionModel.cpp
// Authors:
//   Laurent Gueguen
// Created: 2009-02-08 00:00:00
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

#include <Bpp/Numeric/Matrix/EigenValue.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>

#include "TripletSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/WordAlphabet.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/******************************************************************************/

TripletSubstitutionModel::TripletSubstitutionModel(
  shared_ptr<const CodonAlphabet> palph,
  unique_ptr<NucleotideSubstitutionModelInterface> pmod) :
  AbstractParameterAliasable("Triplet."),
  WordSubstitutionModel(palph, make_shared<CanonicalStateMap>(palph, false), "Triplet.")
{
  unsigned int i;
  addParameters_(pmod->getParameters());

  Vrate_.resize(3);
  for (i = 0; i < 3; i++)
  {
    VSubMod_.push_back(std::move(pmod));
    VnestedPrefix_.push_back(VSubMod_[i]->getNamespace());
    Vrate_[i] = 1. / 3.;
  }

  // relative rates
  for (i = 0; i < 2; i++)
  {
    addParameter_(new Parameter("Triplet.relrate" + TextTools::toString(i + 1), 1.0 / (3 - i), Parameter::PROP_CONSTRAINT_EX));
  }

  WordSubstitutionModel::updateMatrices_();
}

TripletSubstitutionModel::TripletSubstitutionModel(
  shared_ptr<const CodonAlphabet> palph,
  unique_ptr<NucleotideSubstitutionModelInterface> pmod1,
  unique_ptr<NucleotideSubstitutionModelInterface> pmod2,
  unique_ptr<NucleotideSubstitutionModelInterface> pmod3) :
  AbstractParameterAliasable("Triplet."),
  WordSubstitutionModel(palph, make_shared<CanonicalStateMap>(palph, false), "Triplet.")
{
  string st = "Triplet.";

  VSubMod_.push_back(std::move(pmod1));
  VnestedPrefix_.push_back(VSubMod_[0]->getNamespace());
  VSubMod_[0]->setNamespace(st + "1_" + VnestedPrefix_[0]);
  addParameters_(VSubMod_[0]->getParameters());

  VSubMod_.push_back(std::move(pmod2));
  VnestedPrefix_.push_back(VSubMod_[1]->getNamespace());
  VSubMod_[1]->setNamespace(st + "2_" + VnestedPrefix_[1]);
  addParameters_(VSubMod_[1]->getParameters());

  VSubMod_.push_back(std::move(pmod3));
  VnestedPrefix_.push_back(VSubMod_[2]->getNamespace());
  VSubMod_[2]->setNamespace(st + "3_" + VnestedPrefix_[2]);
  addParameters_(VSubMod_[2]->getParameters());

  Vrate_.resize(3);
  for (unsigned int i = 0; i < 3; ++i)
  {
    Vrate_[i] = 1.0 / 3;
  }

  // relative rates
  for (unsigned int i = 0; i < 2; ++i)
  {
    addParameter_(new Parameter(st + "relrate" + TextTools::toString(i + 1), 1.0 / (3 - i), Parameter::PROP_CONSTRAINT_EX));
  }

  WordSubstitutionModel::updateMatrices_();
}

string TripletSubstitutionModel::getName() const
{
  string s = "TripletSubstitutionModel model:";
  for (auto& vi :  VSubMod_)
  {
    s += " " + vi->getName();
  }

  return s;
}

