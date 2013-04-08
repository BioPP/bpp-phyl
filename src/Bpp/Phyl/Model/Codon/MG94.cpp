//
// File: MG94.cpp
// Created by:  Laurent Gueguen
// Created on: July 2009
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)
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

#include "MG94.h"
#include "../FrequenciesSet/CodonFrequenciesSet.h"
#include "../Nucleotide/K80.h"

using namespace bpp;

using namespace std;

/******************************************************************************/

MG94::MG94(const GeneticCode* gc, FrequenciesSet* codonFreqs) :
  AbstractBiblioSubstitutionModel("MG94."),
  pmodel_(new CodonDistancePhaseFrequenciesSubstitutionModel(gc, new K80(dynamic_cast<const CodonAlphabet*>(gc->getSourceAlphabet())->getNucleicAlphabet()), codonFreqs))
{
  addParameter_(new Parameter("MG94.rho", 1, &Parameter::R_PLUS_STAR));

  pmodel_->setNamespace("MG94.");
  addParameters_(codonFreqs->getParameters());

  lParPmodel_.addParameters(pmodel_->getParameters());
  
  vector<std::string> v=pmodel_->getFrequenciesSet()->getParameters().getParameterNames();
  for (unsigned int i=0;i<v.size();i++)
    mapParNamesFromPmodel_[v[i]]=getParameterNameWithoutNamespace(v[i]);

  mapParNamesFromPmodel_["MG94.beta"]="rho";
  
  updateMatrices();
}

MG94::MG94(const MG94& mg94) :
  AbstractBiblioSubstitutionModel(mg94),
  pmodel_(new CodonDistancePhaseFrequenciesSubstitutionModel(*mg94.pmodel_))
{}

MG94& MG94::operator=(const MG94& mg94)
{
  AbstractBiblioSubstitutionModel::operator=(mg94);
  pmodel_.reset(new CodonDistancePhaseFrequenciesSubstitutionModel(*mg94.pmodel_));
  return *this;
}

MG94::~MG94() {}


