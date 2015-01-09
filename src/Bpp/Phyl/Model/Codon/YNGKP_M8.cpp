//
// File: YNGKP_M8.cpp
// Created by:  Laurent Gueguen
// Created on: May 2010
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

#include "YNGKP_M8.h"
#include "YN98.h"

#include <Bpp/Numeric/Prob/MixtureOfDiscreteDistributions.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/BetaDiscreteDistribution.h>
#include <Bpp/Text/TextTools.h>

using namespace bpp;

using namespace std;

/******************************************************************************/

YNGKP_M8::YNGKP_M8(const GeneticCode* gc, FrequenciesSet* codonFreqs, unsigned int nclass) :
  AbstractBiblioMixedSubstitutionModel("YNGKP_M8."),
  pmixmodel_(0),
  synfrom_(),
  synto_()
{
  if (nclass <= 0)
    throw Exception("Bad number of classes for model YNGKP_M8: " + TextTools::toString(nclass));

  // build the submodel

  BetaDiscreteDistribution* pbdd = new BetaDiscreteDistribution(nclass, 2, 2);

  vector<double> val; val.push_back(2);
  vector<double> prob; prob.push_back(1);
  SimpleDiscreteDistribution* psdd = new SimpleDiscreteDistribution(val, prob);

  vector<DiscreteDistribution*> v_distr;
  v_distr.push_back(pbdd); v_distr.push_back(psdd);
  prob.clear(); prob.push_back(0.5); prob.push_back(0.5);

  MixtureOfDiscreteDistributions* pmodd = new MixtureOfDiscreteDistributions(v_distr, prob);

  map<string, DiscreteDistribution*> mpdd;
  mpdd["omega"] = pmodd;

  YN98* yn98 = new YN98(gc, codonFreqs);
  yn98->setNamespace("YNGKP_M8.");

  pmixmodel_.reset(new MixtureOfASubstitutionModel(gc->getSourceAlphabet(), yn98, mpdd));
  delete pbdd;

  pmixmodel_->setNamespace("YNGKP_M8.");

  vector<int> supportedChars = yn98->getAlphabetStates();
  
  // mapping the parameters

  ParameterList pl = pmixmodel_->getParameters();
  for (size_t i = 0; i < pl.size(); i++)
  {
    lParPmodel_.addParameter(Parameter(pl[i]));
  }

  vector<std::string> v = dynamic_cast<YN98*>(pmixmodel_->getNModel(0))->getFrequenciesSet()->getParameters().getParameterNames();

  for (size_t i = 0; i < v.size(); i++)
  {
    mapParNamesFromPmodel_[v[i]] = getParameterNameWithoutNamespace(v[i]);
  }

  mapParNamesFromPmodel_["YNGKP_M8.kappa"] = "kappa";
  mapParNamesFromPmodel_["YNGKP_M8.omega_Mixture.theta1"] = "p0";
  mapParNamesFromPmodel_["YNGKP_M8.omega_Mixture.1_Beta.alpha"] = "p";
  mapParNamesFromPmodel_["YNGKP_M8.omega_Mixture.1_Beta.beta"] = "q";
  mapParNamesFromPmodel_["YNGKP_M8.omega_Mixture.2_Simple.V1"] = "omegas";

  // specific parameters

  string st;
  for (map<string, string>::iterator it = mapParNamesFromPmodel_.begin(); it != mapParNamesFromPmodel_.end(); it++)
  {
    st = pmixmodel_->getParameterNameWithoutNamespace(it->first);
    if (it->second != "omegas")
      addParameter_(new Parameter("YNGKP_M8." + it->second, pmixmodel_->getParameterValue(st),
                              pmixmodel_->getParameter(st).hasConstraint() ? pmixmodel_->getParameter(st).getConstraint()->clone() : 0, true));
  }

  addParameter_(new Parameter("YNGKP_M8.omegas", 2., new IntervalConstraint(1, 1, false), true));

  // look for synonymous codons
  for (synfrom_ = 1; synfrom_ < supportedChars.size(); synfrom_++)
  {
    for (synto_ = 0; synto_ < synfrom_; synto_++)
    {
      if ((gc->areSynonymous(supportedChars[synfrom_], supportedChars[synto_]))
          && (pmixmodel_->getNModel(0)->Qij(synfrom_, synto_) != 0)
          && (pmixmodel_->getNModel(1)->Qij(synfrom_, synto_) != 0))
        break;
    }
    if (synto_ < synfrom_)
      break;
  }

  if (synto_ == gc->getSourceAlphabet()->getSize())
    throw Exception("Impossible to find synonymous codons");

  // update Matrices
  updateMatrices();
}

YNGKP_M8::YNGKP_M8(const YNGKP_M8& mod2) : AbstractBiblioMixedSubstitutionModel(mod2),
  pmixmodel_(new MixtureOfASubstitutionModel(*mod2.pmixmodel_)),
  synfrom_(mod2.synfrom_),
  synto_(mod2.synto_)
{}

YNGKP_M8& YNGKP_M8::operator=(const YNGKP_M8& mod2)
{
  AbstractBiblioMixedSubstitutionModel::operator=(mod2);

  pmixmodel_.reset(new MixtureOfASubstitutionModel(*mod2.pmixmodel_));
  synfrom_ = mod2.synfrom_;
  synto_ = mod2.synto_;

  return *this;
}

YNGKP_M8::~YNGKP_M8() {}

void YNGKP_M8::updateMatrices()
{
  AbstractBiblioSubstitutionModel::updateMatrices();

  // homogeneization of the synonymous substittion rates

  Vdouble vd;

  for (unsigned int i = 0; i < pmixmodel_->getNumberOfModels(); i++)
  {
    vd.push_back(1 / pmixmodel_->getNModel(i)->Qij(synfrom_, synto_));
  }

  pmixmodel_->setVRates(vd);
}

