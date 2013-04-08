//
// File: YNGKP_M3.cpp
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

#include "YNGKP_M3.h"
#include "YN98.h"

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

using namespace bpp;

using namespace std;

/******************************************************************************/

YNGKP_M3::YNGKP_M3(const GeneticCode* gc, FrequenciesSet* codonFreqs, unsigned int nbOmega) :
  AbstractBiblioMixedSubstitutionModel("YNGKP_M3."),
  pmixmodel_(0),
  synfrom_(-1),
  synto_(-1)
{
  if (nbOmega < 1)
    throw Exception("At least one omega is necessary in the YNGKP_M3 model");

  // build the submodel

  vector<double> v1, v2;
  v1.push_back(0.5);
  for (unsigned int i = 1; i < nbOmega; i++)
  {
    v1.push_back(0.5 + 0.5 * i);
  }

  for (unsigned int i = 0; i < nbOmega; i++)
  {
    v2.push_back(1. / nbOmega);
  }

  SimpleDiscreteDistribution* psdd = new SimpleDiscreteDistribution(v1, v2, NumConstants::MILLI());

  map<string, DiscreteDistribution*> mpdd;
  mpdd["omega"] = psdd;

  pmixmodel_.reset(new MixtureOfASubstitutionModel(gc->getSourceAlphabet(), new YN98(gc, codonFreqs), mpdd));
  delete psdd;

  // mapping the parameters

  ParameterList pl = pmixmodel_->getParameters();
  for (size_t i = 0; i < pl.size(); i++)
  {
    lParPmodel_.addParameter(Parameter(pl[i]));
  }

  vector<std::string> v = dynamic_cast<YN98*>(pmixmodel_->getNModel(0))->getFrequenciesSet()->getParameters().getParameterNames();
  for (size_t i = 0; i < v.size(); i++)
  {
    mapParNamesFromPmodel_[v[i]] = getParameterNameWithoutNamespace("YNGKP_M3." + v[i].substr(5));
  }

  mapParNamesFromPmodel_["YN98.kappa"] = "kappa";

  for (unsigned int i = 1; i < nbOmega; i++)
  {
    mapParNamesFromPmodel_["YN98.omega_Simple.theta" + TextTools::toString(i)] = "theta" + TextTools::toString(i);
  }


  mapParNamesFromPmodel_["YN98.omega_Simple.V1"] = "omega0";
  for (unsigned int i = 1; i < nbOmega; i++)
  {
    mapParNamesFromPmodel_["YN98.omega_Simple.V" + TextTools::toString(i + 1)] = "delta" + TextTools::toString(i);
  }

  // specific parameters

  string st;
  for (map<string, string>::iterator it = mapParNamesFromPmodel_.begin(); it != mapParNamesFromPmodel_.end(); it++)
  {
    st = pmixmodel_->getParameterNameWithoutNamespace(it->first);
    if (it->second.substr(0, 5) != "delta")
      addParameter_(new Parameter("YNGKP_M3." + it->second, pmixmodel_->getParameterValue(st),
                              pmixmodel_->getParameter(st).hasConstraint() ? pmixmodel_->getParameter(st).getConstraint()->clone() : 0, true));
  }

  for (unsigned int i = 1; i < nbOmega; i++)
  {
    addParameter_(new Parameter("YNGKP_M3.delta" + TextTools::toString(i), 0.5, new IntervalConstraint(NumConstants::MILLI(), 999, true, true, NumConstants::MILLI()), true));
  }

  // look for synonymous codons
  for (synfrom_ = 1; synfrom_ < (int)gc->getSourceAlphabet()->getSize(); synfrom_++)
  {
    for (synto_ = 0; synto_ < synfrom_; synto_++)
    {
      if ((gc->areSynonymous(synfrom_, synto_))
          && (pmixmodel_->getNModel(0)->Qij(synfrom_, synto_) != 0)
          && (pmixmodel_->getNModel(1)->Qij(synfrom_, synto_) != 0))
        break;
    }
    if (synto_ < synfrom_)
      break;
  }

  if (synto_ == (int)gc->getSourceAlphabet()->getSize())
    throw Exception("Impossible to find synonymous codons");

  // update Matrices
  updateMatrices();
}

YNGKP_M3::YNGKP_M3(const YNGKP_M3& mod2) : AbstractBiblioMixedSubstitutionModel(mod2),
  pmixmodel_(new MixtureOfASubstitutionModel(*mod2.pmixmodel_)),
  synfrom_(mod2.synfrom_),
  synto_(mod2.synto_)
{}

YNGKP_M3& YNGKP_M3::operator=(const YNGKP_M3& mod2)
{
  AbstractBiblioMixedSubstitutionModel::operator=(mod2);

  pmixmodel_.reset(new MixtureOfASubstitutionModel(*mod2.pmixmodel_));
  synfrom_ = mod2.synfrom_;
  synto_ = mod2.synto_;

  return *this;
}

YNGKP_M3::~YNGKP_M3() {}

void YNGKP_M3::updateMatrices()
{
  for (unsigned int i = 0; i < lParPmodel_.size(); i++)
  {
    if (mapParNamesFromPmodel_.find(lParPmodel_[i].getName()) != mapParNamesFromPmodel_.end())
    {
      if (lParPmodel_[i].getName()[18] == 'V')
      {
        unsigned int ind = TextTools::toInt(lParPmodel_[i].getName().substr(19));
        double x = getParameterValue("omega0");
        for (unsigned j = 1; j < ind; j++)
        {
          x += getParameterValue("delta" + TextTools::toString(j));
        }
        lParPmodel_[i].setValue(x);
      }
      else
      {
        lParPmodel_[i].setValue(getParameter(getParameterNameWithoutNamespace(mapParNamesFromPmodel_[lParPmodel_[i].getName()])).getValue());
      }
    }
  }

  pmixmodel_->matchParametersValues(lParPmodel_);

  // homogeneization of the synonymous substitution rates


  Vdouble vd;

  for (unsigned int i = 0; i < pmixmodel_->getNumberOfModels(); i++)
  {
    vd.push_back(1 / pmixmodel_->getNModel(i)->Qij(synfrom_, synto_));
  }

  pmixmodel_->setVRates(vd);
}


