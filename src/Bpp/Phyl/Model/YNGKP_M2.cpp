//
// File: YNGKP_M2.cpp
// Created by:  Laurent Gueguen
// Created on: May 2010
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

#include "YNGKP_M2.h"
#include "YN98.h"
#include "FrequenciesSet.h"

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

using namespace bpp;

using namespace std;

/******************************************************************************/

YNGKP_M2::YNGKP_M2(const GeneticCode* gc, FrequenciesSet* codonFreqs) :
  MixedSubstitutionModel(gc->getSourceAlphabet(), "YNGKP_M2."), pmixmodel_(0),
  mapParNamesFromPmodel_(), lParPmodel_()
{
  // build the submodel

  vector<double> v1, v2;
  v1.push_back(0.5); v1.push_back(1); v1.push_back(2);
  v2.push_back(0.333333); v2.push_back(0.333333); v2.push_back(0.333334);

  SimpleDiscreteDistribution* psdd=new SimpleDiscreteDistribution(v1,v2);

  map<string, DiscreteDistribution*> mpdd;
  mpdd["omega"]=psdd;

  pmixmodel_= new MixtureOfASubstitutionModel(gc->getSourceAlphabet(),
                                              new YN98(gc, codonFreqs),
                                              mpdd);
  delete psdd;

  // mapping the parameters
  
  ParameterList pl=pmixmodel_->getParameters();
  for (unsigned int i=0;i<pl.size();i++)
    lParPmodel_.addParameter(Parameter(pl[i]));

  vector<std::string> v=dynamic_cast<YN98*>(pmixmodel_->getNModel(0))->getFreq().getParameters().getParameterNames();

  for (unsigned int i=0;i<v.size();i++)
    mapParNamesFromPmodel_[v[i]]=getParameterNameWithoutNamespace("YNGKP_M2."+v[i].substr(5));

  mapParNamesFromPmodel_["YN98.kappa"]="kappa";
  mapParNamesFromPmodel_["YN98.omega_Simple.V1"]="omega0";
  mapParNamesFromPmodel_["YN98.omega_Simple.theta1"]="theta1";
  mapParNamesFromPmodel_["YN98.omega_Simple.V3"]="omega2";
  mapParNamesFromPmodel_["YN98.omega_Simple.theta2"]="theta2";

  // specific parameters
  
  string st;
  for (map<string,string>::iterator it=mapParNamesFromPmodel_.begin(); it!= mapParNamesFromPmodel_.end(); it++){
    st=pmixmodel_->getParameterNameWithoutNamespace(it->first);
    if (it->second.substr(0,5)!="omega")
      addParameter_(Parameter("YNGKP_M2."+it->second, pmixmodel_->getParameterValue(st),
                              pmixmodel_->getParameter(st).hasConstraint()? pmixmodel_->getParameter(st).getConstraint()->clone():0,true));
  }
  
  addParameter_(Parameter("YNGKP_M2.omega0", 0.5, &Parameter::PROP_CONSTRAINT_EX));
  addParameter_(Parameter("YNGKP_M2.omega2", 2, new ExcludingInterval(1,999), true));

  //update Matrices
  
  updateMatrices();
}

YNGKP_M2::YNGKP_M2(const YNGKP_M2& mod2) : MixedSubstitutionModel(mod2),
                                           pmixmodel_(new MixtureOfASubstitutionModel(*mod2.pmixmodel_)),
                                           mapParNamesFromPmodel_(mod2.mapParNamesFromPmodel_),
                                           lParPmodel_(mod2.lParPmodel_)
{
  
}

YNGKP_M2& YNGKP_M2::operator=(const YNGKP_M2& mod2)
{
  MixedSubstitutionModel::operator=(mod2);

  pmixmodel_=new MixtureOfASubstitutionModel(*mod2.pmixmodel_);
  mapParNamesFromPmodel_=mod2.mapParNamesFromPmodel_;
  lParPmodel_=mod2.lParPmodel_;
  
  return *this;
}

YNGKP_M2::~YNGKP_M2()
{
  if (pmixmodel_)
    delete pmixmodel_;
}

void YNGKP_M2::updateMatrices()
{
  for (unsigned int i=0;i<lParPmodel_.size();i++)
    if (hasParameter(mapParNamesFromPmodel_[lParPmodel_[i].getName()]))
      lParPmodel_[i].setValue(getParameter(mapParNamesFromPmodel_[lParPmodel_[i].getName()]).getValue());
  
  pmixmodel_->matchParametersValues(lParPmodel_);
}
