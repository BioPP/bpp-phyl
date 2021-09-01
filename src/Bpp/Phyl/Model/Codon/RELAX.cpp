//
// File: RELAX.cpp
// Created by:  Keren Halabi
// Created on: August 2018
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

#include "RELAX.h"
#include "YN98.h"
#include "../MixtureOfASubstitutionModel.h"
#include <math.h>                     /* pow */

#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/SimpleDiscreteDistribution.h>

using namespace bpp;

using namespace std;

/******************************************************************************/

RELAX::RELAX(const GeneticCode* gc, std::shared_ptr<FrequencySet> codonFreqs) :
  YNGP_M("RELAX.") // RELAX currenly inherits from YNGP_M as well, since it uses kappa and instead of the 5 GTR parameters
{
  // set the initial omegas distribution
  vector<double> omega_initials, omega_frequencies_initials;
  omega_initials.push_back(0.5); omega_initials.push_back(1); omega_initials.push_back(2);
  omega_frequencies_initials.push_back(0.333333); omega_frequencies_initials.push_back(0.333333); omega_frequencies_initials.push_back(0.333334);

  SimpleDiscreteDistribution* psdd = new SimpleDiscreteDistribution(omega_initials, omega_frequencies_initials);

  map<string, DiscreteDistribution*> mpdd;
  mpdd["omega"] = psdd;

  // build the submodel as a basic Yang Nielsen model (with kappa instead of 5 GTR nucleotide substituion rate parameters) 
  unique_ptr<YN98> yn98(new YN98(gc, codonFreqs));

  // initialize the site model with the initial omegas distribution
  pmixmodel_.reset(new MixtureOfASubstitutionModel(gc->getSourceAlphabet(), yn98.get(), mpdd));
  pmixsubmodel_=dynamic_cast<const MixtureOfASubstitutionModel*>(&getMixedModel());      

  delete psdd; // delete the initial omegas distibution, that is already embedded in the mixture model
  
  vector<int> supportedChars = yn98->getAlphabetStates();

  // mapping the parameters
  ParameterList pl = pmixmodel_->getParameters();
  for (size_t i = 0; i < pl.size(); i++)
  {
    lParPmodel_.addParameter(Parameter(pl[i])); // add the parameter to the biblio wrapper instance - see Laurent's response in https://groups.google.com/forum/#!searchin/biopp-help-forum/likelihood$20ratio$20test|sort:date/biopp-help-forum/lH8MYit_Mr8/2CBND79B11YJ
  }

  // v consists of 9 shared theta parameters, that are used for the F3X4 estimation of codon frequencies  
  vector<std::string> v = dynamic_cast<YN98*>(pmixmodel_->getNModel(0))->getFrequencySet()->getParameters().getParameterNames();

  for (size_t i = 0; i < v.size(); i++)
  {
    mapParNamesFromPmodel_[v[i]] = v[i].substr(5);
  }

  // map the parameters of RELAX to the parameters of the sub-models
  mapParNamesFromPmodel_["YN98.kappa"] = "kappa";
  mapParNamesFromPmodel_["YN98.omega_Simple.V1"] = "p";            // omega0=p*omega1 (p is re-parameterization of omega0)
  mapParNamesFromPmodel_["YN98.omega_Simple.V2"] = "omega1";     
  mapParNamesFromPmodel_["YN98.omega_Simple.theta1"] = "theta1";  // frequency of omega1 (denoted in YNGP_M2's documentation as p0)
  mapParNamesFromPmodel_["YN98.omega_Simple.V3"] = "omega2";
  mapParNamesFromPmodel_["YN98.omega_Simple.theta2"] = "theta2";  // theta2 = (p1/(p1+p2))
  /* codon frequencies parameterization using F3X4: for each _Full.theta, corresponding to a a codon position over {0,1,2}:
  getFreq_(0) = theta1 * (1. - theta);
  getFreq_(1) = (1 - theta2) * theta;
  getFreq_(2) = theta2 * theta;
  getFreq_(3) = (1 - theta1) * (1. - theta); */

  string st;
  for (map<string, string>::iterator it = mapParNamesFromPmodel_.begin(); it != mapParNamesFromPmodel_.end(); it++)
  {
    st = pmixmodel_->getParameterNameWithoutNamespace(it->first);
    if (it->second.substr(0, 5) != "omega" && it->second.substr(0, 5) !="p")
    {
      addParameter_(new Parameter("RELAX." + it->second, pmixmodel_->getParameterValue(st),
                              pmixmodel_->getParameter(st).hasConstraint() ? std::shared_ptr<Constraint>(pmixmodel_->getParameter(st).getConstraint()->clone()) : 0));
	}
  }

  /* set the below parameters that are used for parameterizing the omega parameters of the sumodels of type YN98 as autoparameters to supress exceptions when constraints of the YN98 omega parameters are exceeded
  YN98_0.omega = (RELAX.p * RELAX.omega1) ^ RELAX.k
  YN98_1.omega = RELAX.omega1 ^ RELAX.k  
  YN98_2.omega = RELAX.omega2 ^ RELAX.k */
  // reparameterization of omega0: RELAX.omega0 = RELAX.p*RELAX.omega1
  addParameter_(new Parameter("RELAX.p", 0.5, std::make_shared<IntervalConstraint>(0.01, 1, true, true)));
  addParameter_(new Parameter("RELAX.omega1", 1, std::make_shared<IntervalConstraint>(0.1, 1, true, true)));

  // the upper bound of omega3 in its submodel is 999, so I must restrict upperBound(RELAX.omega2)^upperBound(RELAX.k)<=999 -> set maximal omega to 5  
  addParameter_(new Parameter("RELAX.omega2", 2, std::make_shared<IntervalConstraint>(1, 999, true, true)));

  // add a selection intensity parameter k, which is 1 in the null case
  addParameter_(new Parameter("RELAX.k", 1, std::make_shared<IntervalConstraint>(0, 10, false, true))); // selection intensity parameter for purifying and neutral selection parameters 

  // look for synonymous codons
  // assumes that the states number follow the map in the genetic code and thus:
  // synfrom_ = index of source codon
  // synto_ = index of destination codon
  for (synfrom_ = 1; synfrom_ < supportedChars.size(); ++synfrom_)
  {
    for (synto_ = 0; synto_ < synfrom_; ++synto_)
    {
      if (gc->areSynonymous(supportedChars[synfrom_], supportedChars[synto_])
          && (pmixsubmodel_->getSubNModel(0)->Qij(synfrom_, synto_) != 0)
          && (pmixsubmodel_->getSubNModel(1)->Qij(synfrom_, synto_) != 0))
        break;
    }
    if (synto_ < synfrom_)
      break;
  }

  if (synto_ == supportedChars.size())
    throw Exception("Impossible to find synonymous codons");

  // update the 3 rate matrices of the model (strict BG or strict FG)
  computeFrequencies(false);
  updateMatrices();
}



void RELAX::updateMatrices()
{
  // update the values of the sub-model parameters, that are used in the 3 rate matrices
  for (unsigned int i = 0; i < lParPmodel_.size(); i++)
  {
    // first update the values of the non omega patrameters
    const string& np = lParPmodel_[i].getName();
    if (np.size()>19 && np[18] == 'V')
    {
      double k = getParameterValue("k"); // get the value of k
      int ind = -1 ;
      if (np.size()>19) {
        ind = TextTools::to<int>(np.substr(19)) - 1; // ind corresponds to the index of the omega that belongs to submodel ind+1
      }
      // change ind not to be -1 to allow names omega0, omega1, omega2
      double omega;
      if (ind == 0)
      {                         // handle omega0 differently due to reparameterization via RELAX.p
        omega = getParameterValue("p") * getParameterValue("omega1");
        omega = pow(omega,k);
        if (omega < 0.002)
        {
          omega = 0.002;
        }
      }
      else if (ind == 1)
      {
        omega = getParameterValue("omega1");
        omega = pow (omega, k);  
        if (omega <= 0.002)
        {
          omega = 0.002;
        }
      }
      else
      {
        omega = getParameterValue("omega2");
        omega = pow (omega, k);  
        if (omega > 999)
        {
          omega = 999;
        }
      }
      lParPmodel_[i].setValue(omega);
    }
    else
    {
      lParPmodel_[i].setValue(getParameter(getParameterNameWithoutNamespace(mapParNamesFromPmodel_[np])).getValue());
    }
  }


  pmixmodel_->matchParametersValues(lParPmodel_);

  // normalize the synonymous substitution rate in all the Q matrices of the 3 submodels to be the same
  Vdouble vd;

  for (unsigned int i = 0; i < pmixmodel_->getNumberOfModels(); i++)
  {
    vd.push_back(1 / pmixsubmodel_->getSubNModel(i)->Qij(synfrom_, synto_));
  }

  pmixmodel_->setVRates(vd);
}
