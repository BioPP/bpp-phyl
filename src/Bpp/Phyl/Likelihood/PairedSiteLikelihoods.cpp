//
// File: PairedSiteLikelihoods.cpp
// Created by: Nicolas Rochette
// Created on: January 6, 2011
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

// From the STL
#include <vector>
#include <string>
#include <numeric>
#include <cmath>

#include "PairedSiteLikelihoods.h"
#include "TreeLikelihood.h"

using namespace std;
using namespace bpp;

/***
 * Constructors:
 ***/

PairedSiteLikelihoods::PairedSiteLikelihoods() :
  logLikelihoods_(),
  modelNames_()
{}

PairedSiteLikelihoods::PairedSiteLikelihoods(
  const vector<vector<double> >& siteLogLikelihoods,
  const vector<string>& modelNames) throw (Exception) :
  logLikelihoods_(siteLogLikelihoods),
  modelNames_(modelNames)
{
  if (modelNames_.size() != getNumberOfModels())
  {
    if (modelNames_.size() == 0)
      modelNames_.assign(getNumberOfModels(), string());
    else
      throw Exception("PairedSiteLikelihoods: There should be as many model names as model site-loglikelihoods records.");
  }

  if (this->getNumberOfModels() > 0)
  {
    for (vector<vector<double> >::const_iterator siteLLiks = siteLogLikelihoods.begin();
         siteLLiks != siteLogLikelihoods.end();
         ++siteLLiks)
    {
      if (siteLLiks->size() != getNumberOfSites())
        throw Exception("PairedSiteLikelihoods: Models site-loglikelihoods records do not have the same number of elements.");
    }
  }
}

void PairedSiteLikelihoods::appendModel(
  const vector<double>& siteLogLikelihoods,
  const string& modelName
  ) throw (Exception)
{
  if (getNumberOfModels() > 0 && siteLogLikelihoods.size() != getNumberOfSites())
    throw Exception("PairedSiteLikelihoods::appendModel: Model site-loglikelihoods record does not have the correct number of elements");

  logLikelihoods_.push_back(siteLogLikelihoods);
  modelNames_.push_back(modelName);
}

void PairedSiteLikelihoods::appendModel(const bpp::TreeLikelihood& treeLikelihood) throw (Exception)
{
  const vector<double>& siteLogLikelihoods = treeLikelihood.getLogLikelihoodForEachSite();
  const string& modelName = treeLikelihood.getTree().getName();

  PairedSiteLikelihoods::appendModel(siteLogLikelihoods, modelName);
}

void PairedSiteLikelihoods::appendModels(const PairedSiteLikelihoods& psl)  throw (Exception)
{
  if (getNumberOfModels() > 0 && psl.getNumberOfModels() > 0 && psl.getNumberOfSites() != getNumberOfSites())
    throw Exception("PairedSiteLikelihoods::appendModels: The two PairedSiteLikelihood objects have different number of sites.");

  logLikelihoods_.insert(logLikelihoods_.end(),
                        psl.logLikelihoods_.begin(),
                        psl.logLikelihoods_.end()
                        );

  modelNames_.insert(modelNames_.end(),
                    psl.modelNames_.begin(),
                    psl.modelNames_.end()
                    );
}


pair<vector<string>, vector<double> > PairedSiteLikelihoods::computeExpectedLikelihoodWeights (int replicates) const
{
  // Initialize the weights
  vector<double> weights(getNumberOfModels(), 0);

  // Sum the model weights over replicates
  for (int r = 0; r < replicates; ++r)
  {
    // Draw the pseudoreplicate
    vector<int> siteCounts = bootstrap(getNumberOfSites());

    // Compute the loglikelihood of each model for this replicate
    vector<double> models_logliks (getNumberOfModels(), 0);
    for (size_t m = 0; m < getNumberOfModels(); ++m)
    {
      const vector<double>& modelSiteLLiks = logLikelihoods_.at(m);
      double Y = 0;
      for (size_t s = 0; s < getNumberOfSites(); ++s)
      {
        Y += modelSiteLLiks.at(s) * siteCounts.at(s);
      }
      models_logliks.at(m) = Y;
    }

    // Get the best loglikelihood
    double Ymax = *max_element(models_logliks.begin(), models_logliks.end());

    // Get the exponentials of the loglikelihood differences
    // and the sum of these values
    vector<double> exp_logliks_diffs (getNumberOfModels(), 0);
    for (size_t m = 0; m < getNumberOfModels(); ++m)
    {
      exp_logliks_diffs.at(m) = exp(models_logliks.at(m) - Ymax);
    }

    double sumELLD = accumulate(exp_logliks_diffs.begin(), exp_logliks_diffs.end(), 0.0);

    // Get the models weights for this replicate
    for (size_t m = 0; m < getNumberOfModels(); ++m)
    {
      double w = exp_logliks_diffs.at(m) / sumELLD;
      weights.at(m) += w;
    }
  } //for replicates

  // Divide all weights by the number of replicates.
  for (vector<double>::iterator w = weights.begin(); w != weights.end(); ++w)
  {
    *w /= replicates;
  }

  return make_pair(modelNames_, weights);
}

std::vector<int> PairedSiteLikelihoods::bootstrap(std::size_t length, double scaling)
{
  vector<int> v(length, 0);

  for (size_t i = 0; i < static_cast<size_t>(static_cast<double>(length) * scaling + 0.5); ++i)
  {
    ++v[RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(length)];
  }

  return v;
}

