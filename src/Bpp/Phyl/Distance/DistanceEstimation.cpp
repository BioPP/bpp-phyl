// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../SitePatterns.h"
#include "../Tree/Tree.h"
#include "../PatternTools.h"
#include "../Io/Newick.h"

#include "DistanceEstimation.h"

// From bpp-core:
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/AutoParameter.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/DistanceMatrix.h>

using namespace bpp;

// From the STL:
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

/******************************************************************************/

void DistanceEstimation::computeMatrix()
{
  size_t n = sites_->getNumberOfSequences();
  vector<string> names = sites_->getSequenceNames();
  dist_ = std::shared_ptr<DistanceMatrix>(new DistanceMatrix(names));
  optimizer_->setVerbose(static_cast<unsigned int>(max(static_cast<int>(verbose_) - 2, 0)));

  // const SiteContainer* sc = dynamic_cast<const SiteContainer*>(sites_);
  // const VectorProbabilisticSiteContainer* psc = dynamic_cast<const VectorProbabilisticSiteContainer*>(sites_);
  Newick reader;

  for (size_t i = 0; i < n; ++i)
  {
    (*dist_)(i, i) = 0;
    if (verbose_ == 1)
    {
      ApplicationTools::displayGauge(i, n - 1, '=');
    }
    for (size_t j = i + 1; j < n; ++j)
    {
      if (verbose_ > 1)
      {
        ApplicationTools::displayGauge(j - i - 1, n - i - 2, '=');
      }

      Context context;
  
      auto phyloTree = std::shared_ptr<bpp::PhyloTree>(reader.parenthesisToPhyloTree("(" + names[i] + ":0.01," + names[j] + ":0.01);", false, "", false, false));

      auto process = std::make_shared<RateAcrossSitesSubstitutionProcess>(model_, rateDist_, phyloTree);
      
      auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites_, process);

      auto llh = std::make_shared<SingleProcessPhyloLikelihood>(context, lik);

      // size_t d = sc ?
      //            SymbolListTools::getNumberOfDistinctPositions(sc->getSequence(i), sc->getSequence(j)) :
      //            SymbolListTools::getNumberOfDistinctPositions(*psc->getSequence(i), *psc->getSequence(j));
      // size_t g = sc ?
      //            SymbolListTools::getNumberOfPositionsWithoutGap(sc->getSequence(i), sc->getSequence(j)) :
      //            SymbolListTools::getNumberOfPositionsWithoutGap(*psc->getSequence(i), *psc->getSequence(j));

      // llh.setParameterValue("BrLen", g == 0 ? lik->getMinimumBranchLength() : std::max(lik->getMinimumBranchLength(), static_cast<double>(d) / static_cast<double>(g)));
      // Optimization:
      optimizer_->setFunction(llh);
      optimizer_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
      ParameterList params = llh->getBranchLengthParameters();
      params.addParameters(parameters_);
      optimizer_->init(params);
      optimizer_->optimize();
      // Store results:
      (*dist_)(i, j) = (*dist_)(j, i) = llh->getParameterValue("BrLen0") + llh->getParameterValue("BrLen1");

    }
    if (verbose_ > 1 && ApplicationTools::message)
      ApplicationTools::message->endLine();
  }
}

/******************************************************************************/
