//
// File: DistanceEstimation.cpp
// Authors:
//   Julien Dutheil
//   Vincent Ranwez
// Created: 2005-06-08 10:39:00
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
  vector<string> names = sites_->getSequencesNames();
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
    for (size_t j = i + 1; j < n; j++)
    {
      if (verbose_ > 1)
      {
        ApplicationTools::displayGauge(j - i - 1, n - i - 2, '=');
      }

      Context context;
  
      auto phyloTree = std::shared_ptr<bpp::PhyloTree>(reader.parenthesisToPhyloTree("("+names[i]+":0.01,"+names[j]+":0.01);", false, "", false, false));

      auto process=std::make_shared<RateAcrossSitesSubstitutionProcess>(model_, rateDist_, phyloTree);
      
      auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, *sites_, *process);

      SingleProcessPhyloLikelihood llh(context, lik);

      // size_t d = sc ?
      //            SymbolListTools::getNumberOfDistinctPositions(sc->getSequence(i), sc->getSequence(j)) :
      //            SymbolListTools::getNumberOfDistinctPositions(*psc->getSequence(i), *psc->getSequence(j));
      // size_t g = sc ?
      //            SymbolListTools::getNumberOfPositionsWithoutGap(sc->getSequence(i), sc->getSequence(j)) :
      //            SymbolListTools::getNumberOfPositionsWithoutGap(*psc->getSequence(i), *psc->getSequence(j));

      // llh.setParameterValue("BrLen", g == 0 ? lik->getMinimumBranchLength() : std::max(lik->getMinimumBranchLength(), static_cast<double>(d) / static_cast<double>(g)));
      // Optimization:
      optimizer_->setFunction(&llh);
      optimizer_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
      ParameterList params = llh.getBranchLengthParameters();
      params.addParameters(parameters_);
      optimizer_->init(params);
      optimizer_->optimize();
      // Store results:
      (*dist_)(i, j) = (*dist_)(j, i) = llh.getParameterValue("BrLen0") + llh.getParameterValue("BrLen1");

    }
    if (verbose_ > 1 && ApplicationTools::message)
      ApplicationTools::message->endLine();
  }
}

/******************************************************************************/
