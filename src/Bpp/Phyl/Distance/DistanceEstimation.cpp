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

// bpp-phyl
#include <Bpp/Phyl/Likelihood/AutonomousSubstitutionProcess.h>

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
  Context context;

  size_t n = sites_->getNumberOfSequences();
  vector<string> names = sites_->getSequenceNames();
  dist_ = std::shared_ptr<DistanceMatrix>(new DistanceMatrix(names));
  optimizer_->setVerbose(static_cast<unsigned int>(max(static_cast<int>(verbose_) - 2, 0)));

  Newick reader;

  auto autoProc = dynamic_pointer_cast<AutonomousSubstitutionProcessInterface>(process_);
  auto procMb = dynamic_pointer_cast<SubstitutionProcessCollectionMember>(process_);
  if (!autoProc && !procMb)
    throw Exception("DistanceMatrix::computeMatrix : unknown process type. Ask developpers.");
  size_t maxTNb = 0;
  if (procMb){
    const auto& vTn = procMb->getCollection()->getTreeNumbers();
    maxTNb  = *std::max_element(vTn.begin(),vTn.end())+1;
  }

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

      auto phyloTree = make_shared<bpp::ParametrizablePhyloTree>(*reader.parenthesisToPhyloTree("(" + names[j] + ":0.01," + names[i] + ":0.01);", false, "", false, false));

      if (autoProc)
        autoProc->setPhyloTree(*phyloTree);
      else
        if (procMb)
        {
          auto& coll = procMb->collection();
          if (!coll.hasTreeNumber(maxTNb))
          {
            coll.addTree(phyloTree, maxTNb);
            procMb->setTreeNumber(maxTNb, false);
          }
          else
          {
            coll.replaceTree(phyloTree, maxTNb);
          }
        }

      auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, sites_, process_);

      auto llh = std::make_shared<SingleProcessPhyloLikelihood>(context, lik);

      // Optimization:
      optimizer_->setFunction(llh);
      optimizer_->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
      ParameterList params = llh->getBranchLengthParameters();
      params.addParameters(parameters_);
      optimizer_->init(params);
      optimizer_->optimize();

      // Store results:
      if (autoProc)
        (*dist_)(i, j) = (*dist_)(j, i) = llh->getParameterValue("BrLen0") + llh->getParameterValue("BrLen1");
      else
        (*dist_)(i, j) = (*dist_)(j, i) = llh->getParameterValue("BrLen0_"+TextTools::toString(maxTNb)) + llh->getParameterValue("BrLen1_"+TextTools::toString(maxTNb));

    }
    if (verbose_ > 1 && ApplicationTools::message)
      ApplicationTools::message->endLine();
  }
}

/******************************************************************************/
