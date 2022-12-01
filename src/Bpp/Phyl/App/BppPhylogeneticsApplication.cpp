//
// File: BppPhylogeneticsApplication.cpp
// Authors:
//   Laurent GuÃÂ©guen, Julien Dutheil
// Created: 2021-06-15 15:14:00
//

/*
  Copyright or ÃÂ© or Copr. Development Team, (November 17, 2021)
  
  This software is a computer program whose purpose is to provide basal and
  utilitary classes. This file belongs to the Bio++ Project.
  
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



// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>

#include "BppPhylogeneticsApplication.h"
#include "PhylogeneticsApplicationTools.h"

using namespace std;
using namespace bpp;

/******************************************************************************/

map<size_t, std::shared_ptr<PhyloTree> > BppPhylogeneticsApplication::getPhyloTreesMap(
  const map<size_t, std::shared_ptr<const AlignmentDataInterface> >& mSites,
  map<string, string>& unparsedParams,
  const std::string& prefix,
  const std::string& suffix,
  bool suffixIsOptional) const
{
  map<size_t, std::shared_ptr<PhyloTree> > mpTree = PhylogeneticsApplicationTools::getPhyloTrees(params_, mSites, unparsedParams, prefix, suffix, suffixIsOptional, verbose_, warn_);

  // Scaling of trees:
  double scale = ApplicationTools::getDoubleParameter("input.tree.scale", params_, 1, suffix, suffixIsOptional, warn_);

  if (scale != 1)
  {
    if (verbose_)
    {
      ApplicationTools::displayResult("Trees are scaled by", scale);
    }

    for (auto it : mpTree)
    {
      it.second->scaleTree(scale);
    }
  }

  return mpTree;
}

unique_ptr<SubstitutionProcessCollection> BppPhylogeneticsApplication::getCollection(
  std::shared_ptr<const Alphabet> alphabet,
  std::shared_ptr<const GeneticCode> gCode,
  const map<size_t, std::shared_ptr<const AlignmentDataInterface> >& mSites,
  map<string, string>& unparsedParams,
  const std::string& prefix,
  const std::string& suffix,
  bool suffixIsOptional) const
{
  auto mpTree = getPhyloTreesMap(mSites, unparsedParams, prefix, suffix, suffixIsOptional);
  auto SPC = getCollection(alphabet, gCode, mSites, mpTree, unparsedParams, prefix, suffix, suffixIsOptional);
  return SPC;
}

std::unique_ptr<SubstitutionProcessCollection> BppPhylogeneticsApplication::getCollection(
  std::shared_ptr<const Alphabet> alphabet,
  std::shared_ptr<const GeneticCode> gCode,
  const map<size_t, std::shared_ptr<const AlignmentDataInterface> >& mSites,
  const map<size_t, std::shared_ptr<PhyloTree> >& mpTree,
  map<string, string>& unparsedParams,
  const std::string& prefix,
  const std::string& suffix,
  bool suffixIsOptional) const
{
  auto mDist = PhylogeneticsApplicationTools::getRateDistributions(params_, suffix, suffixIsOptional, verbose_);
  auto mModU = PhylogeneticsApplicationTools::getBranchModels(alphabet, gCode, mSites, params_, unparsedParams, suffix, suffixIsOptional, verbose_, warn_);
  auto mMod = PhylogeneticsApplicationTools::uniqueToSharedMap<BranchModelInterface>(mModU);
  auto mRootFreqU = PhylogeneticsApplicationTools::getRootFrequencySets(alphabet, gCode, mSites, params_, unparsedParams, suffix, suffixIsOptional, verbose_, warn_);
  auto mRootFreq = PhylogeneticsApplicationTools::uniqueToSharedMap<FrequencySetInterface>(mRootFreqU);
  auto mModelPathU = PhylogeneticsApplicationTools::getModelPaths(params_, mMod, suffix, suffixIsOptional, verbose_, warn_);
  auto mModelPath = PhylogeneticsApplicationTools::uniqueToSharedMap<ModelPath>(mModelPathU);
  auto mScenarioU = PhylogeneticsApplicationTools::getModelScenarios(params_, mModelPath, mMod, suffix, suffixIsOptional, verbose_, warn_);
  auto mScenario = PhylogeneticsApplicationTools::uniqueToSharedMap<ModelScenario>(mScenarioU);

  auto SPC = PhylogeneticsApplicationTools::getSubstitutionProcessCollection(alphabet, gCode, mpTree, mMod, mRootFreq, mDist, mScenario, params_, unparsedParams, suffix, suffixIsOptional, verbose_, warn_);

  return SPC;
}


map<size_t, std::unique_ptr<SequenceEvolution> > BppPhylogeneticsApplication::getProcesses(
  SubstitutionProcessCollection& collection,
  map<string, string>& unparsedParams,
  const std::string& suffix,
  bool suffixIsOptional) const
{
  return PhylogeneticsApplicationTools::getSequenceEvolutions(
    collection, params_, unparsedParams, suffix, suffixIsOptional, verbose_, warn_);
}


std::unique_ptr<PhyloLikelihoodContainer> BppPhylogeneticsApplication::getPhyloLikelihoods(
  Context& context,
  map<size_t, shared_ptr<SequenceEvolution> > mSeqEvol,
  SubstitutionProcessCollection& collection,
  const map<size_t, shared_ptr<const AlignmentDataInterface> >& mSites,
  const std::string& suffix,
  bool suffixIsOptional) const
{
  return PhylogeneticsApplicationTools::getPhyloLikelihoodContainer(
    context, collection, mSeqEvol, mSites, params_, suffix, suffixIsOptional, verbose_, warn_);
}


void BppPhylogeneticsApplication::fixLikelihood(
  shared_ptr<const Alphabet> alphabet,
  shared_ptr<const GeneticCode> gCode,
  shared_ptr<PhyloLikelihoodInterface> phylolik,
  const std::string& suffix,
  bool suffixIsOptional) const
{
  double logL = phylolik->getValue();

  if (!std::isnormal(logL))
  {
    // This may be due to null branch lengths, leading to null likelihood!
    ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
    ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
    ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
    ParameterList pl = phylolik->getBranchLengthParameters();

    for (size_t i = 0; i < pl.size(); i++)
    {
      if (pl[i].getValue() < 0.000001)
        pl[i].setValue(0.001);
    }
    phylolik->matchParametersValues(pl);
    logL = phylolik->getValue();
  }

  ApplicationTools::displayMessage("");
  ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
  if (!std::isnormal(logL))
  {
    ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");

    map<size_t, shared_ptr<SingleDataPhyloLikelihoodInterface> > mSD;

    if (dynamic_pointer_cast<SingleDataPhyloLikelihoodInterface>(phylolik))
      mSD[1] = dynamic_pointer_cast<SingleDataPhyloLikelihoodInterface>(phylolik);
    else
    {
      auto sOAP = dynamic_pointer_cast<SetOfPhyloLikelihood>(phylolik);
      if (sOAP)
      {
        const vector<size_t>& nSD = sOAP->getNumbersOfPhyloLikelihoods();

        for (size_t iSD = 0; iSD < nSD.size(); ++iSD)
        {
          auto pASDP = dynamic_pointer_cast<SingleDataPhyloLikelihoodInterface>(sOAP->getPhyloLikelihood(nSD[iSD]));

          if (pASDP)
            mSD[nSD[iSD]] = pASDP;
        }
      }
    }

    for (auto& itm : mSD)
    {
      ApplicationTools::displayWarning("Checking for phyloLikelihood " + TextTools::toString(itm.first));

      if (!std::isnormal(itm.second->getValue()))
      {
        auto sDP = itm.second;

        auto vData = std::shared_ptr<AlignmentDataInterface>(sDP->getData()->clone());

        auto vSC = std::dynamic_pointer_cast<SiteContainerInterface>(vData);
        auto pSC = std::dynamic_pointer_cast<ProbabilisticSiteContainerInterface>(vData);

        if (AlphabetTools::isCodonAlphabet(alphabet.get()))
        {
          bool f = false;
          size_t s;
          for (size_t i = 0; i < vData->getNumberOfSites(); ++i)
          {
            if (!std::isnormal(sDP->getLogLikelihoodForASite(i)))
            {
              if (vSC)
              {
                const Site& site = vSC->getSite(i);
                s = site.size();
                for (size_t j = 0; j < s; ++j)
                {
                  if (gCode->isStop(site.getValue(j)))
                  {
                    (*ApplicationTools::error << "Stop Codon at site " << site.getCoordinate() << " in sequence " << vData->getSequence(j).getName()).endLine();
                    f = true;
                  }
                }
              }
              else
              {
                const ProbabilisticSite& site = pSC->getSite(i);
                s = site.size();
                for (size_t j = 0; j < s; ++j)
                {
                  bool g = false;
                  for (int st = 0; !g && st < static_cast<int>(alphabet->getSize()); ++st)
                  {
                    g = (site.getStateValueAt(j, st) != 0 && !gCode->isStop(st));
                  }

                  if (!g)
                  {
                    (*ApplicationTools::error << "Only stop Codons at site " << site.getCoordinate() << " in sequence " << vData->getSequence(j).getName()).endLine();
                    f = true;
                  }
                }
              }
            }
          }
          if (f)
            exit(-1);
        }

        // Then remove saturated positions

        bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites", params_, false, suffix, suffixIsOptional, warn_);

        if (removeSaturated)
        {
          ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);
          for (size_t i = vData->getNumberOfSites(); i > 0; --i)
          {
            if (!std::isnormal(sDP->getLogLikelihoodForASite(i - 1)))
            {
              ApplicationTools::displayResult("Ignore saturated site", vData->getSite(i - 1).getCoordinate());
              vData->deleteSites(i - 1, 1);
            }
          }
          ApplicationTools::displayResult("Number of sites retained", vData->getNumberOfSites());

          sDP->setData(vData);
        }

        logL = sDP->getValue();

        if (!std::isnormal(logL))
        {
          ApplicationTools::displayError("!!! Looking at each site:");
          for (unsigned int i = 0; i < vData->getNumberOfSites(); i++)
          {
            (*ApplicationTools::error << "Site " << vData->getSite(i).getCoordinate() << "\tlog likelihood = " << sDP->getLogLikelihoodForASite(i)).endLine();
          }
          ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
          exit(1);
        }
        else
          ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
      }
    }
  }
}


void BppPhylogeneticsApplication::displayParameters(const PhyloLikelihoodInterface& tl, bool displaylL) const
{
  // Write parameters to screen:
  if (displaylL)
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl.getValue(), 15));

  if (tl.getNumberOfParameters() - tl.getBranchLengthParameters().size() >= 30)
    ApplicationTools::displayMessage("Too many parameters for screen output!");
  else
  {
    ParameterList parameters = tl.getParameters();
    parameters.deleteParameters(tl.getBranchLengthParameters().getParameterNames(), false);
    for (size_t i = 0; i < parameters.size(); i++)
    {
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
  }
}
