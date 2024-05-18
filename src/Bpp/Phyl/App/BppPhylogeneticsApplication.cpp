// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>

#include "BppPhylogeneticsApplication.h"
#include "PhylogeneticsApplicationTools.h"
#include <Bpp/Seq/Container/SiteContainerTools.h>

using namespace std;
using namespace bpp;

/******************************************************************************/

map<size_t, std::shared_ptr<PhyloTree>> BppPhylogeneticsApplication::getPhyloTreesMap(
    const map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mSites,
    map<string, string>& unparsedParams,
    const std::string& prefix,
    const std::string& suffix,
    bool suffixIsOptional) const
{
  map<size_t, std::shared_ptr<PhyloTree>> mpTree = PhylogeneticsApplicationTools::getPhyloTrees(params_, mSites, unparsedParams, prefix, suffix, suffixIsOptional, verbose_, warn_);

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
    const map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mSites,
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
    const map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mSites,
    const map<size_t, std::shared_ptr<PhyloTree>>& mpTree,
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


map<size_t, std::unique_ptr<SequenceEvolution>> BppPhylogeneticsApplication::getProcesses(
    shared_ptr<SubstitutionProcessCollection> collection,
    map<string, string>& unparsedParams,
    const std::string& suffix,
    bool suffixIsOptional) const
{
  return PhylogeneticsApplicationTools::getSequenceEvolutions(
        collection, params_, unparsedParams, suffix, suffixIsOptional, verbose_, warn_);
}


std::shared_ptr<PhyloLikelihoodContainer> BppPhylogeneticsApplication::getPhyloLikelihoods(
    Context& context,
    map<size_t, shared_ptr<SequenceEvolution>> mSeqEvol,
    shared_ptr<SubstitutionProcessCollection> collection,
    const map<size_t, shared_ptr<const AlignmentDataInterface>>& mSites,
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
    ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.001.");
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

    map<size_t, shared_ptr<SingleDataPhyloLikelihoodInterface>> mSD;

    if (dynamic_pointer_cast<SingleDataPhyloLikelihoodInterface>(phylolik))
      mSD[1] = dynamic_pointer_cast<SingleDataPhyloLikelihoodInterface>(phylolik);
    else
    {
      auto sOAP = dynamic_pointer_cast<PhyloLikelihoodSetInterface>(phylolik);
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

        if (AlphabetTools::isCodonAlphabet(*alphabet))
        {
          bool f = false;
          size_t s;
          for (size_t i = 0; i < vData->getNumberOfSites(); ++i)
          {
            if (!std::isnormal(sDP->getLogLikelihoodForASite(i)))
            {
              if (vSC)
              {
                const Site& site = vSC->site(i);
                s = site.size();
                for (size_t j = 0; j < s; ++j)
                {
                  if (gCode->isStop(site.getValue(j)))
                  {
                    (*ApplicationTools::error << "Stop Codon at site " << site.getCoordinate() << " in sequence " << vData->sequence(j).getName()).endLine();
                    f = true;
                  }
                }
              }
              else
              {
                const ProbabilisticSite& site = pSC->site(i);
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
                    (*ApplicationTools::error << "Only stop Codons at site " << site.getCoordinate() << " in sequence " << vData->sequence(j).getName()).endLine();
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
              ApplicationTools::displayResult("Ignore saturated site", vData->site(i - 1).getCoordinate());
              vData->deleteSites(i - 1, 1);
            }
          }
          ApplicationTools::displayResult("Number of sites retained", vData->getNumberOfSites());

          sDP->setData(vData);
        }

        logL = sDP->getValue();

        vector<size_t> vsiteok, vsitemin;

        if (!std::isnormal(logL))
        {
          ApplicationTools::displayError("!!! Removing problem sites:");
          for (unsigned int i = 0; i < vData->getNumberOfSites(); i++)
          {
            auto x = sDP->getLogLikelihoodForASite(i);
            if (!std::isnormal(x))
            {
              (*ApplicationTools::error << "Site " << vData->site(i).getCoordinate() << "\tlog likelihood = " << x).endLine();
              vsitemin.push_back(i);
            }
            else
              vsiteok.push_back(i);
          }

          shared_ptr<AlignmentDataInterface>  vDataok = SiteContainerTools::getSelectedSites(*vData, vsiteok);
//          auto vDatamin = SiteContainerTools::getSelectedSites(*vData, vsitemin); Not taken into account yet

          sDP->setData(vDataok);
          logL = sDP->getValue();
          ApplicationTools::displayResult("Filtered log likelihood", TextTools::toString(-logL, 15));

          // auto phylo2 = itm.second->clone();  To be finished
          // phylo2->setData(*vDatamin);
          // auto logL2 = phylo2->getValue();

          // ApplicationTools::displayResult("Left log likelihood", TextTools::toString(-logL2, 15));
        }
        else
          ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
      }
      else
        ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-itm.second->getValue(), 15));
    }
  }
}


void BppPhylogeneticsApplication::displayParameters(const PhyloLikelihoodInterface& tl, bool displaylL) const{

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
