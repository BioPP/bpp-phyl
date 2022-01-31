//
// File: SubstitutionProcessCollection.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: jeudi 13 juin 2013, ÃÂ  11h 57
//

/*
  Copyright or <A9> or Copr. CNRS, (November 16, 2004)
  
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

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/VectorTools.h>

#include "ParametrizablePhyloTree.h"
#include "SubstitutionProcessCollection.h"
#include "SubstitutionProcessCollectionMember.h"

using namespace bpp;
using namespace std;

SubstitutionProcessCollection::SubstitutionProcessCollection(const SubstitutionProcessCollection& set) :
  AbstractParameterAliasable(set),
  modelColl_(set.modelColl_),
  mModelToSubPro_(set.mModelToSubPro_),
  freqColl_(set.freqColl_),
  mFreqToSubPro_(set.mFreqToSubPro_),
  distColl_(set.distColl_),
  mVConstDist_(set.mVConstDist_),
  mDistToSubPro_(set.mDistToSubPro_),
  treeColl_(set.treeColl_),
  mTreeToSubPro_(set.mTreeToSubPro_),
  mSubProcess_()
{
  map<size_t, SubstitutionProcessCollectionMember*>::const_iterator it;

  for (it = set.mSubProcess_.begin(); it != set.mSubProcess_.end(); it++)
  {
    mSubProcess_[it->first] = dynamic_cast<SubstitutionProcessCollectionMember*>(it->second->clone());
  }
}

SubstitutionProcessCollection& SubstitutionProcessCollection::operator=(const SubstitutionProcessCollection& set)
{
  clear();

  AbstractParameterAliasable::operator=(set);
  modelColl_ = set.modelColl_;
  mModelToSubPro_ = set.mModelToSubPro_;
  freqColl_ = set.freqColl_;
  mFreqToSubPro_ = set.mFreqToSubPro_;
  distColl_ = set.distColl_;
  mDistToSubPro_ = set.mDistToSubPro_;
  mVConstDist_ = set.mVConstDist_;
  treeColl_ = set.treeColl_;
  mTreeToSubPro_ = set.mTreeToSubPro_;

  map<size_t, SubstitutionProcessCollectionMember*>::const_iterator it;

  for (it = set.mSubProcess_.begin(); it != set.mSubProcess_.end(); it++)
  {
    mSubProcess_[it->first] = dynamic_cast<SubstitutionProcessCollectionMember*>(it->second->clone());
  }

  return *this;
}

void SubstitutionProcessCollection::clear()
{
  resetParameters_();

  modelColl_.clear();
  mModelToSubPro_.clear();
  freqColl_.clear();
  mFreqToSubPro_.clear();
  distColl_.clear();
  mDistToSubPro_.clear();
  mVConstDist_.clear();
  treeColl_.clear();
  mTreeToSubPro_.clear();

  map<size_t, SubstitutionProcessCollectionMember*>::const_iterator it;

  for (it = mSubProcess_.begin(); it != mSubProcess_.end(); it++)
  {
    delete it->second;
  }

  mSubProcess_.clear();
}

void SubstitutionProcessCollection::addParametrizable(std::shared_ptr<Parametrizable> parametrizable, size_t parametrizableIndex, bool withParameters)
{
  if (parametrizableIndex<1)
    throw BadIntegerException("SubstitutionProcessCollection::addParametrizable: parametrizableIndex should be at least 1.",(int)parametrizableIndex);
  
  ParameterList pl;
  if (std::dynamic_pointer_cast<BranchModel>(parametrizable))
  {
    modelColl_.addObject(std::dynamic_pointer_cast<BranchModel>(parametrizable), parametrizableIndex);
    pl = modelColl_.getParametersForObject(parametrizableIndex);
  }
  else if (std::dynamic_pointer_cast<FrequencySet>(parametrizable))
  {
    freqColl_.addObject(std::dynamic_pointer_cast<FrequencySet>(parametrizable), parametrizableIndex);
    pl = freqColl_.getParametersForObject(parametrizableIndex);
  }
  else if (std::dynamic_pointer_cast<DiscreteDistribution>(parametrizable))
  {
    distColl_.addObject(std::dynamic_pointer_cast<DiscreteDistribution>(parametrizable), parametrizableIndex);
    pl = distColl_.getParametersForObject(parametrizableIndex);
  }
  else if (std::dynamic_pointer_cast<ParametrizablePhyloTree>(parametrizable))
  {
    treeColl_.addObject(std::dynamic_pointer_cast<ParametrizablePhyloTree>(parametrizable), parametrizableIndex);
    pl = treeColl_.getParametersForObject(parametrizableIndex);
  }

  else
    throw Exception("Unknown parametrizable object in SubstitutionProcessCollection::addParametrizable.");

  if (withParameters)
    addParameters_(pl);
}


ParameterList SubstitutionProcessCollection::getNonDerivableParameters() const
{
  ParameterList pl = distColl_.getIndependentParameters();
  pl.addParameters(modelColl_.getIndependentParameters());
  pl.addParameters(freqColl_.getIndependentParameters());

  pl.includeParameters(getAliasedParameters(pl));
  pl.includeParameters(getFromParameters(pl));

  return pl;
}

void SubstitutionProcessCollection::fireParameterChanged(const ParameterList& parameters)
{
  ParameterList gAP = getAliasedParameters(parameters);

  gAP.addParameters(parameters);

  modelColl_.clearChanged();
  modelColl_.matchParametersValues(gAP);

  freqColl_.clearChanged();
  freqColl_.matchParametersValues(gAP);

  // map of the SubProcess to be fired

  map<size_t, bool> toFire;

  distColl_.clearChanged();
  distColl_.matchParametersValues(gAP);
  const vector<size_t>& vD = distColl_.hasChanged();

  for (size_t i = 0; i < vD.size(); i++)
  {
    const vector<size_t>& vs = mDistToSubPro_[vD[i]];
    for (size_t j = 0; j < vs.size(); j++)
    {
      toFire[vs[j]] = true;
    }

    if (mVConstDist_.find(vD[i]) != mVConstDist_.end())
    {
      auto dd = getRateDistribution(vD[i]);
      vector<size_t>& vv = mVConstDist_[vD[i]];

      for (size_t j = 0; j < vv.size(); j++)
      {
        gAP.addParameter(new Parameter("Constant.value_" + TextTools::toString(10000 * (vD[i] + 1) + vv[j]), dd->getCategory(j)));
        std::dynamic_pointer_cast<ConstantDistribution>(distColl_[10000 * (vD[i] + 1) + vv[j]])->setParameterValue("value", dd->getCategory(j));
        const vector<size_t>&  vs2 = mDistToSubPro_[10000 * (vD[i] + 1) + vv[j]];
        for (size_t k = 0; k < vs2.size(); k++)
        {
          toFire[vs2[k]] = true;
        }
      }
    }
  }


  treeColl_.clearChanged();
  treeColl_.matchParametersValues(gAP);

  const vector<size_t>& vT = treeColl_.hasChanged();

  for (size_t i = 0; i < vT.size(); i++)
  {
    const vector<size_t>& vs = mTreeToSubPro_[vT[i]];
    for (size_t j = 0; j < vs.size(); j++)
    {
      toFire[vs[j]] = true;
    }
  }


  // send the message to subprocesses

  map<size_t, bool>::const_iterator it;

  for (it = toFire.begin(); it != toFire.end(); it++)
  {
    mSubProcess_[it->first]->fireParameterChanged(gAP);
  }
}


void SubstitutionProcessCollection::setNamespace(const string& prefix)
{
  AbstractParameterAliasable::setNamespace(prefix);
  for (std::map<size_t, SubstitutionProcessCollectionMember*>::iterator it = mSubProcess_.begin(); it != mSubProcess_.end(); it++)
  {
    it->second->setNamespace(prefix);
  }
}

void SubstitutionProcessCollection::aliasParameters(const std::string& p1, const std::string& p2)
{
  AbstractParameterAliasable::aliasParameters(p1, p2);
  for (std::map<size_t, SubstitutionProcessCollectionMember*>::iterator it = mSubProcess_.begin(); it != mSubProcess_.end(); it++)
  {
    if (it->second->hasParameter(p2))
    {
      string p = p2;
      it->second->deleteParameter_(p);
    }
  }
}

void SubstitutionProcessCollection::unaliasParameters(const std::string& p1, const std::string& p2)
{
  AbstractParameterAliasable::unaliasParameters(p1, p2);
  for (std::map<size_t, SubstitutionProcessCollectionMember*>::iterator it = mSubProcess_.begin(); it != mSubProcess_.end(); it++)
  {
    it->second->updateParameters();
  }
}

void SubstitutionProcessCollection::aliasParameters(std::map<std::string, std::string>& unparsedParams, bool verbose)
{
  AbstractParameterAliasable::aliasParameters(unparsedParams, verbose);
  for (std::map<std::string, std::string>::iterator itp = unparsedParams.begin(); itp != unparsedParams.end(); itp++)
  {
    string p2 = itp->second;

    for (std::map<size_t, SubstitutionProcessCollectionMember*>::iterator it = mSubProcess_.begin(); it != mSubProcess_.end(); it++)
    {
      if (it->second->hasParameter(p2))
        it->second->deleteParameter_(p2);
    }
  }
}


void SubstitutionProcessCollection::addSubstitutionProcess(size_t nProc, std::map<size_t, std::vector<unsigned int> > mModBr, size_t nTree, size_t nRate, size_t nFreq)
{
  addSubstitutionProcess(nProc, mModBr, nTree, nRate);

  if (!freqColl_.hasObject(nFreq))
    throw BadIntegerException("Wrong Root Frequencies Set number", (int)nFreq);

  SubstitutionProcessCollectionMember& pSMS = dynamic_cast<SubstitutionProcessCollectionMember&>(getSubstitutionProcess(nProc));
  pSMS.setRootFrequencies(nFreq);

  mFreqToSubPro_[nFreq].push_back(nProc);
}

void SubstitutionProcessCollection::addSubstitutionProcess(size_t nProc, std::map<size_t, std::vector<unsigned int> > mModBr, size_t nTree, size_t nRate)
{
  if (mSubProcess_.find(nProc) != mSubProcess_.end())
    throw BadIntegerException("Already assigned substitution process", (int)nProc);

  if (!treeColl_.hasObject(nTree))
    throw BadIntegerException("Wrong Tree number", (int)nTree);
  if (!distColl_.hasObject(nRate))
    throw BadIntegerException("Wrong Rate distribution number", (int)nRate);

  SubstitutionProcessCollectionMember* pSMS = new SubstitutionProcessCollectionMember(this, nProc, nTree, nRate);

  std::map<size_t, std::vector<unsigned int> >::iterator it;
  for (it = mModBr.begin(); it != mModBr.end(); it++)
  {
    pSMS->addModel(it->first, it->second);
    mModelToSubPro_[it->first].push_back(nProc);
  }

  pSMS->isFullySetUp();

  mTreeToSubPro_[nTree].push_back(nProc);
  mDistToSubPro_[nRate].push_back(nProc);

  mSubProcess_[nProc] = pSMS;
}

ParameterList SubstitutionProcessCollection::getSubstitutionProcessParameters() const
{
  ParameterList pl = modelColl_.getParameters();
  pl.addParameters(treeColl_.getParameters());
  pl.addParameters(freqColl_.getParameters());
  pl.addParameters(distColl_.getParameters());

  return pl;
}


void SubstitutionProcessCollection::addOnePerBranchSubstitutionProcess(size_t nProc, size_t nMod, size_t nTree, size_t nRate, size_t nFreq, const vector<string>& sharedParameterNames)
{
  addOnePerBranchSubstitutionProcess(nProc, nMod, nTree, nRate, sharedParameterNames);

  if (!freqColl_.hasObject(nFreq))
    throw BadIntegerException("Wrong Root Frequencies Set number", (int)nFreq);

  SubstitutionProcessCollectionMember& pSMS = dynamic_cast<SubstitutionProcessCollectionMember&>(getSubstitutionProcess(nProc));
  pSMS.setRootFrequencies(nFreq);

  mFreqToSubPro_[nFreq].push_back(nProc);
}


void SubstitutionProcessCollection::addOnePerBranchSubstitutionProcess(size_t nProc, size_t nMod, size_t nTree, size_t nRate, const vector<string>& sharedParameterNames)
{
  if (mSubProcess_.find(nProc) != mSubProcess_.end())
    throw BadIntegerException("Already assigned substitution process", (int)nProc);

  if (!treeColl_.hasObject(nTree))
    throw BadIntegerException("Wrong Tree number", (int)nTree);
  if (!distColl_.hasObject(nRate))
    throw BadIntegerException("Wrong Rate distribution number", (int)nRate);

  // Build new models and assign to a map

  auto tree = getTree(nTree);
  const auto model = getModel(nMod);

  vector<uint> ids = tree->getAllEdgesIndexes();
  sort(ids.begin(), ids.end());
  vector<size_t> vModN = getModelNumbers();

  size_t maxMod = *max_element(vModN.begin(), vModN.end());

  std::map<size_t, std::vector<unsigned int> > mModBr;
  mModBr[nMod] = vector<uint>(1, ids[0]);

  for (auto it = ids.begin() + 1; it != ids.end(); it++)
  {
    size_t mNb = maxMod + *it;
    addModel(std::shared_ptr<BranchModel>(model->clone()), mNb);
    mModBr[mNb] = vector<uint>(1, *it);
  }

  addSubstitutionProcess(nProc, mModBr, nTree, nRate);

  /////////////////////////////////////
  // Look for aliases in parameters

  ParameterList sharedParameters = model->getIndependentParameters();

  vector<string> sharedParameterNames2; // vector of the complete names of global parameters

  // First get correct parameter names
  size_t i, j;

  for (i = 0; i < sharedParameterNames.size(); i++)
  {
    if (sharedParameterNames[i].find("*") != string::npos)
    {
      for (j = 0; j < sharedParameters.size(); j++)
      {
        StringTokenizer stj(sharedParameterNames[i], "*", true, false);
        size_t pos1, pos2;
        string parn = sharedParameters[j].getName();
        bool flag(true);
        string g = stj.nextToken();
        pos1 = parn.find(g);
        if (pos1 != 0)
          flag = false;
        pos1 += g.length();
        while (flag && stj.hasMoreToken())
        {
          g = stj.nextToken();
          pos2 = parn.find(g, pos1);
          if (pos2 == string::npos)
          {
            flag = false;
            break;
          }
          pos1 = pos2 + g.length();
        }
        if (flag &&
            ((g.length() == 0) || (pos1 == parn.length()) || (parn.rfind(g) == parn.length() - g.length())))
          sharedParameterNames2.push_back(parn);
      }
    }
    else if (!sharedParameters.hasParameter(sharedParameterNames[i]))
      throw Exception("SubstitutionProcessCollection::addOnePerBranchSubstitutionProcess. Parameter '" + sharedParameterNames[i] + "' is not valid.");
    else
      sharedParameterNames2.push_back(sharedParameterNames[i]);
  }

  // Now alias all shared parameters on all nodes:

  string pnum = TextTools::toString(mModBr.begin()->first);

  for (const auto& pname : sharedParameterNames2)
  {
    for (auto it = mModBr.begin(); it != mModBr.end(); it++)
    {
      if (it != mModBr.begin())
        aliasParameters(pname + "_" + pnum, pname + "_" + TextTools::toString(it->first));
    }
  }
}
