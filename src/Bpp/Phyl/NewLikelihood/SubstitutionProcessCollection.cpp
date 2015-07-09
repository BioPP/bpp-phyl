//
// File: SubstitutionProcessCollection.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 13 juin 2013, à 11h 57
//

/*
  Copyright or <A9> or Copr. CNRS, (November 16, 2004)

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

#include "SubstitutionProcessCollection.h"

#include "SubstitutionProcessCollectionMember.h"

#include "../Tree/TreeTemplateTools.h"

#include <Bpp/Numeric/Prob/ConstantDistribution.h>


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
    mSubProcess_[it->first]=dynamic_cast<SubstitutionProcessCollectionMember*>(it->second->clone());
}

SubstitutionProcessCollection& SubstitutionProcessCollection::operator=(const SubstitutionProcessCollection& set)
{
  clear();

  AbstractParameterAliasable::operator=(set);
  modelColl_=set.modelColl_;
  mModelToSubPro_=set.mModelToSubPro_;
  freqColl_=set.freqColl_;
  mFreqToSubPro_=set.mFreqToSubPro_;
  distColl_=set.distColl_;
  mDistToSubPro_=set.mDistToSubPro_;
  mVConstDist_=set.mVConstDist_;
  treeColl_=set.treeColl_;
  mTreeToSubPro_=set.mTreeToSubPro_;

  map<size_t, SubstitutionProcessCollectionMember*>::const_iterator it;
  
  for (it = set.mSubProcess_.begin(); it != set.mSubProcess_.end(); it++)
    mSubProcess_[it->first]=dynamic_cast<SubstitutionProcessCollectionMember*>(it->second->clone());

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
    delete it->second;

  mSubProcess_.clear();
}

void SubstitutionProcessCollection::addParametrizable(Parametrizable* parametrizable, size_t parametrizableIndex, bool withParameters)
{
  ParameterList pl;
  if (dynamic_cast<SubstitutionModel*>(parametrizable)){
    modelColl_.addObject(dynamic_cast<SubstitutionModel*>(parametrizable), parametrizableIndex);
    pl=modelColl_.getParametersForObject(parametrizableIndex);
  }
  else
    if (dynamic_cast<FrequenciesSet*>(parametrizable)){
      freqColl_.addObject(dynamic_cast<FrequenciesSet*>(parametrizable), parametrizableIndex);
      pl=freqColl_.getParametersForObject(parametrizableIndex);
    }
    else
      if (dynamic_cast<DiscreteDistribution*>(parametrizable)){
        distColl_.addObject(dynamic_cast<DiscreteDistribution*>(parametrizable), parametrizableIndex);
        pl=distColl_.getParametersForObject(parametrizableIndex);
      }
      else
        if (dynamic_cast<ParametrizableTree*>(parametrizable)){
          treeColl_.addObject(dynamic_cast<ParametrizableTree*>(parametrizable), parametrizableIndex);
          pl=treeColl_.getParametersForObject(parametrizableIndex);
        }

        else
          throw Exception("Unknown parametrizable object in SubstitutionProcessCollection::addParametrizable.");

  if (withParameters)
    addParameters_(pl);
}

ParameterList SubstitutionProcessCollection::getNonDerivableParameters() const
{
  // patch, to be fixed properly later
  return getIndependentParameters();

  ParameterList pl=distColl_.getIndependentParameters();
  pl.addParameters(modelColl_.getIndependentParameters());
  pl.addParameters(freqColl_.getIndependentParameters());
  
  return pl;
}
  
void SubstitutionProcessCollection::fireParameterChanged(const ParameterList& parameters)
{
  AbstractParameterAliasable::fireParameterChanged(parameters);

  ParameterList gAP=getAliasedParameters(parameters);

  gAP.addParameters(parameters);

  modelColl_.clearChanged();
  modelColl_.matchParametersValues(gAP);

  
  const vector<size_t>& vM=modelColl_.hasChanged();
  for (size_t i=0; i<vM.size(); i++)
  {
    const vector<size_t>& vs=mModelToSubPro_[vM[i]];
    for (size_t j=0; j<vs.size(); j++){
      mSubProcess_[vs[j]]->changedModel(vM[i]);
      if (mSubProcess_[vs[j]]->isStationary())
        mSubProcess_[vs[j]]->changedRoot();
    }
  }
  
  freqColl_.clearChanged();
  freqColl_.matchParametersValues(gAP);

  vector<size_t> keys=freqColl_.keys();

  const vector<size_t> vMf=freqColl_.hasChanged();
  for (size_t i=0; i<vMf.size(); i++)
  {
    const vector<size_t>& vs=mFreqToSubPro_[vMf[i]];
    for (size_t j=0; j<vs.size(); j++)
      mSubProcess_[vs[j]]->changedRoot();
  }
  

  // map of the SubProcess to be fired
  
  map<size_t, bool> toFire;

  distColl_.clearChanged();
  distColl_.matchParametersValues(gAP);
  const vector<size_t>& vD=distColl_.hasChanged();
  
  for (size_t i=0; i<vD.size(); i++)
  {
    const vector<size_t>& vs=mDistToSubPro_[vD[i]];
    for (size_t j=0; j<vs.size(); j++)
      toFire[vs[j]]=true;

    if (mVConstDist_.find(vD[i])!=mVConstDist_.end()){
      const DiscreteDistribution& dd=getRateDistribution(vD[i]);
      vector<size_t>& vv=mVConstDist_[vD[i]];
      
      for (size_t j=0;j<vv.size();j++){
        gAP.addParameter(new Parameter("Constant.value_"+TextTools::toString(10000*(vD[i]+1)+vv[j]),dd.getCategory(j)));
        dynamic_cast<ConstantDistribution*>(distColl_[10000*(vD[i]+1)+vv[j]])->setParameterValue("value",dd.getCategory(j));
        const vector<size_t>&  vs2=mDistToSubPro_[10000*(vD[i]+1)+vv[j]];
        for (size_t k=0; k<vs2.size(); k++)
          toFire[vs2[k]]=true;
      }
    }
  }

  
  treeColl_.clearChanged();
  treeColl_.matchParametersValues(gAP);
  
  const vector<size_t>& vT=treeColl_.hasChanged();
  for (size_t i=0; i<vT.size(); i++)
  {
    const vector<size_t>& vs=mTreeToSubPro_[vT[i]];      
    for (size_t j=0; j<vs.size(); j++)
      toFire[vs[j]]=true;
  }


  // send the message to subprocesses

  map<size_t, bool>::const_iterator it;
  
  for (it=toFire.begin(); it != toFire.end(); it++)
    mSubProcess_[it->first]->fireParameterChanged(gAP);
  
}


void SubstitutionProcessCollection::setNamespace(const string& prefix)
{
  AbstractParameterAliasable::setNamespace(prefix);
  for (std::map<size_t, SubstitutionProcessCollectionMember*>::iterator it=mSubProcess_.begin(); it != mSubProcess_.end(); it++)
    it->second->setNamespace(prefix);
}

void SubstitutionProcessCollection::aliasParameters(const std::string& p1, const std::string& p2) throw (ParameterNotFoundException, Exception) 
{
  AbstractParameterAliasable::aliasParameters(p1, p2);
  for (std::map<size_t, SubstitutionProcessCollectionMember*>::iterator it=mSubProcess_.begin(); it != mSubProcess_.end(); it++)
    if (it->second->hasParameter(p2))
    {
      string p=p2;
      it->second->deleteParameter_(p);
    }
}

void SubstitutionProcessCollection::unaliasParameters(const std::string& p1, const std::string& p2) throw (ParameterNotFoundException, Exception) 
{
  AbstractParameterAliasable::unaliasParameters(p1, p2);
  for (std::map<size_t, SubstitutionProcessCollectionMember*>::iterator it=mSubProcess_.begin(); it != mSubProcess_.end(); it++)
    it->second->updateParameters();
}

void SubstitutionProcessCollection::aliasParameters(std::map<std::string, std::string>& unparsedParams, bool verbose) throw (ParameterNotFoundException, Exception) 
{
  AbstractParameterAliasable::aliasParameters(unparsedParams, verbose);
  for (std::map<std::string, std::string>::iterator itp=unparsedParams.begin(); itp!=unparsedParams.end(); itp++)
  {
    string p2=itp->second;
    
    for (std::map<size_t, SubstitutionProcessCollectionMember*>::iterator it=mSubProcess_.begin(); it != mSubProcess_.end(); it++)
      if (it->second->hasParameter(p2))
        it->second->deleteParameter_(p2);
  }
}


void SubstitutionProcessCollection::addSubstitutionProcess(size_t nProc, std::map<size_t, std::vector<int> > mModBr, size_t nTree, size_t nRate, size_t nFreq)
{
  if (mSubProcess_.find(nProc)!=mSubProcess_.end())
    throw BadIntegerException("Already assigned substitution process",(int)nProc);
  
  if (!freqColl_.hasObject(nFreq))
    throw BadIntegerException("Wrong Root Frequencies Set number",(int)nFreq);
  if (!treeColl_.hasObject(nTree))
    throw BadIntegerException("Wrong Tree number",(int)nTree);
  if (!distColl_.hasObject(nRate))
    throw BadIntegerException("Wrong Rate distribution number",(int)nRate);

  SubstitutionProcessCollectionMember* pSMS=new SubstitutionProcessCollectionMember(this, nProc, nTree, nRate);
  pSMS->setRootFrequencies(nFreq);

  std::map<size_t, std::vector<int> >::iterator it;
  for (it=mModBr.begin(); it!=mModBr.end(); it++){
    pSMS->addModel(it->first, it->second);
    mModelToSubPro_[it->first].push_back(nProc);
  }
  
  pSMS->isFullySetUp();

  mTreeToSubPro_[nTree].push_back(nProc);
  mDistToSubPro_[nRate].push_back(nProc);
  mFreqToSubPro_[nFreq].push_back(nProc);
  
  mSubProcess_[nProc]=pSMS;
}

void SubstitutionProcessCollection::addSubstitutionProcess(size_t nProc, std::map<size_t, std::vector<int> > mModBr, size_t nTree, size_t nRate)
{
  if (mSubProcess_.find(nProc)!=mSubProcess_.end())
    throw BadIntegerException("Already assigned substitution process",(int)nProc);

  if (!treeColl_.hasObject(nTree))
    throw BadIntegerException("Wrong Tree number",(int)nTree);
  if (!distColl_.hasObject(nRate))
    throw BadIntegerException("Wrong Rate distribution number",(int)nRate);

  SubstitutionProcessCollectionMember* pSMS=new SubstitutionProcessCollectionMember(this, nProc, nTree, nRate);

  std::map<size_t, std::vector<int> >::iterator it;
  for (it=mModBr.begin(); it!=mModBr.end(); it++){
    pSMS->addModel(it->first, it->second);
    mModelToSubPro_[it->first].push_back(nProc);
  }

  pSMS->isFullySetUp();
  
  mTreeToSubPro_[nTree].push_back(nProc);
  mDistToSubPro_[nRate].push_back(nProc);
  
  mSubProcess_[nProc]=pSMS;
}

bool SubstitutionProcessCollection::hasBranchLengthParameter(const std::string& name) const
{
  return treeColl_.hasParameter(name);
}

ParameterList SubstitutionProcessCollection::getSubstitutionProcessParameters() const
{
  ParameterList pl=modelColl_.getParameters();
  pl.addParameters(treeColl_.getParameters());
  pl.addParameters(freqColl_.getParameters());
  pl.addParameters(distColl_.getParameters());

  return pl;
}



