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

#include "../TreeTemplateTools.h"

//#include "AbstractSubstitutionModel.h"

using namespace bpp;
using namespace std;

SubstitutionProcessCollection::SubstitutionProcessCollection(const SubstitutionProcessCollection& set) :
  AbstractParameterAliasable(set),
  modelColl_(set.modelColl_),
  mModelToSubPro_(set.mModelToSubPro_),
  freqColl_(set.freqColl_),
  mFreqToSubPro_(set.mFreqToSubPro_),
  distColl_(set.distColl_),
  mDistToSubPro_(set.mDistToSubPro_),
  treeColl_(set.treeColl_),
  mTreeToSubPro_(set.mTreeToSubPro_),
  vSubProcess_()
{
  for (size_t i=0;i<set.vSubProcess_.size(); i++)
    vSubProcess_.push_back(dynamic_cast<SubstitutionProcessCollectionMember*>(set.vSubProcess_[i]->clone()));
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
  treeColl_=set.treeColl_;
  mTreeToSubPro_=set.mTreeToSubPro_;
  for (size_t i=0;i<set.vSubProcess_.size(); i++)
    vSubProcess_.push_back(dynamic_cast<SubstitutionProcessCollectionMember*>(set.vSubProcess_[i]->clone()));

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
  treeColl_.clear();
  mTreeToSubPro_.clear();
  for (size_t i=0;i<vSubProcess_.size(); i++)
    delete vSubProcess_[i];

  vSubProcess_.clear();
}

void SubstitutionProcessCollection::addParametrizable(Parametrizable* parametrizable, size_t parametrizableIndex)
{
  ParameterList pl;
  if (dynamic_cast<SubstitutionModel*>(parametrizable)){
    modelColl_.addObject(dynamic_cast<SubstitutionModel*>(parametrizable), parametrizableIndex);
    pl=modelColl_[parametrizableIndex]->getIndependentParameters();
  }
  else
    if (dynamic_cast<FrequenciesSet*>(parametrizable)){
      freqColl_.addObject(dynamic_cast<FrequenciesSet*>(parametrizable), parametrizableIndex);
      pl=freqColl_[parametrizableIndex]->getParameters();
    }
    else
      if (dynamic_cast<DiscreteDistribution*>(parametrizable)){
        distColl_.addObject(dynamic_cast<DiscreteDistribution*>(parametrizable), parametrizableIndex);
        pl=distColl_[parametrizableIndex]->getIndependentParameters();
      }
      else
        if (dynamic_cast<ParametrizableTree*>(parametrizable)){
          treeColl_.addObject(dynamic_cast<ParametrizableTree*>(parametrizable), parametrizableIndex);
          pl=treeColl_[parametrizableIndex]->getParameters();
        }

        else
          throw Exception("Unknown parametrizable object in SubstitutionProcessCollection::addParametrizable.");

  for (size_t i=0; i<pl.size();i++)
    pl[i].setName(pl[i].getName()+"_"+TextTools::toString(parametrizableIndex));
  
  addParameters_(pl);
}

ParameterList SubstitutionProcessCollection::getNonDerivableParameters() const
{
  ParameterList pl=distColl_.getIndependentParameters();
  pl.addParameters(modelColl_.getIndependentParameters());
  pl.addParameters(freqColl_.getIndependentParameters());

  return pl;
}
  


void SubstitutionProcessCollection::fireParameterChanged(const ParameterList& parameters)
{
  modelColl_.matchParametersValues(parameters);
  const vector<size_t>& vM=modelColl_.hasChanged();
  for (size_t i=0; i<vM.size(); i++)
  {
    const vector<size_t>& vs=mModelToSubPro_[vM[i]];
    for (size_t j=0; j<vs.size(); j++)
      vSubProcess_[vs[j]]->changedModel(vM[i]);
  }
  
  freqColl_.matchParametersValues(parameters);
  
  vector<bool> toFire(vSubProcess_.size(), false);

  distColl_.matchParametersValues(parameters);
  
  const vector<size_t>& vD=distColl_.hasChanged();
  for (size_t i=0; i<vD.size(); i++)
  {
    const vector<size_t>& vs=mDistToSubPro_[vD[i]];
    for (size_t j=0; j<vs.size(); j++)
      toFire[vs[j]]=true;
  }

  treeColl_.matchParametersValues(parameters);
  const vector<size_t>& vT=treeColl_.hasChanged();
  for (size_t i=0; i<vT.size(); i++)
  {
    const vector<size_t>& vs=mTreeToSubPro_[vT[i]];      
    for (size_t j=0; j<vs.size(); j++)
      toFire[vs[j]]=true;
  }

  for (size_t j=0; j<toFire.size();j++)
    if (toFire[j])
      vSubProcess_[j]->fireParameterChanged(parameters);

}


void SubstitutionProcessCollection::addSubstitutionProcess(std::map<size_t, std::vector<int> > mModBr, size_t nTree, size_t nRate, size_t nFreq)
{
  if (!freqColl_.hasObject(nFreq))
    throw BadIntegerException("Wrong Root Frequencies Set number",(int)nFreq);
  if (!treeColl_.hasObject(nTree))
    throw BadIntegerException("Wrong Tree number",(int)nTree);
  if (!distColl_.hasObject(nRate))
    throw BadIntegerException("Wrong Rate distribution number",(int)nRate);

  SubstitutionProcessCollectionMember* nSMS=new SubstitutionProcessCollectionMember(this, nTree, nRate);
  nSMS->setRootFrequencies(nFreq);

  std::map<size_t, std::vector<int> >::iterator it;
  for (it=mModBr.begin(); it!=mModBr.end(); it++){
    nSMS->addModel(it->first, it->second);
    mModelToSubPro_[it->first].push_back(vSubProcess_.size());
  }
  
  nSMS->isFullySetUp();

  mTreeToSubPro_[nTree].push_back(vSubProcess_.size());
  mDistToSubPro_[nRate].push_back(vSubProcess_.size());
  mFreqToSubPro_[nFreq].push_back(vSubProcess_.size());
  
  vSubProcess_.push_back(nSMS);
}

void SubstitutionProcessCollection::addSubstitutionProcess(std::map<size_t, std::vector<int> > mModBr, size_t nTree, size_t nRate)
{
  if (!treeColl_.hasObject(nTree))
    throw BadIntegerException("Wrong Tree number",(int)nTree);
  if (!distColl_.hasObject(nRate))
    throw BadIntegerException("Wrong Rate distribution number",(int)nRate);

  SubstitutionProcessCollectionMember* nSMS=new SubstitutionProcessCollectionMember(this, nTree, nRate);

  std::map<size_t, std::vector<int> >::iterator it;
  for (it=mModBr.begin(); it!=mModBr.end(); it++){
    nSMS->addModel(it->first, it->second);
    mModelToSubPro_[it->first].push_back(vSubProcess_.size());
  }

  nSMS->isFullySetUp();
  
  mTreeToSubPro_[nTree].push_back(vSubProcess_.size());
  mDistToSubPro_[nRate].push_back(vSubProcess_.size());
  
  vSubProcess_.push_back(nSMS);
}

bool SubstitutionProcessCollection::hasBranchLengthsParameter(const std::string& name) const
{
  return treeColl_.hasParameter(name);
}

ParameterList SubstitutionProcessCollection::getSubstitutionProcessParameters() const
{
  ParameterList pl=modelColl_.getParameters();

  pl.addParameters(freqColl_.getParameters());
  pl.addParameters(distColl_.getParameters());

  return pl;
}

SubstitutionProcess* SubstitutionProcessCollection::getSubstitutionProcess(size_t  i)
{
  return dynamic_cast<SubstitutionProcess*>(vSubProcess_[i]);
}

const SubstitutionProcess* SubstitutionProcessCollection::getSubstitutionProcess(size_t  i) const
{
  return dynamic_cast<const SubstitutionProcess*>(vSubProcess_[i]);
}
