//
// File: ModelPath.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 27 novembre 2019, à 09h 09
// From: MixedSubstitutionModelSet
//

/*
   Copyright or <A9> or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "ModelPath.h"

using namespace bpp;
using namespace std;

#include <Bpp/Phyl/Model/AbstractBiblioMixedTransitionModel.h>
#include <Bpp/Phyl/Model/MixtureOfTransitionModels.h>
#include <Bpp/Phyl/Model/MixtureOfATransitionModel.h>

ModelPath::ModelPath(const ModelPath& hn) :
  mModPath_(hn.mModPath_),
  proba_(hn.proba_)
{}


ModelPath& ModelPath::operator=(const ModelPath& hn)
{
  mModPath_=hn.mModPath_;
  proba_ = hn.proba_;

  return *this;
}

void ModelPath::setModel(std::shared_ptr<MixedTransitionModel> mMod, const Vuint& vnS)
{
  if (vnS.size()==0)
    return;
  mModPath_[mMod] = PathNode();
  mModPath_[mMod].insertN(vnS);
  
  if (mModPath_[mMod].back() >= mMod->getNumberOfModels())
    throw IndexOutOfBoundsException("ModelPath::setModel. Bad submodel number in mixed model", mModPath_[mMod].back(), 0, mMod->getNumberOfModels()-1);
}

void ModelPath::changeModel(std::shared_ptr<MixedTransitionModel> mMod1,
                            std::shared_ptr<MixedTransitionModel> mMod2)
{
  if (mModPath_.find(mMod1)==mModPath_.end())
    throw Exception("ModelPath::changeModel : Unknown model " + mMod1->getName());

  if (leadMod_==mMod1)
    leadMod_=mMod2;
  
  const auto& np=mModPath_[mMod1];
  setModel(mMod2, np);
  mModPath_.erase(mMod1);
}
    
void ModelPath::addToModel(std::shared_ptr<MixedTransitionModel> mMod, const Vuint& vnS)
{
  if (mModPath_.find(mMod)==mModPath_.end())
    mModPath_[mMod] = PathNode();
  
  mModPath_[mMod].insertN(vnS);
  
  if (mModPath_.size()>0 && mModPath_[mMod].back() >= mMod->getNumberOfModels())
    throw IndexOutOfBoundsException("ModelPath::addToModel. Bad submodel number in mixed model", mModPath_[mMod].back(), 0, mMod->getNumberOfModels()-1);
}

bool ModelPath::operator<=(const ModelPath& hn) const
{
  const auto& mpath2=hn.mModPath_;
  
  for (const auto& ipath : mModPath_)
  {
    if (mpath2.find(ipath.first)!=mpath2.end() &&
        !(ipath.second <= mpath2.at(ipath.first)))
      return false;
  }

  return true;
}

bool ModelPath::operator>=(const ModelPath& hn) const
{
  return hn <= *this;
}

bool ModelPath::intersects(const ModelPath& hn) const
{
  const auto& mpath2=hn.mModPath_;

  for (const auto& ipath : mModPath_)
  {
    if (mpath2.find(ipath.first)==mpath2.end() ||
        ipath.second.intersects(mpath2.at(ipath.first)))
      return true;
  }

  for (const auto& ipath : mpath2)
  {
    if (mModPath_.find(ipath.first)==mModPath_.end())
      return true;
  }

  return false;
}

ModelPath& ModelPath::operator+=(const ModelPath& hn)
{
  const auto& mpath2=hn.mModPath_;
  
  for (const auto& ipath : mpath2)
    addToModel(ipath.first,ipath.second);
  
  return *this;
}

ModelPath& ModelPath::operator-=(const ModelPath& hn)
{
  const auto& mpath2=hn.mModPath_;
  
  for (const auto& ipath : mpath2)
  {
    if (mModPath_.find(ipath.first)!=mModPath_.end())
    {
      mModPath_[ipath.first]-=ipath.second;
      if (mModPath_[ipath.first].size()==0)
        mModPath_.erase(ipath.first);
    }
  }
  
  return *this;
}

std::vector<std::shared_ptr<MixedTransitionModel>> ModelPath::getModels() const
{
  std::vector<std::shared_ptr<MixedTransitionModel>> models;
  
  std::transform(
    mModPath_.begin(),
    mModPath_.end(),
    std::back_inserter(models),
    [](const std::map<std::shared_ptr<MixedTransitionModel>, PathNode>::value_type &pair){return pair.first;});
  
  return models;
};

std::string ModelPath::to_string() const
{
  std::string output;
  bool deb=true;
  for (const auto& mod:mModPath_)
  {
    if (!deb)
      output += "&";

    auto model=mod.first;

    if (dynamic_cast<const AbstractBiblioMixedTransitionModel*>(model.get()) == NULL)
    {
      std::string name="";
      auto pMS = dynamic_cast<const MixtureOfTransitionModels*>(model.get());
      if (pMS)
      {
        name= "Mixture[";
        bool com=false;
        for (auto nb:mod.second)
        {
          name = name + (com?", ":"") + pMS->AbstractMixedTransitionModel::getNModel(nb-1)->getName();
        }
        name += "]";
      }
      else
      {
        auto pMT = dynamic_cast<const MixtureOfATransitionModel*>(model.get());
        name = "MixedModel";
        
        const TransitionModel* eM = pMT->getModel(0);

        name += "." + eM->getName() + "[" + mod.second.to_string() + "]";
      }
      output += name;
    }
    else
      output += model->getName() + "[" + mod.second.to_string() + "]";
    
    deb=false;
  }
  return output;
}

/**********************************************************/
/******************** NODE ********************************/
/***********************************************************/

void ModelPath::PathNode::insertN(const Vuint& vn)
{
  for (const auto& it2 : vn)
  {
    vector<uint>::const_iterator it(begin());
    
    for (; it != end(); it++)
    {
      if (*it >= it2)
        break;
    }
    if (it == end())
      push_back(it2);
    else if (*it != it2)
      insert(it, it2);
  }
}

void ModelPath::PathNode::removeN(const Vuint& vn)
{
  erase(std::remove_if(begin(), end(),
                       [vn](const uint x) -> bool {return std::find(vn.begin(), vn.end(), x)!=vn.end();}),end());
}

bool ModelPath::PathNode::operator<=(const PathNode& n) const
{
  vector<uint>::const_iterator it2(n.begin());

  for (const auto& it : *this)
  {
    while (it2 != n.end()  && (*it2 < it))
      it2++;
    if (it2 == n.end() || (*it2 > it))
      return false;
  }
  return true;
}

bool ModelPath::PathNode::operator>=(const PathNode& n) const
{
  return n <= *this;
}

bool ModelPath::PathNode::intersects(const PathNode& n) const
{
  vector<uint>::const_iterator it2(n.begin());

  for (const auto& it : *this)
  {
    while (it2 != n.end()  && (*it2 < it))
      it2++;

    if (it2 == n.end())
      return false;
    if (*it2 == it)
      return true;
  }
  return false;
}

std::string ModelPath::PathNode::to_string() const
{
  std::string output;
  bool deb=true;
  for (auto ind:*this)
  {
    if (!deb)
      output += ",";
    output += std::to_string(ind+1);
    deb=false;
  }
  return output;
}
