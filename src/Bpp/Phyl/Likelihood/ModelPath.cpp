// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
  mModPath_ = hn.mModPath_;
  proba_ = hn.proba_;

  return *this;
}

void ModelPath::setModel(std::shared_ptr<MixedTransitionModelInterface> mMod, const Vuint& vnS)
{
  if (vnS.size() == 0)
    return;
  mModPath_[mMod] = PathNode();
  mModPath_[mMod].insertN(vnS);

  if (mModPath_[mMod].back() >= mMod->getNumberOfModels())
    throw IndexOutOfBoundsException("ModelPath::setModel. Bad submodel number in mixed model", mModPath_[mMod].back(), 0, mMod->getNumberOfModels() - 1);
}

void ModelPath::changeModel(shared_ptr<MixedTransitionModelInterface> mMod1,
                            shared_ptr<MixedTransitionModelInterface> mMod2)
{
  if (mModPath_.find(mMod1) == mModPath_.end())
    throw Exception("ModelPath::changeModel : Unknown model " + mMod1->getName());

  if (leadMod_ == mMod1)
    leadMod_ = mMod2;

  const auto& np = mModPath_[mMod1];
  setModel(mMod2, np);
  mModPath_.erase(mMod1);
}

void ModelPath::addToModel(std::shared_ptr<MixedTransitionModelInterface> mMod, const Vuint& vnS)
{
  if (mModPath_.find(mMod) == mModPath_.end())
    mModPath_[mMod] = PathNode();

  mModPath_[mMod].insertN(vnS);

  if (mModPath_.size() > 0 && mModPath_[mMod].back() >= mMod->getNumberOfModels())
    throw IndexOutOfBoundsException("ModelPath::addToModel. Bad submodel number in mixed model", mModPath_[mMod].back(), 0, mMod->getNumberOfModels() - 1);
}

bool ModelPath::operator<=(const ModelPath& hn) const
{
  const auto& mpath2 = hn.mModPath_;

  for (const auto& ipath : mModPath_)
  {
    if (mpath2.find(ipath.first) != mpath2.end() &&
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
  const auto& mpath2 = hn.mModPath_;

  for (const auto& ipath : mModPath_)
  {
    if (mpath2.find(ipath.first) == mpath2.end() ||
        ipath.second.intersects(mpath2.at(ipath.first)))
      return true;
  }

  for (const auto& ipath : mpath2)
  {
    if (mModPath_.find(ipath.first) == mModPath_.end())
      return true;
  }

  return false;
}

ModelPath& ModelPath::operator+=(const ModelPath& hn)
{
  const auto& mpath2 = hn.mModPath_;

  for (const auto& ipath : mpath2)
  {
    addToModel(ipath.first, ipath.second);
  }

  return *this;
}

ModelPath& ModelPath::operator-=(const ModelPath& hn)
{
  const auto& mpath2 = hn.mModPath_;

  for (const auto& ipath : mpath2)
  {
    if (mModPath_.find(ipath.first) != mModPath_.end())
    {
      mModPath_[ipath.first] -= ipath.second;
      if (mModPath_[ipath.first].size() == 0)
        mModPath_.erase(ipath.first);
    }
  }

  return *this;
}

vector< shared_ptr<MixedTransitionModelInterface> > ModelPath::getModels() const
{
  vector< shared_ptr<MixedTransitionModelInterface> > models;

  std::transform(
    mModPath_.begin(),
    mModPath_.end(),
    std::back_inserter(models),
    [](const map<shared_ptr<MixedTransitionModelInterface>, PathNode>::value_type& pair){
    return pair.first;
  });

  return models;
}

std::string ModelPath::toString() const
{
  std::string output;
  bool deb = true;
  for (const auto& mod : mModPath_)
  {
    if (!deb)
      output += "&";

    auto model = mod.first;

    if (dynamic_cast<const AbstractBiblioMixedTransitionModel*>(model.get()) == nullptr)
    {
      std::string name = "";
      auto pMS = dynamic_cast<const MixtureOfTransitionModels*>(model.get());
      if (pMS)
      {
        name = "Mixture[";
        bool com = false;
        for (auto nb:mod.second)
        {
          name = name + (com ? ", " : "") + pMS->AbstractMixedTransitionModel::nModel(nb - 1).getName();
        }
        name += "]";
      }
      else
      {
        auto pMT = dynamic_cast<const MixtureOfATransitionModel*>(model.get());
        name = "MixedModel";

        const TransitionModelInterface& eM = pMT->model(0);

        name += "." + eM.getName() + "[" + mod.second.to_string() + "]";
      }
      output += name;
    }
    else
      output += model->getName() + "[" + mod.second.to_string() + "]";

    deb = false;
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

    for ( ; it != end(); it++)
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
                       [vn](const uint x) -> bool {
    return std::find(vn.begin(), vn.end(), x) != vn.end();
  }), end());
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
  bool deb = true;
  for (auto ind:*this)
  {
    if (!deb)
      output += ",";
    output += std::to_string(ind + 1);
    deb = false;
  }
  return output;
}
