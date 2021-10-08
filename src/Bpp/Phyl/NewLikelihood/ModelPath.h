//
// File: ModelPath.h
// Created by: Laurent Guéguen
// Created on: mardi 26 novembre 2019, à 13h 42
// From: MixedSubstitutionModelSet
//

/*
   Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _MODEL_PATH_H_
#define _MODEL_PATH_H_

#include "../Model/MixedTransitionModel.h"

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/VectorTools.h>


namespace bpp
{
/**
 * @brief Organization of submodels in mixed substitution models in a
 * path. See class ModelScenario for a thorough description.
 *
 * A path is defined through an hypergraph, ie a list of hypernodes.
 *
 */

class ModelPath
{
public:
  /**
   * @brief A vector<int> where all elements are different and in
   * INCREASING ORDER. So inclusion should be done through dedicated
   * methods.
   *
   */

  class PathNode : public Vuint
  {
public:
    PathNode() {}
    PathNode(const PathNode& n) : Vuint(n){}

    ~PathNode(){}

    /**
     * @brief Insert elements
     *
     */

    void insertN(const Vuint& vn);

    /**
     * @brief Remove elements
     *
     */

    void removeN(const Vuint& vn);


    PathNode& operator=(const Vuint& vn)
    {
      clear();
      insertN(vn);
      return *this;
    }

    /**
     * @brief Cumulates the elements of the given PathNode into this one.
     *
     */
    PathNode& operator+=(const PathNode& n)
    {
      insertN(n);
      return *this;
    }

    /**
     * @brief Remove the elements of the given PathNode from this one.
     *
     */
    PathNode& operator-=(const PathNode& n)
    {
      removeN(n);
      return *this;
    }

    /**
     * @brief checks if this PathNode is included in another one.
     *
     */

    bool operator<=(const PathNode&) const;

    /**
     * @brief checks if this PathNode includes another one.
     *
     */

    bool operator>=(const PathNode&) const;

    /**
     * @brief checks if this PathNode intersects another one.
     *
     */

    bool intersects(const PathNode&) const;

    /**
     * @brief Output
     *
     */

    std::string to_string() const;
  };

private:
  /*
   * Which submodels of MixedTransitionModels are used.
   *
   */

  std::map<std::shared_ptr<MixedTransitionModel>, PathNode> mModPath_;

  /*
   * The leading model (ie which that decides of the submodels probabilities).
   *
   */

  std::shared_ptr<MixedTransitionModel> leadMod_;

  /**
   * @brief probability of this ModelPath.
   *
   */

  double proba_;

public:
  ModelPath() : mModPath_(), leadMod_(0), proba_(0) {}
  ModelPath(const ModelPath&);
  ModelPath& operator=(const ModelPath&);
  ~ModelPath(){}

  size_t size() const
  {
    return mModPath_.size();
  }

  /*
   * @brief gets the leader model.
   *
   */
  std::shared_ptr<MixedTransitionModel> getLeadModel()
  {
    return leadMod_;
  }

  const std::shared_ptr<MixedTransitionModel> getLeadModel() const
  {
    return leadMod_;
  }

  /*
   * @brief sets the leader model.
   *
   * The model must be in the map before.
   *
   */
  void setLeadModel(std::shared_ptr<MixedTransitionModel> model)
  {
    if (mModPath_.find(model) == mModPath_.end())
      throw Exception("ModelPath::setLeadModel: Unknown model " + model->getName());
    leadMod_ = model;
  }

  /**
   * @brief sets submodel numbers in the mixed model. Checks if all
   *  the numbers are valid.
   *
   * @param mMod the mixed model
   * @param vnS vector of indexes of the submodels
   */

  void setModel(std::shared_ptr<MixedTransitionModel> mMod, const Vuint& vnS);

  /**
   * @brief change from a model to another
   *
   * @param mMod1 the ancient mixed model
   * @param mMod2 the new mixed model
   *
   * if mMod1 is the lead model, mMod2 becomes the lead model
   *
   */

  void changeModel(std::shared_ptr<MixedTransitionModel> mMod1,
                   std::shared_ptr<MixedTransitionModel> mMod2);

  /**
   * @brief Remove a model
   *
   * if the model is the leadmodel, the leadmodel is set to 0.
   */
  void removeModel(std::shared_ptr<MixedTransitionModel> mMod)
  {
    if (mMod == leadMod_)
      leadMod_ = 0;

    if (mModPath_.find(mMod) != mModPath_.end())
      mModPath_.erase(mMod);
  }

  /**
   * @brief adds submodel numbers to the mixed model. Checks if all
   *  the numbers are valid.
   *
   * @param mMod the mixed model
   * @param vnS vector of numbers of the submodel
   */

  void addToModel(std::shared_ptr<MixedTransitionModel> mMod, const Vuint& vnS);

  /**
   * @brief Cumulates the PathNodes of the given ModelPath into this one.
   *
   */

  ModelPath& operator+=(const ModelPath&);

  /**
   * @brief Remove from the PathNodes of this objet the matching ones of the ModelPath.
   *
   */

  ModelPath& operator-=(const ModelPath&);

  /**
   * @brief checks if this ModelPath is included in another one.
   *  Which means that all submodels of this path are in the other part.
   *
   */

  bool operator<=(const ModelPath&) const;

  /**
   * @brief checks if this ModelPath includes another one.
   *
   */

  bool operator>=(const ModelPath&) const;

  /**
   * @brief checks if this ModelPath intersects another one. Which
   * means that one submodel explicitely declared in a ModelPath is
   * in the other.
   *
   */

  bool intersects(const ModelPath&) const;

  /**
   * @brief returns the probability
   *
   */
  double getProbability() const {return proba_; }

  /**
   * @brief sets the probability
   *
   */
  void setProbability(double x) { proba_ = x; }

  /**
   * @brief checks if there is a pathnode associated with a model
   *
   *
   */
  bool hasModel(std::shared_ptr<MixedTransitionModel> mMod) const
  { return mModPath_.find(mMod) != mModPath_.end(); }

  bool hasModel(const MixedTransitionModel* mMod) const
  {
    for (const auto& mn:mModPath_)
    {
      if (mn.first.get() == mMod)
        return true;
    }
    return false;
  }

  /**
   * @brief gets the pathnode associated with a model
   *
   *
   */
  const PathNode& getPathNode(std::shared_ptr<MixedTransitionModel> mMod) const
  { return mModPath_.at(mMod); }

  const PathNode& getPathNode(const MixedTransitionModel* mMod) const
  {
    for (const auto& mn:mModPath_)
    {
      if (mn.first.get() == mMod)
        return mn.second;
    }
    throw Exception("ModelPath::getPathNode : unknown model " + mMod->getName());
  }

  /**
   * @brief gets the MixedTransitionModel used in the ModelPath
   *
   */

  std::vector<std::shared_ptr<MixedTransitionModel> > getModels() const;

  /**
   * @brief string description
   *
   */

  std::string to_string() const;
};
} // end of namespace bpp.

#endif// _MODEL_PATH_H_
