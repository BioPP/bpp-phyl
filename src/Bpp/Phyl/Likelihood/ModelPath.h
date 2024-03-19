// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_LIKELIHOOD_MODELPATH_H
#define BPP_PHYL_LIKELIHOOD_MODELPATH_H

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/VectorTools.h>

#include "../Model/MixedTransitionModel.h"

namespace bpp
{
/**
 * @brief Organization of submodels in mixed substitution models in a
 * path. See class ModelScenario for a thorough description.
 *
 * A path is defined through an hypergraph, ie a list of hypernodes.
 */
class ModelPath
{
public:
  /**
   * @brief A vector<int> where all elements are different and in
   * INCREASING ORDER. So inclusion should be done through dedicated
   * methods.
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
     */
    PathNode& operator-=(const PathNode& n)
    {
      removeN(n);
      return *this;
    }

    /**
     * @brief checks if this PathNode is included in another one.
     */
    bool operator<=(const PathNode&) const;

    /**
     * @brief checks if this PathNode includes another one.
     */
    bool operator>=(const PathNode&) const;

    /**
     * @brief checks if this PathNode intersects another one.
     */
    bool intersects(const PathNode&) const;

    /**
     * @brief Output
     */
    std::string to_string() const;
  };

private:
  /**
   * Which submodels of MixedTransitionModels are used.
   */
  std::map<std::shared_ptr<MixedTransitionModelInterface>, PathNode> mModPath_;

  /**
   * The leading model (ie which that decides of the submodels probabilities).
   */
  std::shared_ptr<MixedTransitionModelInterface> leadMod_;

  /**
   * @brief probability of this ModelPath.
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

  /**
   * @brief gets the leader model.
   */
  std::shared_ptr<MixedTransitionModelInterface> getLeadModel()
  {
    return leadMod_;
  }

  const std::shared_ptr<MixedTransitionModelInterface> getLeadModel() const
  {
    return leadMod_;
  }

  /**
   * @brief sets the leader model.
   *
   * The model must be in the map before.
   */
  void setLeadModel(std::shared_ptr<MixedTransitionModelInterface> model)
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
  void setModel(std::shared_ptr<MixedTransitionModelInterface> mMod, const Vuint& vnS);

  /**
   * @brief change from a model to another
   *
   * @param mMod1 the ancient mixed model
   * @param mMod2 the new mixed model
   *
   * if mMod1 is the lead model, mMod2 becomes the lead model
   */
  void changeModel(std::shared_ptr<MixedTransitionModelInterface> mMod1,
      std::shared_ptr<MixedTransitionModelInterface> mMod2);

  /**
   * @brief Remove a model
   *
   * if the model is the leadmodel, the leadmodel is set to 0.
   */
  void removeModel(std::shared_ptr<MixedTransitionModelInterface> mMod)
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
  void addToModel(std::shared_ptr<MixedTransitionModelInterface> mMod, const Vuint& vnS);

  /**
   * @brief Cumulates the PathNodes of the given ModelPath into this one.
   */
  ModelPath& operator+=(const ModelPath&);

  /**
   * @brief Remove from the PathNodes of this objet the matching ones of the ModelPath.
   */
  ModelPath& operator-=(const ModelPath&);

  /**
   * @brief checks if this ModelPath is included in another one.
   *  Which means that all submodels of this path are in the other part.
   */
  bool operator<=(const ModelPath&) const;

  /**
   * @brief checks if this ModelPath includes another one.
   */
  bool operator>=(const ModelPath&) const;

  /**
   * @brief checks if this ModelPath intersects another one. Which
   * means that one submodel explicitely declared in a ModelPath is
   * in the other.
   */
  bool intersects(const ModelPath&) const;

  /**
   * @brief returns the probability
   */
  double getProbability() const {return proba_; }

  /**
   * @brief sets the probability
   */
  void setProbability(double x) { proba_ = x; }

  /**
   * @brief checks if there is a pathnode associated with a model
   */
  bool hasModel(std::shared_ptr<MixedTransitionModelInterface> mMod) const
  { return mModPath_.find(mMod) != mModPath_.end(); }

  bool hasModel(std::shared_ptr<const MixedTransitionModelInterface> mMod) const
  {
    for (const auto& mn:mModPath_)
    {
      if (mn.first == mMod)
        return true;
    }
    return false;
  }

  /**
   * @brief gets the pathnode associated with a model
   */
  const PathNode& getPathNode(std::shared_ptr<MixedTransitionModelInterface> mMod) const
  { return mModPath_.at(mMod); }

  const PathNode& getPathNode(std::shared_ptr<const MixedTransitionModelInterface> mMod) const
  {
    for (const auto& mn:mModPath_)
    {
      if (mn.first == mMod)
        return mn.second;
    }
    throw Exception("ModelPath::getPathNode : unknown model " + mMod->getName());
  }

  /**
   * @brief gets the MixedTransitionModel used in the ModelPath
   */
  std::vector<std::shared_ptr<MixedTransitionModelInterface>> getModels() const;

  /**
   * @brief string description
   */
  std::string toString() const;
};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_MODELPATH_H
