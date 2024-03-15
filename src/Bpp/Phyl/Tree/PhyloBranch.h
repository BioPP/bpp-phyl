// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_PHYLOBRANCH_H
#define BPP_PHYL_TREE_PHYLOBRANCH_H

#include <Bpp/Clonable.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Utils/MapTools.h>

#include "PhyloTreeExceptions.h"

namespace bpp
{
class PhyloBranchParam;

class PhyloBranch
{
protected:
  bool isLengthDefined_;
  double length_;
  mutable std::map<std::string, Clonable*> properties_;

public:
  /**
   * @brief Constructors.
   *
   * @warning phyloTree_ does not know the edge exists.
   *
   */

  PhyloBranch() :
    isLengthDefined_(false),
    length_(0),
    properties_()
  {}

  PhyloBranch(double length) :
    isLengthDefined_(true),
    length_(length),
    properties_()
  {}

  PhyloBranch(const PhyloBranchParam& branch);

  /**
   * @brief Copy constructor.
   *
   * @param branch The branch to copy.
   */
  PhyloBranch(const PhyloBranch& branch);

  /**
   * @brief Assignation operator.
   *
   * @param branch the branch to copy.
   * @return A reference toward this branch.
   */
  PhyloBranch& operator=(const PhyloBranch& branch);

  PhyloBranch* clone() const { return new PhyloBranch(*this); }

  /**
   * @brief destructor. In Graph, nothing is changed.
   *
   */
  virtual ~PhyloBranch()
  {
    deleteProperties();
  }

  /**
   * @brief Has the length been set?
   * @return true if a length has been defined
   */
  bool hasLength() const
  {
    return isLengthDefined_;
  }


  /**
   * @brief Delete length
   */
  void deleteLength()
  {
    isLengthDefined_ = false;
  }

  /**
   * @brief What is the branch length?
   * @return a double representing the branch length, 0 if length is
   * not defined.
   *
   */
  double getLength() const
  {
    if (!isLengthDefined_)
      throw PhyloBranchPException("PhyloBranch::getLength: length is not defined.", this);

    return length_;
  }


  /**
   * Define a new branch length
   * @param newLength a double repserenting the new length of the branch
   */
  void setLength(double newLength)
  {
    length_ = newLength;
    isLengthDefined_ = true;
  }

  /**
   * @brief Set/add a branch property.
   *
   * If no property with the same name is found, the new property will be added to the list.
   * Conversely, the property will be deleted and replaced by the new one.
   * If you want to keep a copy of the old property, consider using the removeProperty function before.
   *
   * @param name The name of the property to set.
   * @param property The property object (will be cloned).
   */
  void setProperty(const std::string& name, const Clonable& property)
  {
    if (hasProperty(name))
      delete properties_[name];
    properties_[name] = property.clone();
  }

  Clonable* getProperty(const std::string& name)
  {
    if (hasProperty(name))
      return properties_[name];
    else
      throw PhyloBranchPropertyNotFoundException("", name, this);
  }

  const Clonable* getProperty(const std::string& name) const
  {
    if (hasProperty(name))
      return const_cast<const Clonable*>(properties_[name]);
    else
      throw PhyloBranchPropertyNotFoundException("", name, this);
  }

  Clonable* removeProperty(const std::string& name)
  {
    if (hasProperty(name))
    {
      Clonable* removed = properties_[name];
      properties_.erase(name);
      return removed;
    }
    else
      throw PhyloBranchPropertyNotFoundException("", name, this);
  }

  void deleteProperty(const std::string& name)
  {
    if (hasProperty(name))
    {
      delete properties_[name];
      properties_.erase(name);
    }
    else
      throw PhyloBranchPropertyNotFoundException("", name, this);
  }

  /**
   * @brief Remove all branch properties.
   *
   * Attached objects will not be deleted.
   */
  void removeProperties()
  {
    properties_.clear();
  }

  /**
   * @brief Delete all branch properties.
   */
  void deleteProperties()
  {
    for (std::map<std::string, Clonable*>::iterator i = properties_.begin(); i != properties_.end(); i++)
    {
      delete i->second;
    }
    properties_.clear();
  }

  bool hasProperty(const std::string& name) const { return properties_.find(name) != properties_.end(); }

  std::vector<std::string> getPropertyNames() const { return MapTools::getKeys(properties_); }

  bool hasBootstrapValue() const
  {
    return properties_.find("bootstrap") != properties_.end();
  }


  double getBootstrapValue() const
  {
    if (!hasBootstrapValue())
      throw PhyloBranchPropertyNotFoundException("", "bootstrap", this);

    return (dynamic_cast<const Number<double>*>(properties_.find("bootstrap")->second))->getValue();
  }
}; // end of class PhyloBranch
} // end of namespace bpp.
#endif // BPP_PHYL_TREE_PHYLOBRANCH_H
