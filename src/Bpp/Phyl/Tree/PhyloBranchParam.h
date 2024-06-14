// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_TREE_PHYLOBRANCHPARAM_H
#define BPP_PHYL_TREE_PHYLOBRANCHPARAM_H

#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/NumConstants.h>


namespace bpp
{
class PhyloBranch;

/**
 * Basic branch in which length is parameterized.
 *
 * Note : Branch Length must be the first parameter (index 0).
 *
 **/

class PhyloBranchParam :
  public AbstractParametrizable
{
public:
  /**
   * @brief Constructors.
   *
   * @warning phyloTree_ does not know the edge exists.
   *
   */

  PhyloBranchParam(const std::string& prefix = "") :
    AbstractParametrizable(prefix)
  {
    addParameter_(new Parameter("BrLen", NumConstants::SMALL(), Parameter::R_PLUS_STAR));
  }

  PhyloBranchParam(const PhyloBranch& branch);

  /**
   * @brief Copy constructor.
   *
   * @param branch The branch to copy.
   */
  PhyloBranchParam(const PhyloBranchParam& branch) :
    AbstractParametrizable(branch)
  {}


  /**
   * @brief Assignation operator.
   *
   * @warning This operator copies all fields, excepted father and
   * son branch pointers. Without specific id, the PhyloTree does
   * not know the existence of the new edge.
   *
   * @param branch the branch to copy.
   * @return A reference toward this branch.
   */
  PhyloBranchParam& operator=(const PhyloBranchParam& branch)
  {
    AbstractParametrizable::operator=(* this);
    return *this;
  }


  PhyloBranchParam* clone() const { return new PhyloBranchParam(*this); }


  /**
   * @brief What is the branch length?
   * @return a double representing the branch length, 0 if length is
   * not defined.
   *
   */
  double getLength() const
  {
    return getParameter_(0).getValue();
  }


  /**
   * Define a new branch length
   * @param newLength a double repserenting the new length of the branch
   */
  void setLength(double newLength)
  {
    getParameter_(0).setValue(newLength);
  }


  friend class ParametrizablePhyloTree;
}; // end of class PhyloBranchParam
} // end of namespace bpp.
#endif // BPP_PHYL_TREE_PHYLOBRANCHPARAM_H
