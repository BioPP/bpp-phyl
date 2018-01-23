//
// File: new_dataflow.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-02-23 00:00:00
// Last modified: 2017-06-08
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

  This software is a computer program whose purpose is to provide classes
  for numerical calculus. This file is part of the Bio++ project.

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

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <Bpp/NewPhyl/DataFlowInternal.h>
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/NumericalDerivation.h>
#include <fstream>

using namespace bpp::DF;

struct NumericallyDerivable : public Value<double>
{
  using Dependencies = FunctionOfValues<double, double>;

  NumericallyDerivable(NodeRefVec&& deps)
    : Value<double>(std::move(deps))
  {
    checkDependencies(*this);
  }

  void compute() final
  {
    callWithValues(*this, [](double& r, double a, double b) { r = a * a + b * b; });
  }
};

TEST_CASE("shift delta")
{
  auto delta = makeNode<Constant<double>>(0.0001);
  auto x = makeNode<Mutable<double>>(1);

  auto shift_0 = makeNode<NumericalDerivationShiftDelta<double>>({delta, x}, 0);
  CHECK(shift_0 == x); // Simplification

  auto shift_1 = makeNode<NumericalDerivationShiftDelta<double>>({delta, x}, 1);
  CHECK(shift_1->getValue() == doctest::Approx(x->getValue() + delta->getValue()));

  auto shift_2 = makeNode<NumericalDerivationShiftDelta<double>>({delta, shift_1}, 1);
  CHECK(shift_2->dependency(1) == x); // Simplification

  // Print DF graph
  std::ofstream fd("df_debug");
  bpp::debugDag(fd, x, DebugOptions::FollowUpwardLinks | DebugOptions::DetailedNodeInfo);
}

TEST_CASE("AAA")
{
  using namespace bpp::DF;

  auto a = makeNode<Mutable<double>>(1);
  auto b = makeNode<Mutable<double>>(-1);
  auto f = makeNode<NumericallyDerivable>({a, b});

  // Print DF graph
  std::ofstream fd("df_debug");
  bpp::debugDag(fd, f);
}
