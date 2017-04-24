//
// File: new_phyl_dataflow.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-19
// Last modified: 2017-04-19
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

#include <Bpp/NewPhyl/DataFlow.h>
#include <fstream>
#include <iostream>

using bpp::DF::Node;

struct A : public Node::Impl
{
  void compute() override {}
  std::string name() const override { return "A"; }

  void addDep(Node n)
  {
    n.get().registerNode(this);
    dependencyNodes_.emplace_back(std::move(n));
  }
};

TEST_CASE("test")
{
  auto a = Node::create<A>();
  auto b = Node::create<A>();
  auto c = Node::create<A>();
  auto d = Node::create<A>();
  dynamic_cast<A&>(a.get()).addDep(b);
  dynamic_cast<A&>(a.get()).addDep(c);
  dynamic_cast<A&>(d.get()).addDep(a);
  auto file = std::ofstream("df_debug");
  bpp::DF::debugDag(file, a);
}
