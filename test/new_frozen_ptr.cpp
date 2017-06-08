//
// File: new_frozen_ptr.cpp
// Authors:
// Created: 2017-06-08
// Last modified: 2017-06-08
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

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

#include <Bpp/NewPhyl/FrozenPtr.h>

struct Base
{
  virtual ~Base() = default;
};
struct Derived : Base
{
};

TEST_CASE("FreezableUniquePtr")
{
  auto empty = bpp::FreezableUniquePtr<int>{};
  CHECK(!empty);

  auto p = bpp::make_freezable_unique<int>(42);
  CHECK(p);
  CHECK(*p == 42);
  *p = 3;
  CHECK(*p == 3);

  auto p2 = std::move(p);
  CHECK(p2);
  CHECK(!p);
  CHECK(*p2 == 3);

  auto derived = bpp::make_freezable_unique<Derived>();
  bpp::FreezableUniquePtr<Base> base{std::move(derived)};
}

TEST_CASE("FrozenSharedPtr")
{
  auto p = bpp::make_freezable_unique<int>(44);
  auto sp = std::move(p).freeze();
  CHECK(!p);
  CHECK(sp);
  CHECK(*sp == 44);

  auto spcpy = sp;
  CHECK(sp);
  CHECK(spcpy);
  CHECK(*sp == 44);
  CHECK(*spcpy == 44);
  CHECK(sp.get() == spcpy.get());

  bpp::FrozenSharedPtr<Base> basep{bpp::make_freezable_unique<Derived>().freeze()};
}
