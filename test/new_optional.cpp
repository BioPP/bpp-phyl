//
// File: new_optional.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-22 00:00:00
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

#include <Bpp/NewPhyl/Optional.h>
#include <map>
#include <memory> // tests with uniq_ptr
#include <string> // to_string
#include <type_traits>

TEST_CASE("basic operations")
{
  bpp::Optional<int> a;
  CHECK_FALSE(a);

  a = 42;
  CHECK(a);
  CHECK(a.value() == 42);

  a.reset();
  CHECK_FALSE(a);

  a = 4;
  CHECK(a);
  CHECK(*a == 4);
  a = bpp::nullopt;
  CHECK(!a);

  // Copy construction
  a = 54;
  bpp::Optional<int> b{a};
  CHECK(b);
  CHECK(*b == 54);

  // Copy assign
  b = -1;
  a.reset();
  CHECK(!a);
  CHECK(b);
  CHECK(*b == -1);
  a = b;
  CHECK(a);
  CHECK(*a == -1);

  // Null optional copy
  a = bpp::nullopt;
  CHECK(!a);
  b = a;
  CHECK(!b);

  bpp::Optional<int> c{b};
  CHECK(!c);
}

// Non movable nor copyable type, shoud support emplace stuff
struct OnlyConstructible
{
  int a;
  OnlyConstructible(int i)
    : a(i)
  {
  }
  OnlyConstructible(const OnlyConstructible&) = delete;
};

TEST_CASE("non move/copy objects")
{
  bpp::Optional<OnlyConstructible> a{bpp::in_place, 32};
  CHECK(a);
  CHECK(a->a == 32);

  a.emplace(12);
  CHECK(a);
  CHECK(a->a == 12);
}

// Movable only type
TEST_CASE("move only objects")
{
  using UPI = std::unique_ptr<int>;
  bpp::Optional<UPI> p;
  CHECK(!p);

  // Move value in
  p = UPI{new int{42}};
  CHECK(p);
  CHECK(**p == 42);

  // Move construct
  bpp::Optional<UPI> p2{std::move(p)};
  CHECK(p);   // Still defined
  CHECK(!*p); // Moved from
  CHECK(p2);
  CHECK(**p2 == 42);

  // Move assign
  p = std::move(p2);
  CHECK(p2);   // Defined
  CHECK(!*p2); // Moved from
  CHECK(p);
  CHECK(**p == 42);

  // Move value out
  UPI up = *std::move(p);
  CHECK(p);
  CHECK(!*p);
  CHECK(up);
  CHECK(*up == 42);

  // Move value in
  p = std::move(up);
  CHECK(!up);
  CHECK(*p);
  CHECK(**p == 42);
}

TEST_CASE("constness")
{
  // Can change value of const Opt<T>
  const bpp::Optional<int> const_opt{23};
  *const_opt = 42;
  CHECK(*const_opt == 42);

  // Can modify the value only by reconstructing it for Opt<const T>
  bpp::Optional<const int> opt_const;
  CHECK(!opt_const);
  opt_const = 42; // Automatically uses reconstruction
  CHECK(opt_const);
  CHECK(*opt_const == 42);
  opt_const.emplace(33);
  CHECK(*opt_const == 33);
  opt_const.reset();
  CHECK(!opt_const);
}

TEST_CASE("reference optionals")
{
  int a = 42;
  bpp::Optional<int&> ref{a};
  CHECK(ref);
  CHECK(*ref == a);
  *ref = 32;
  CHECK(a == 32);
  ref = nullptr;
  CHECK(!ref);
  ref = a;
  CHECK(ref);
  CHECK(&ref.value() == &a);
  ref = bpp::nullopt;
  CHECK(!ref);
  ref = &a;
  CHECK(ref);
  CHECK(&ref.value() == &a);
}
