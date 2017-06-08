//
// File: new_range.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-05 00:00:00
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

#include <Bpp/NewPhyl/Range.h>
#include <algorithm>
#include <list>
#include <string>
#include <type_traits>
#include <vector>

TEST_CASE("type deduction")
{
  std::vector<int> mut{1, 2, 3, 4};
  auto r_mut = bpp::range(mut);
  using MutRange_Is_MutableIt = std::is_same<bpp::Range::Range<std::vector<int>::iterator>, decltype(r_mut)>;
  CHECK(MutRange_Is_MutableIt::value);

  const std::vector<int> konst{1, 2, 3, 4};
  auto r_konst = bpp::range(konst);
  using KonstRange_Is_ConstIt = std::is_same<bpp::Range::Range<std::vector<int>::const_iterator>, decltype(r_konst)>;
  CHECK(KonstRange_Is_ConstIt::value);

  auto r_int = bpp::range(int(42));
  using IntRange_Is_IntIt = std::is_same<bpp::Range::Range<bpp::Iterator::Integer<int>>, decltype(r_int)>;
  CHECK(IntRange_Is_IntIt::value);
}

TEST_CASE("empty vector")
{
  std::vector<int> empty_vect;
  auto r = bpp::range(empty_vect);
  CHECK(r.empty());
  CHECK(r.size() == 0);
  CHECK_FALSE(r.contains(empty_vect.begin()));
  CHECK_FALSE(r.contains(empty_vect.end()));
}

TEST_CASE("filled vector")
{
  std::vector<int> v1234 = {1, 2, 3, 4};
  auto r = bpp::range(v1234);
  CHECK_FALSE(r.empty());
  CHECK(r.front() == v1234.front());
  CHECK(r.back() == v1234.back());
  CHECK(r.size() == v1234.size());
  CHECK(r.contains(v1234.begin()));
  CHECK(r.contains(v1234.begin() + 2));
  CHECK_FALSE(r.contains(v1234.end()));
  CHECK(std::equal(v1234.begin(), v1234.end(), r.begin()));
  CHECK(*r.at(1) == 2);
  CHECK(*r.at(-3) == 2);
  // Modify
  r[2] = 42;
  CHECK(std::equal(v1234.begin(), v1234.end(), r.begin()));
  // Slice
  auto slice = r.slice(1, 3);
  CHECK(std::equal(v1234.begin() + 1, v1234.begin() + 3, slice.begin()));
  slice[1] = 3;
  CHECK(v1234[2] == 3);
  // Slice from end
  auto rslice = r.slice(1, -1);
  CHECK(std::equal(slice.begin(), slice.end(), rslice.begin()));
}

TEST_CASE("between vector, list, string")
{
  // Test exchange between list, string, and others. List has more costly and restricted api (bidir
  // only ops).
  std::string s("hello world");
  auto v = bpp::range(s).to_container<std::vector<char>>();
  CHECK(std::equal(v.begin(), v.end(), s.begin()));
  auto l = bpp::range(v).pop_back(6).to_container<std::list<char>>();
  auto s2 = bpp::range(l).to_container<std::string>();
  CHECK(s2 == "hello");
}

TEST_CASE("integer range")
{
  auto r = bpp::range(10);
  CHECK(r.size() == 10);
  auto v = r.slice(4, -4).to_container<std::vector<int>>();
  CHECK(v.size() == 2);
  CHECK(v[0] == 4);
  CHECK(v[1] == 5);
}
