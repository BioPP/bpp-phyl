//
// File: new_template_expr.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-08-30
// Last modified: 2017-08-30
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

#include <iostream>

#include <Bpp/NewPhyl/Debug.h>
#include <Bpp/NewPhyl/TemplateExpression.h>

using bpp::Expr::ref;
using bpp::Expr::constant;

void NOTICEME(int& r, const int& a, const int& b)
{
  r = (ref(a) + ref(b)).compute();
}
// -> est optimisé correctement

TEST_CASE("test")
{
  int a = 42;
  int b = 10;
  int r = 0;
  NOTICEME(r, a, b);
  std::cout << r << '\n';

  auto p = bpp::Expr::make_simplified_assignment<int>(ref(b) + constant(0));
  int r2 = 0;
  p->update(r2);
  std::cout << "r2= " << r2 << '\n';
  std::cout << "r2_expr= " << bpp::demangle(typeid(*p).name()) << '\n';
}
