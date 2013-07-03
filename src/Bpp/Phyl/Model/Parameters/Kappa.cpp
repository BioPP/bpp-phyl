//
// File: Kappa.cpp
// Created by: Laurent Guéguen
// Created on: mardi 16 avril 2013, à 22h 18
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "Kappa.h"

#include "../../Mapping/SubstitutionRegister.h"

#include <Bpp/Numeric/Matrix/Matrix.h>


using namespace bpp;

using namespace std;

/******************************************************************************/

Kappa::Kappa(const Alphabet* alpha, const string& name, double kappa) :
  RegisteredParameter(alpha, name, kappa, &Parameter::R_PLUS)
{
  RowMatrix<size_t> mat(getAlphabet()->getSize(), getAlphabet()->getSize());
  if (dynamic_cast<const NucleicAlphabet*>(alpha))
  {
    for (size_t i = 0; i < 4; i++)
    {
      for (size_t j = 0; j < 4; j++)
      {
        switch (i - j)
        {
        case 0:
          mat(i, j) = 0;
          break;
        case 2:
        case -2:
          mat(i, j) = 1;
          break;
        default:
          mat(i, j) = 2;
        }
      }
    }
  }
  else if (dynamic_cast<const CodonAlphabet*>(alpha))
  {
    for (size_t i = 0; i < 64; i++)
    {
      for (size_t j = 0; j < 64; j++)
      {
        if (i == j)
        {
          mat(i, j) = 0;
          continue;
        }
        bool flag = false;
        size_t i2 = i;
        size_t j2 = j;
        unsigned int ph;
        for (ph = 0; ph < 3; ph++)
        {
          if (i2 % 4 != j2 % 4)
          {
            if (flag)
              break;
            else
            {
              flag = true;
            }
          }
          i2 /= 4;
          j2 /= 4;
        }
        if (ph == 3)
          mat(i, j) = (abs(static_cast<int>(i) - static_cast<int>(j)) == 2 ? 1 : 2);
      }
    }
  }
  else
    throw Exception("Bad alphabet for Kappa parameter: " + alpha->getAlphabetType());

  setSubstitutionRegister(new GeneralSubstitutionRegister(alpha, mat));
}

double Kappa::computeEstimate(std::vector<double> values) const
{
  if (values[1] == 0.)
    return 0.;
  else
    return values[0] / values[1];
}

