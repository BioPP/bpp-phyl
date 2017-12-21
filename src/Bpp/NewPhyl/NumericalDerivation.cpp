//
// File: NumericalDerivation.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-12-19
// Last modified: 2017-12-19
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

#include <Bpp/NewPhyl/LinearAlgebra.h>
#include <Bpp/NewPhyl/NumericalDerivation.h>

namespace bpp {
namespace DF {

	template <typename T> class NumericalDerivationShiftDelta : public Value<T> {
	public:
		using Dependencies = FunctionOfValues<double, T>;

    NumericalDerivationShiftDelta (NodeRefVec && deps);


	private:
		int n_;

		void compute () override final {
			callWithValues (*this, [this](T & result, double & delta, const T & arg) {
				result = static_cast<double> (this->n_ * delta) * arg;
			});
		}
	};

	template <> class NumericalDerivationShiftDelta<double>;
	template <> class NumericalDerivationShiftDelta<VectorDouble>;
	template <> class NumericalDerivationShiftDelta<MatrixDouble>;
} // namespace DF
} // namespace bpp
