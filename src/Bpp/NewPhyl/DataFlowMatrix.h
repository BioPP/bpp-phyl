//
// File: DataFlowMatrix.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-11
// Last modified: 2017-10-11
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

#ifndef BPP_NEWPHYL_DATAFLOWMATRIX_H
#define BPP_NEWPHYL_DATAFLOWMATRIX_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Signed.h>
#include <Eigen/Core>

namespace bpp {
namespace DF {
	// Explicit template declaration TODO test usefulness
	// extern template class Value<Eigen::MatrixXd>;

	using VectorDouble = Eigen::VectorXd;
	using MatrixDouble = Eigen::MatrixXd;

	struct ConstantMatrixDouble : public Value<MatrixDouble> {
		template <typename Derived>
		ConstantMatrixDouble (const Eigen::EigenBase<Derived> & expr)
		    : Value<MatrixDouble> (noDependency, expr) {}
		void compute () override final;
		NodeRef derive (const Node & node) override final;
		template <typename Derived>
		static std::shared_ptr<ConstantMatrixDouble> create (const Eigen::EigenBase<Derived> & expr) {
			return std::make_shared<ConstantMatrixDouble> (expr);
		}
	};

	struct AddMatrixDouble : public Value<MatrixDouble> {
		using Dependencies = ReductionOfValue<MatrixDouble>;
		AddMatrixDouble (NodeRefVec && deps, SizeType nbRows, SizeType nbCols);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
		static ValueRef<MatrixDouble> create (NodeRefVec && deps, SizeType nbRows, SizeType nbCols);
	};

	struct MulMatrixDouble : public Value<MatrixDouble> {
		using Dependencies = FunctionOfValues<MatrixDouble, MatrixDouble>;
		MulMatrixDouble (NodeRefVec && deps, SizeType nbRows, SizeType nbCols);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
		static ValueRef<MatrixDouble> create (NodeRefVec && deps, SizeType nbRows, SizeType nbCols);
	};

	// TODO
	struct CWiseMulMatrixDouble : public Value<MatrixDouble> {
		using Dependencies = ReductionOfValue<MatrixDouble>;
		CWiseMulMatrixDouble (NodeRefVec && deps, SizeType nbRows, SizeType nbCols);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
		static ValueRef<MatrixDouble> create (NodeRefVec && deps, SizeType nbRows, SizeType nbCols);
	};
} // namespace DF
} // namespace bpp

#endif // BPP_NEWPHYL_DATAFLOWMATRIX_H
