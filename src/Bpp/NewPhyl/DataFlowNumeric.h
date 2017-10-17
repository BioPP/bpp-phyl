//
// File: DataFlowNumeric.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-09-15 00:00:00
// Last modified: 2017-10-10
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

#ifndef BPP_NEWPHYL_DATAFLOWNUMERIC_H
#define BPP_NEWPHYL_DATAFLOWNUMERIC_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/Signed.h>
#include <Eigen/Core>
#include <string>
#include <typeinfo>
#include <utility>

namespace bpp {
namespace DF {

	// Typedefs TODO use wrapper types that can be forward declared
	using VectorDouble = Eigen::VectorXd;
	using MatrixDouble = Eigen::MatrixXd;

	// Debug info
	std::string debugInfoFor (const double & d);
	std::string debugInfoFor (const MatrixDouble & m);

	// Dimensions
	struct MatrixDimension {
		SizeType rows;
		SizeType cols;
		constexpr MatrixDimension (SizeType rows_, SizeType cols_) noexcept
		    : rows (rows_), cols (cols_) {}
		std::string toString () const;
	};
	constexpr bool operator== (const MatrixDimension & lhs, const MatrixDimension & rhs) noexcept {
		return lhs.rows == rhs.rows && lhs.cols == rhs.cols;
	}
	constexpr bool operator!= (const MatrixDimension & lhs, const MatrixDimension & rhs) noexcept {
		return !(lhs == rhs);
	}
	inline MatrixDimension dimensions (const MatrixDouble & m) noexcept {
		return {m.rows (), m.cols ()};
	}
	inline MatrixDimension dimensions (const Value<MatrixDouble> & node) noexcept {
		return dimensions (node.accessValue ());
	}

	/* Double nodes.
	 */
	struct ConstantDouble : public Value<double> {
		ConstantDouble (double d);
		void compute () override final;
		std::string debugInfo () const override final;
		bool isConstant () const override final;
		NodeRef derive (const Node &) override final;
		static std::shared_ptr<ConstantDouble> zero;
		static std::shared_ptr<ConstantDouble> one;
		static std::shared_ptr<ConstantDouble> create (double d);
	};

	struct AddDouble : public Value<double> {
		using Dependencies = ReductionOfValue<double>;
		AddDouble (NodeRefVec && deps);
		void compute () override final;
		std::string debugInfo () const override final;
		NodeRef derive (const Node & node) override final;
		static ValueRef<double> create (NodeRefVec && deps);
	};

	struct MulDouble : public Value<double> {
		using Dependencies = ReductionOfValue<double>;
		MulDouble (NodeRefVec && deps);
		void compute () override final;
		std::string debugInfo () const override final;
		NodeRef derive (const Node & node) override final;
		static ValueRef<double> create (NodeRefVec && deps);
	};

	/* Matrix nodes.
	 */
	struct ConstantMatrixDouble : public Value<MatrixDouble> {
		template <typename Derived>
		ConstantMatrixDouble (const Eigen::EigenBase<Derived> & expr)
		    : Value<MatrixDouble> (noDependency, expr) {
			this->makeValid ();
		}
		void compute () override final;
		std::string debugInfo () const override final;
		NodeRef derive (const Node & node) override final;
		template <typename Derived>
		static std::shared_ptr<ConstantMatrixDouble> create (const Eigen::EigenBase<Derived> & expr) {
			return std::make_shared<ConstantMatrixDouble> (expr);
		}
		static std::shared_ptr<ConstantMatrixDouble> createZero (MatrixDimension dim);
	};

	struct AddMatrixDouble : public Value<MatrixDouble> {
		using Dependencies = ReductionOfValue<MatrixDouble>;
		AddMatrixDouble (NodeRefVec && deps, MatrixDimension dim);
		void compute () override final;
		std::string debugInfo () const override final;
		NodeRef derive (const Node & node) override final;
		static ValueRef<MatrixDouble> create (NodeRefVec && deps, MatrixDimension dim);
	};

	struct MulMatrixDouble : public Value<MatrixDouble> {
		using Dependencies = FunctionOfValues<MatrixDouble, MatrixDouble>;
		MulMatrixDouble (NodeRefVec && deps, MatrixDimension dim);
		void compute () override final;
		std::string debugInfo () const override final;
		NodeRef derive (const Node & node) override final;
		static ValueRef<MatrixDouble> create (NodeRefVec && deps, MatrixDimension dim);
		static ValueRef<MatrixDouble> create (ValueRef<MatrixDouble> lhs, ValueRef<MatrixDouble> rhs,
		                                      MatrixDimension dim);
	};

	struct CWiseMulMatrixDouble : public Value<MatrixDouble> {
		using Dependencies = ReductionOfValue<MatrixDouble>;
		CWiseMulMatrixDouble (NodeRefVec && deps, MatrixDimension dim);
		void compute () override final;
		std::string debugInfo () const override final;
		NodeRef derive (const Node & node) override final;
		static ValueRef<MatrixDouble> create (NodeRefVec && deps, MatrixDimension dim);
	};
} // namespace DF
} // namespace bpp
#endif // BPP_NEWPHYL_DATAFLOWNUMERIC_H
