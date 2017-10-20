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
#include <Bpp/NewPhyl/DataFlowTemplates.h> // TODO rm when constant is moved away
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

	// Debug info defined as overrides by types TODO move DF.h
	template <> std::string Value<double>::debugInfo () const;
	template <> std::string Value<VectorDouble>::debugInfo () const;
	template <> std::string Value<MatrixDouble>::debugInfo () const;

	// Additional Constant<T> specialisation (TODO move to DFT.h, with forward declaration)
	template <> NodeRef Constant<VectorDouble>::derive (const Node & node);
	template <> struct Builder<Constant<VectorDouble>> {
		template <typename EigenVector>
		static std::shared_ptr<Constant<VectorDouble>> make (const EigenVector & v) {
			return std::make_shared<Constant<VectorDouble>> (v);
		}
		static std::shared_ptr<Constant<VectorDouble>> makeZero (SizeType size);
	};

	struct MatrixDimension;
	template <> NodeRef Constant<MatrixDouble>::derive (const Node & node);
	template <> struct Builder<Constant<MatrixDouble>> {
		template <typename EigenMatrix>
		static std::shared_ptr<Constant<MatrixDouble>> make (const EigenMatrix & m) {
			return std::make_shared<Constant<MatrixDouble>> (m);
		}
		static std::shared_ptr<Constant<MatrixDouble>> makeZero (const MatrixDimension & dim);
		static std::shared_ptr<Constant<MatrixDouble>> makeOne (const MatrixDimension & dim);
	};

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

	inline SizeType dimensions (const VectorDouble & v) noexcept { return v.rows (); }
	inline MatrixDimension dimensions (const MatrixDouble & m) noexcept {
		return {m.rows (), m.cols ()};
	}
	template <typename T>
	auto dimensions (const Value<T> & node) noexcept -> decltype (dimensions (node.accessValue ())) {
		return dimensions (node.accessValue ());
	}

	/* Double nodes.
	 */
	struct AddDouble : public Value<double> {
		using Dependencies = ReductionOfValue<double>;
		AddDouble (NodeRefVec && deps);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
	};
	template <> struct Builder<AddDouble> { static ValueRef<double> make (NodeRefVec && deps); };

	struct MulDouble : public Value<double> {
		using Dependencies = ReductionOfValue<double>;
		MulDouble (NodeRefVec && deps);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
	};
	template <> struct Builder<MulDouble> { static ValueRef<double> make (NodeRefVec && deps); };

	struct ScalarProdDouble : public Value<double> {
		using Dependencies = FunctionOfValues<VectorDouble, VectorDouble>;
		ScalarProdDouble (NodeRefVec && deps);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
	};
	template <> struct Builder<ScalarProdDouble> {
		static ValueRef<double> make (NodeRefVec && deps);
	};

	/* Vector nodes.
	 */
	struct AddVectorDouble : public Value<VectorDouble> {
		using Dependencies = ReductionOfValue<VectorDouble>;
		AddVectorDouble (NodeRefVec && deps, SizeType size);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
	};
	template <> struct Builder<AddVectorDouble> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, SizeType type);
	};

	struct CWiseInverseVectorDouble : public Value<VectorDouble> {
		using Dependencies = FunctionOfValues<VectorDouble>;
		CWiseInverseVectorDouble (NodeRefVec && deps, SizeType size);
		void compute () override final;
		// TODO NodeRef derive (const Node & node) override final;
	};
	template <> struct Builder<CWiseInverseVectorDouble> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, SizeType type);
	};

	/* Matrix nodes.
	 */
	struct AddMatrixDouble : public Value<MatrixDouble> {
		using Dependencies = ReductionOfValue<MatrixDouble>;
		AddMatrixDouble (NodeRefVec && deps, MatrixDimension dim);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
	};
	template <> struct Builder<AddMatrixDouble> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, MatrixDimension dim);
	};

	struct MulMatrixDouble : public Value<MatrixDouble> {
		using Dependencies = FunctionOfValues<MatrixDouble, MatrixDouble>;
		MulMatrixDouble (NodeRefVec && deps, MatrixDimension dim);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
	};
	template <> struct Builder<MulMatrixDouble> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, MatrixDimension dim);
	};

	struct CWiseMulMatrixDouble : public Value<MatrixDouble> {
		using Dependencies = ReductionOfValue<MatrixDouble>;
		CWiseMulMatrixDouble (NodeRefVec && deps, MatrixDimension dim);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
	};
	template <> struct Builder<CWiseMulMatrixDouble> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, MatrixDimension dim);
	};

	/* Combinations
	 */
	struct MulTransposedMatrixVectorDouble : public Value<VectorDouble> {
		using Dependencies = FunctionOfValues<MatrixDouble, VectorDouble>;
		MulTransposedMatrixVectorDouble (NodeRefVec && deps, SizeType size);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
	};
	template <> struct Builder<MulTransposedMatrixVectorDouble> {
		static ValueRef<VectorDouble> make (NodeRefVec && deps, SizeType type);
	};

	struct MulScalarMatrixDouble : public Value<MatrixDouble> {
		using Dependencies = FunctionOfValues<double, MatrixDouble>;
		MulScalarMatrixDouble (NodeRefVec && deps, MatrixDimension dim);
		void compute () override final;
		NodeRef derive (const Node & node) override final;
	};
	template <> struct Builder<MulScalarMatrixDouble> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, MatrixDimension dim);
	};
} // namespace DF
} // namespace bpp
#endif // BPP_NEWPHYL_DATAFLOWNUMERIC_H
