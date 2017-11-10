//
// File: LinearAlgebra.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-10-24
// Last modified: 2017-10-24
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

#ifndef BPP_NEWPHYL_LINEARALGEBRA_H
#define BPP_NEWPHYL_LINEARALGEBRA_H

#include <Bpp/NewPhyl/DataFlow.h>          // Value<T> override declarations
#include <Bpp/NewPhyl/DataFlowTemplates.h> // Constant<T> override declarations
#include <Bpp/NewPhyl/LinearAlgebraFwd.h>
#include <Eigen/Core>
#include <memory>

namespace bpp {
/** Define wrapper types to eigen matrix / vector.
 * No typedef, as typedef cannot be forward declared.
 */
class VectorDouble : public Eigen::VectorXd {
public:
	using Eigen::VectorXd::Base;
	using Eigen::VectorXd::Matrix; // Constructors
};
class MatrixDouble : public Eigen::MatrixXd {
public:
	using Eigen::MatrixXd::Base;
	using Eigen::MatrixXd::Matrix; // Constructors
};
} // namespace bpp

namespace Eigen {
namespace internal {
	// Explain to eigen that wrapper types are just wrappers.
	// Follows guide at https://eigen.tuxfamily.org/dox/TopicNewExpressionType.html

	// Traits
	template <> struct traits<bpp::VectorDouble> : traits<VectorXd> {};
	template <> struct traits<bpp::MatrixDouble> : traits<MatrixXd> {};

	// Evaluators
	template <> struct evaluator<bpp::VectorDouble> : evaluator<VectorXd> {};
	template <> struct evaluator<bpp::MatrixDouble> : evaluator<MatrixXd> {};
} // namespace internal
} // namespace Eigen

namespace bpp {
namespace DF {
	/* Declare overrides of Value<T>, Constant<T>.
	 *
	 * These overrides cannot be declared without a complete type.
	 * Thus no nice declaration with forward declared types in the headers of bpp::DF.
	 */

	// Value<T>.
	template <> std::string Value<VectorDouble>::debugInfo () const;
	template <> std::string Value<MatrixDouble>::debugInfo () const;

	// Constant<T>
	template <> NodeRef Constant<VectorDouble>::derive (const Node & node);
	template <> bool Constant<VectorDouble>::isDerivable (const Node & node);
	template <> struct Builder<Constant<VectorDouble>> {
		template <typename EigenVector>
		static std::shared_ptr<Constant<VectorDouble>> make (const EigenVector & v) {
			return std::make_shared<Constant<VectorDouble>> (v);
		}
		static std::shared_ptr<Constant<VectorDouble>> makeZero (const VectorDimension & dim);
	};

	template <> NodeRef Constant<MatrixDouble>::derive (const Node & node);
	template <> bool Constant<MatrixDouble>::isDerivable (const Node & node);
	template <> struct Builder<Constant<MatrixDouble>> {
		template <typename EigenMatrix>
		static std::shared_ptr<Constant<MatrixDouble>> make (const EigenMatrix & m) {
			return std::make_shared<Constant<MatrixDouble>> (m);
		}
		static std::shared_ptr<Constant<MatrixDouble>> makeZero (const MatrixDimension & dim);
		static std::shared_ptr<Constant<MatrixDouble>> makeOne (const MatrixDimension & dim);
	};

} // namespace DF
} // namespace bpp

#endif // BPP_NEWPHYL_LINEARALGEBRA_H
