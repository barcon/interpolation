#ifndef INTERPOLATION_TYPES_HPP_
#define INTERPOLATION_TYPES_HPP_

#include "utils.hpp"
#include "logger.hpp"
#include "eilig.hpp"
#include "basis.hpp"
#include "nodes.hpp"
#include "kdtree.hpp"

#include <variant>
#include <unordered_map>
#include <memory>

namespace interpolation
{
	using Scalar = eilig::Scalar;
	using Vector = eilig::Vector;
	using Matrix = eilig::Matrix;
	using Status = eilig::Status;

	using Scalars = eilig::Scalars;
	using Vectors = eilig::Vectors;

	using String = utils::String;

	using IBasisPtr = basis::IBasisPtr;
	using INodePtr = nodes::INodePtr;

	using Nodes = nodes::Nodes;
	using TreePtr = kdtree::TreePtr;
	
	const String headerInterpolation = "INTERPOLATION";

	using Dimension = basis::Dimension;
	using Index = std::size_t;
	using NumberNodes = std::size_t;

	using Type = std::size_t;
	static const Type interpolation_fixed	= 1;
	static const Type interpolation_idw		= 2;
	static const Type interpolation_nn		= 3;
	static const Type interpolation_rbf		= 4;

	class IInterpolation;
	using IInterpolationPtr = std::shared_ptr< IInterpolation >;
	using ConstIInterpolationPtr = std::shared_ptr< const IInterpolation >;

	class IInterpolation
	{
	public:
		virtual ~IInterpolation() = default;

		virtual Type GetType() const = 0;
		virtual Matrix GetValue(Scalar x, Scalar y, Scalar z) const = 0;
		virtual void GetValue(INodePtr output) const = 0;

		virtual void SetNodes(const Nodes& nodes) = 0;
		virtual void SetBasis(IBasisPtr basis) = 0;

	protected:
		IInterpolation() = default;
	};

} // namespace interpolation

#endif /* INTERPOLATION_TYPES_HPP_*/