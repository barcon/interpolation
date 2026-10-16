#ifndef INTERPOLATION_TYPES_HPP_
#define INTERPOLATION_TYPES_HPP_

#include "utils.hpp"
#include "logger.hpp"
#include "eilig.hpp"
#include "basis.hpp"
#include "node.hpp"
#include "kdtree.hpp"

#ifdef EILIG_ENABLE_OPENCL
#include "club.hpp"
#endif

#include <variant>
#include <unordered_map>
#include <memory>

namespace interpolation
{
	using namespace eilig;
	using namespace basis;
	using namespace node;
	using namespace kdtree;

	const String headerInterpolation = "INTERPOLATION";

	enum class InterpolationType
	{
		fixed = 0,
		idw = 1,
		nn = 2,
		rbf = 3
	};

	class IInterpolation;
	using IInterpolationPtr = std::shared_ptr< IInterpolation >;
	using ConstIInterpolationPtr = std::shared_ptr< const IInterpolation >;

	class IInterpolation
	{
	public:
		virtual ~IInterpolation() = default;

		virtual InterpolationType GetType() const = 0;
		virtual Matrix GetValue(const Vector& point) const = 0;
		virtual Matrices GetValue(const Vectors& points) const = 0;
		virtual NumberCoordinates GetNumberCoordinates() const = 0;

		virtual void SetNodes(const Nodes& nodes) = 0;

	protected:
		IInterpolation() = default;
	};

} // namespace interpolation

#endif /* INTERPOLATION_TYPES_HPP_*/