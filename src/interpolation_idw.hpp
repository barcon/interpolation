#ifndef INTERPOLATION_IDW_HPP_
#define INTERPOLATION_IDW_HPP_

#include "interpolation_types.hpp"

namespace interpolation
{
	class InterpolationIDW;
	using InterpolationIDWPtr = std::shared_ptr< InterpolationIDW >;

	InterpolationIDWPtr CreateInterpolationIDW(IBasisPtr basis);

	class InterpolationIDW : public IInterpolation, public std::enable_shared_from_this<InterpolationIDW>
	{
	public:
		virtual ~InterpolationIDW() = default;

		static InterpolationIDWPtr Create(IBasisPtr basis);

		InterpolationType GetType() const override;
		Matrix GetValue(const Vector& point) const override;
		Matrices GetValue(const Vectors& points) const override;
		NumberCoordinates GetNumberCoordinates() const override;

		void SetNodes(const Nodes& nodes) override;
		void SetShape(Scalar value);

	protected:
		InterpolationIDW() = default;

		void SetBasis(IBasisPtr basis);

		InterpolationType type_{ InterpolationType::idw };
		TreePtr tree_{ nullptr };
		IBasisPtr basis_{ nullptr };
		Scalar shape_{ 2.0 };
		NumberNodes numberNodes_{ 3 };
		const NumberCoordinates numberCoordinates_{ 3 };
	};

} // namespace interpolation

#endif /* INTERPOLATION_IDW_HPP_*/