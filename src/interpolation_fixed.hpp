#ifndef INTERPOLATION_FIXED_HPP_
#define INTERPOLATION_FIXED_HPP_

#include "interpolation_types.hpp"

namespace interpolation
{
	class InterpolationFixed;
	using InterpolationFixedPtr = std::shared_ptr< InterpolationFixed >;

	InterpolationFixedPtr CreateInterpolationFixed(IBasisPtr basis);

	class InterpolationFixed : public IInterpolation, public std::enable_shared_from_this<InterpolationFixed>
	{
	public:

		static InterpolationFixedPtr Create(IBasisPtr basis);

		virtual ~InterpolationFixed() = default;

		Type GetType() const override;
		Matrix GetValue(const Vector& point) const override;
		NumberCoordinates GetNumberCoordinates() const override;

		void SetNodes(const Nodes& nodes) override;
		void SetNodeIndex(Index nodeIndex);

	protected:
		InterpolationFixed() = default;

		void SetBasis(IBasisPtr basis);

		Type type_{ interpolation_fixed };
		Nodes nodes_;
		Index nodeIndex_{ 0 };
		IBasisPtr basis_{ nullptr };
	};

} // namespace interpolation

#endif /* INTERPOLATION_FIXED_HPP_*/