#ifndef INTERPOLATION_FIXED_HPP_
#define INTERPOLATION_FIXED_HPP_

#include "interpolation_types.hpp"

namespace interpolation
{
	class InterpolationFixed;
	using InterpolationFixedPtr = std::shared_ptr< InterpolationFixed >;
	using ConstInterpolationFixedPtr = std::shared_ptr< const InterpolationFixed >;

	InterpolationFixedPtr CreateInterpolationFixed();

	class InterpolationFixed : public IInterpolation, public std::enable_shared_from_this<InterpolationFixed>
	{
	public:

		static InterpolationFixedPtr Create();
		InterpolationFixedPtr GetPtr();
		ConstInterpolationFixedPtr GetPtr() const;

		virtual ~InterpolationFixed() = default;

		Type GetType() const override;
		Matrix GetValue(Scalar x, Scalar y, Scalar z) const override;
		void GetValue(NodePtr output) const override;

		void SetNodes(const Nodes& nodes) override;
		void SetBasis(BasisPtr basis) override;
		void SetNodeIndex(Index index);

	protected:
		InterpolationFixed();

		Type type_{ interpolation_fixed };
		Nodes nodes_;
		BasisPtr basis_{ nullptr };

		Index n_{ 0 };
	};

} // namespace interpolation

#endif /* INTERPOLATION_FIXED_HPP_*/