#ifndef INTERPOLATION_IDW_HPP_
#define INTERPOLATION_IDW_HPP_

#include "interpolation_types.hpp"

namespace interpolation
{
	class InterpolationIDW;
	using InterpolationIDWPtr = std::shared_ptr< InterpolationIDW >;
	using ConstInterpolationIDWPtr = std::shared_ptr< const InterpolationIDW >;

	InterpolationIDWPtr CreateInterpolationIDW();

	class InterpolationIDW : public IInterpolation, public std::enable_shared_from_this<InterpolationIDW>
	{
	public:
		virtual ~InterpolationIDW() = default;

		static InterpolationIDWPtr Create();
		InterpolationIDWPtr GetPtr();
		ConstInterpolationIDWPtr GetPtr() const;

		Type GetType() const override;
		Matrix GetValue(Scalar x, Scalar y, Scalar z) const override;
		void GetValue(INodePtr output) const override;

		void SetNodes(const Nodes& nodes) override;
		void SetBasis(IBasisPtr basis) override;
		void SetShape(Scalar value);
		void SetNumberNodes(NumberNodes numberNodes);

	protected:
		InterpolationIDW();

		Type type_{ interpolation_idw };
		TreePtr tree_{ nullptr };
		IBasisPtr basis_{ nullptr };
		Scalar p_{ 2.0 };
		NumberNodes n_{ 3 };
	};

} // namespace interpolation

#endif /* INTERPOLATION_IDW_HPP_*/