#ifndef INTERPOLATION_NN_HPP_
#define INTERPOLATION_NN_HPP_

#include "interpolation_types.hpp"

namespace interpolation
{

	class InterpolationNN;
	using InterpolationNNPtr = std::shared_ptr< InterpolationNN >;
	using ConstInterpolationNNPtr = std::shared_ptr< const InterpolationNN >;

	InterpolationNNPtr CreateInterpolationNN();

	class InterpolationNN : public IInterpolation, public std::enable_shared_from_this<InterpolationNN>
	{
	public:

		static InterpolationNNPtr Create();

		virtual ~InterpolationNN() = default;

		Type GetType() const override;
		Matrix GetValue(const Vector& point) const override;

		void SetNodes(const Nodes& nodes) override;
		void SetBasis(IBasisPtr basis) override;

	protected:
		InterpolationNN();

		Type type_{ interpolation_nn };
		TreePtr tree_{ nullptr };
		IBasisPtr basis_{ nullptr };
	};

} // namespace interpolation

#endif /* INTERPOLATION_NN_HPP_*/