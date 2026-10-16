#ifndef INTERPOLATION_NN_HPP_
#define INTERPOLATION_NN_HPP_

#include "interpolation_types.hpp"

namespace interpolation
{

	class InterpolationNN;
	using InterpolationNNPtr = std::shared_ptr< InterpolationNN >;
	using ConstInterpolationNNPtr = std::shared_ptr< const InterpolationNN >;

	InterpolationNNPtr CreateInterpolationNN(IBasisPtr basis);

	class InterpolationNN : public IInterpolation, public std::enable_shared_from_this<InterpolationNN>
	{
	public:

		static InterpolationNNPtr Create(IBasisPtr basis);

		virtual ~InterpolationNN() = default;

		InterpolationType GetType() const override;
		Matrix GetValue(const Vector& point) const override;
		Matrices GetValue(const Vectors& points) const override;
		NumberCoordinates GetNumberCoordinates() const override;

		void SetNodes(const Nodes& nodes) override;

	protected:
		InterpolationNN() = default;

		void SetBasis(IBasisPtr basis);
		
		InterpolationType type_{ InterpolationType::nn };
		TreePtr tree_{ nullptr };
		IBasisPtr basis_{ nullptr };
	};

} // namespace interpolation

#endif /* INTERPOLATION_NN_HPP_*/