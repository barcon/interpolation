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
		InterpolationNNPtr GetPtr();
		ConstInterpolationNNPtr GetPtr() const;

		virtual ~InterpolationNN() = default;

		Type GetType() const override;
		Matrix GetValue(Scalar x, Scalar y, Scalar z) const override;
		void GetValue(NodePtr output) const override;

		void SetNodes(const Nodes& nodes) override;
		void SetBasis(BasisPtr basis) override;

	protected:
		InterpolationNN();

		Type type_{ interpolation_nn };
		TreePtr tree_{ nullptr };
		BasisPtr basis_{ nullptr };
	};

} // namespace interpolation

#endif /* INTERPOLATION_NN_HPP_*/