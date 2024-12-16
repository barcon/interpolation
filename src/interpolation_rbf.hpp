#ifndef INTERPOLATION_RBF_HPP_
#define INTERPOLATION_RBF_HPP_

#include "interpolation_types.hpp"

namespace interpolation
{
	static const Type interpolation_function_TPS = 1;
	static const Type interpolation_function_MQB = 2;

	class InterpolationRBF;
	using InterpolationRBFPtr = std::shared_ptr< InterpolationRBF >;
	using ConstInterpolationRBFPtr = std::shared_ptr< const InterpolationRBF >;

	InterpolationRBFPtr CreateInterpolationRBF();

	class InterpolationRBF : public IInterpolation, public std::enable_shared_from_this<InterpolationRBF>
	{
	public:
		using FunctionScalar = Scalar(InterpolationRBF::*)(const Vector&, const Vector&) const;

		virtual ~InterpolationRBF() = default;

		static InterpolationRBFPtr Create();
		InterpolationRBFPtr GetPtr();
		ConstInterpolationRBFPtr GetPtr() const;
		
		Type GetType() const override;
		Matrix GetValue(Scalar x, Scalar y, Scalar z) const override;
		void GetValue(NodePtr output) const override;

		void SetNodes(const Nodes& nodes) override;
		void SetBasis(BasisPtr basis) override;
		void SetShape(Scalar value);
		void SetFunction(Type function);

	protected:
		InterpolationRBF();

		Type type_{ interpolation_rbf };
		BasisPtr basis_{ nullptr };
		Nodes nodes_;
		Vectors alpha_;
		Scalar a_{ 1.0e-3 };

		Dimension dim_{ 0 };
		Index nn_{ 0 };
		Index rows1_{ 0 };
		Index cols1_{ 0 };
		Index rows2_{ 0 };
		Index cols2_{ 0 };

		FunctionScalar function_{ &InterpolationRBF::FunctionTPS };

		Scalar FunctionTPS(const Vector& point1, const Vector& point2) const;
		Scalar FunctionMQB(const Vector& point1, const Vector& point2) const;
	};

} // namespace interpolation

#endif /* INTERPOLATION_RBF_HPP_*/