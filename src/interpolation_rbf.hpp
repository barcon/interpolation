#ifndef INTERPOLATION_RBF_HPP_
#define INTERPOLATION_RBF_HPP_

#include "interpolation_types.hpp"

namespace interpolation
{
	static const Type interpolation_function_TPS = 1;
	static const Type interpolation_function_MQB = 2;
	static const Type interpolation_function_GAU = 3;
	static const Type interpolation_function_INQ = 4;

	class InterpolationRBF;
	using InterpolationRBFPtr = std::shared_ptr< InterpolationRBF >;
	using ConstInterpolationRBFPtr = std::shared_ptr< const InterpolationRBF >;

	InterpolationRBFPtr CreateInterpolationRBF();

	class InterpolationRBF : public IInterpolation, public std::enable_shared_from_this<InterpolationRBF>
	{
	public:
		using Function_S_VV = Scalar(InterpolationRBF::*)(const Vector&, const Vector&) const;
		using Function_S_V = Scalar(InterpolationRBF::*)(const Vector&) const;
		using Polynomials = std::vector<Function_S_V>;

		virtual ~InterpolationRBF() = default;

		static InterpolationRBFPtr Create();
		
		Type GetType() const override;
		Matrix GetValue(const Vector& point) const override;

		void SetNodes(const Nodes& nodes) override;
		void SetBasis(IBasisPtr basis) override;
		void SetFunction(Type function, Scalar shape);

	protected:
		InterpolationRBF();

		Type type_{ interpolation_rbf };
		IBasisPtr basis_{ nullptr };
		Nodes nodes_;
		Vectors alpha_;
		Scalar shape_{ 0.1 };

		NumberDof numberDof_{ 0 };
		NumberNodes numberNodes_{ 0 };

		Polynomials polynomials_;
		Function_S_VV function_{ &InterpolationRBF::FunctionTPS };

		Scalar FunctionTPS(const Vector& point1, const Vector& point2) const;
		Scalar FunctionMQB(const Vector& point1, const Vector& point2) const;
		Scalar FunctionGAU(const Vector& point1, const Vector& point2) const;
		Scalar FunctionINQ(const Vector& point1, const Vector& point2) const;
		
		Scalar FunctionP0C(const Vector& point) const;
		Scalar FunctionP1X(const Vector& point) const;
		Scalar FunctionP1Y(const Vector& point) const;
		Scalar FunctionP1Z(const Vector& point) const;
	};

} // namespace interpolation

#endif /* INTERPOLATION_RBF_HPP_*/