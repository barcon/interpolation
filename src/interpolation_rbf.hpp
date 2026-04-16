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

	InterpolationRBFPtr CreateInterpolationRBF(IBasisPtr basis, Type function, Scalar shape);

	class InterpolationRBF : public IInterpolation, public std::enable_shared_from_this<InterpolationRBF>
	{
	public:
		using Function_S_VV = Scalar(InterpolationRBF::*)(const Vector&, const Vector&) const;
		using Function_S_V = Scalar(InterpolationRBF::*)(const Vector&) const;
		using Polynomials = std::vector<Function_S_V>;

		virtual ~InterpolationRBF() = default;

		static InterpolationRBFPtr Create(IBasisPtr basis, Type function, Scalar shape);
		
		Type GetType() const override;
		Matrix GetValue(const Vector& point) const override;
		NumberCoordinates GetNumberCoordinates() const override;

		void SetNodes(const Nodes& nodes) override;

	protected:
		InterpolationRBF() = default;
		
		void SetBasis(IBasisPtr basis);
		void SetFunction(Type function, Scalar shape);

		Type type_{ interpolation_rbf };
		IBasisPtr basis_{ nullptr };
		Nodes nodes_;
		Vectors alpha_;
		Scalar shape_{ 0.1 };

		Polynomials polynomials_;
		Function_S_VV function_{ &InterpolationRBF::FunctionTPS };

		Scalar FunctionTPS(const Vector& point1, const Vector& point2) const;
		Scalar FunctionMQB(const Vector& point1, const Vector& point2) const;
		Scalar FunctionGAU(const Vector& point1, const Vector& point2) const;
		Scalar FunctionINQ(const Vector& point1, const Vector& point2) const;
	};

} // namespace interpolation

#endif /* INTERPOLATION_RBF_HPP_*/