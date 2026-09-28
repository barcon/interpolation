#ifndef INTERPOLATION_RBF_HPP_
#define INTERPOLATION_RBF_HPP_

#include "interpolation_types.hpp"

namespace interpolation
{
	static const Type interpolation_function_TPS = 1;
	static const Type interpolation_function_MQB = 2;
	static const Type interpolation_function_GAU = 3;
	static const Type interpolation_function_INQ = 4;

	inline eilig::CallbackIterative callbackIterative = [](Index iteration, Scalar residual) -> long long int
		{
			Scalar tolerance{ 1e-2 };

			logger::Info(headerInterpolation, utils::string::Format("Iteration: {}, Residual: {}", iteration, residual));

			if (std::isnan(residual))
			{
				return eilig::EILIG_NOT_CONVERGED;
			}

			if (residual < tolerance)
			{
				return eilig::EILIG_SUCCESS;
			}


			return eilig::EILIG_CONTINUE;
		};

	class InterpolationRBF;
	using InterpolationRBFPtr = std::shared_ptr< InterpolationRBF >;
	using ConstInterpolationRBFPtr = std::shared_ptr< const InterpolationRBF >;

	InterpolationRBFPtr CreateInterpolationRBF(IBasisPtr basis, Type function, Scalar shape);

	class InterpolationRBF : public IInterpolation, public std::enable_shared_from_this<InterpolationRBF>
	{
	public:
		using FunctionPtr = Scalar(InterpolationRBF::*)(const Vector&, const Vector&) const;

		virtual ~InterpolationRBF() = default;

		static InterpolationRBFPtr Create(IBasisPtr basis, Type function, Scalar shape);
		
		Type GetType() const override;
		Matrix GetValue(const Vector& point) const override;
		Matrices GetValue(const Vectors& points) const override;
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

		FunctionPtr function_{ &InterpolationRBF::FunctionTPS };

		Scalar FunctionTPS(const Vector& point1, const Vector& point2) const;
		Scalar FunctionMQB(const Vector& point1, const Vector& point2) const;
		Scalar FunctionGAU(const Vector& point1, const Vector& point2) const;
		Scalar FunctionINQ(const Vector& point1, const Vector& point2) const;
	};

} // namespace interpolation

#endif /* INTERPOLATION_RBF_HPP_*/