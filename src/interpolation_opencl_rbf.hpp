#ifndef INTERPOLATION_OPENCL_RBF_HPP_
#define INTERPOLATION_OPENCL_RBF_HPP_

#include "interpolation_rbf.hpp"
#include "interpolation_opencl_kernel.hpp"

namespace interpolation
{
	namespace opencl
	{
		class InterpolationRBF;
		using InterpolationRBFPtr = std::shared_ptr< InterpolationRBF >;
		using ConstInterpolationRBFPtr = std::shared_ptr< const InterpolationRBF >;

		InterpolationRBFPtr CreateInterpolationRBF(IBasisPtr basis, Type function, Scalar shape);

		class InterpolationRBF : public IInterpolation, public std::enable_shared_from_this<InterpolationRBF>
		{
		public:
			using FunctionPtr1 = eilig::opencl::Matrix (InterpolationRBF::*)(NumberNodes, BufferPtr) const;
			using FunctionPtr2 = BufferPtr (InterpolationRBF::*)(NumberNodes, Number, BufferPtr, BufferPtr) const;

			virtual ~InterpolationRBF() = default;

			static InterpolationRBFPtr Create(IBasisPtr basis, Type function, Scalar shape);

			Type GetType() const override;
			Matrix GetValue(const Vector& point) const override;
			Matrices GetValue(const Vectors& points) const override;
			NumberCoordinates GetNumberCoordinates() const override;

			void SetNodes(const Nodes& nodes) override;

		protected:
			InterpolationRBF();

			void SetBasis(IBasisPtr basis);
			void SetFunction(Type function, Scalar shape);			
			
			BufferPtr DotProduct(NumberNodes numberNodes, Number numberPoints, NumberDof numberDof, BufferPtr rbfGPU, BufferPtr alphaGPU) const;			
			BufferPtr CreateBufferPoints(const Vectors& points) const;
			Matrices CreateResult(Number numberPoints, NumberDof numberDof, BufferPtr resultGPU) const;
			
			Type type_{ interpolation_rbf };
			IBasisPtr basis_{ nullptr };
			Nodes nodes_;
			Scalar shape_{ 0.1 };
			
			FunctionPtr1 function1_{ &InterpolationRBF::FunctionTPS1 };
			FunctionPtr2 function2_{ &InterpolationRBF::FunctionTPS2 };
		
			BufferPtr alphaGPU_{ nullptr };
			BufferPtr nodesGPU_{ nullptr };
			DeviceIndex deviceIndex_{ 0 };
			KernelInterpolationRBFPtr kernel_{ nullptr };

			eilig::opencl::Matrix FunctionTPS1(NumberNodes numberNodes, BufferPtr nodes) const;
			eilig::opencl::Matrix FunctionMQB1(NumberNodes numberNodes, BufferPtr nodes) const;
			eilig::opencl::Matrix FunctionGAU1(NumberNodes numberNodes, BufferPtr nodes) const;
			eilig::opencl::Matrix FunctionINQ1(NumberNodes numberNodes, BufferPtr nodes) const;
			
			BufferPtr FunctionTPS2(NumberNodes numberNodes, Number numberPoints, BufferPtr points, BufferPtr nodes) const;
			BufferPtr FunctionMQB2(NumberNodes numberNodes, Number numberPoints, BufferPtr points, BufferPtr nodes) const;
			BufferPtr FunctionGAU2(NumberNodes numberNodes, Number numberPoints, BufferPtr points, BufferPtr nodes) const;
			BufferPtr FunctionINQ2(NumberNodes numberNodes, Number numberPoints, BufferPtr points, BufferPtr nodes) const;
		};
	}
} // namespace interpolation

#endif /* INTERPOLATION_OPENCL_RBF_HPP_*/