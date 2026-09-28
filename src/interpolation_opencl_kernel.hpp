#ifndef INTERPOLATION_OPENCL_KERNEL_HPP_
#define INTERPOLATION_OPENCL_KERNEL_HPP_

#include "interpolation_types.hpp"

namespace interpolation
{
    namespace opencl
    {
        using BufferPtr = club::BufferPtr;
        using PlatformPtr = club::PlatformPtr;
        using ContextPtr = club::ContextPtr;
        using ProgramPtr = club::ProgramPtr;
        using DeviceIndex = club::DeviceIndex;

        class KernelInterpolationRBF;
        using KernelInterpolationRBFPtr = std::shared_ptr<KernelInterpolationRBF>;

        constexpr Index GlobalSize(Index num1, Index num2)
        {
            return num2 * (num1 / num2 + (num1 % num2 != 0));
        }

		void Initialize(const String& fileName, const club::PlatformIndex& platformIndex, const club::DeviceIndices& deviceIndices);
		
        ContextPtr GetContext();

        static const String kFunctionTPS1 = "FunctionTPS1";
        static const String kFunctionMQB1 = "FunctionMQB1";
        static const String kFunctionGAU1 = "FunctionGAU1";
        static const String kFunctionINQ1 = "FunctionINQ1";
              
        static const String kFunctionTPS2 = "FunctionTPS2";
        static const String kFunctionMQB2 = "FunctionMQB2";
        static const String kFunctionGAU2 = "FunctionGAU2";
        static const String kFunctionINQ2 = "FunctionINQ2";
        
        static const String kDotProduct = "DotProduct";
        
        KernelInterpolationRBFPtr CreateKernelInterpolationRBF();

        class KernelInterpolationRBF
        {
        public:
            virtual ~KernelInterpolationRBF() = default;

            static KernelInterpolationRBFPtr Create();

            club::KernelPtr kFunctionTPS1_{ nullptr };
            club::KernelPtr kFunctionMQB1_{ nullptr };
            club::KernelPtr kFunctionGAU1_{ nullptr };
            club::KernelPtr kFunctionINQ1_{ nullptr };
            
            club::KernelPtr kFunctionTPS2_{ nullptr };
            club::KernelPtr kFunctionMQB2_{ nullptr };
            club::KernelPtr kFunctionGAU2_{ nullptr };
            club::KernelPtr kFunctionINQ2_{ nullptr };
            
            club::KernelPtr kDotProduct_{ nullptr };

        protected:
            KernelInterpolationRBF();
            KernelInterpolationRBF(const KernelInterpolationRBF& copy) = delete;
            KernelInterpolationRBF(KernelInterpolationRBF&& move) = delete;

            KernelInterpolationRBF& operator=(const KernelInterpolationRBF& copy) = delete;
            KernelInterpolationRBF& operator=(KernelInterpolationRBF&& copy) = delete;
        };


    } // namespace opencl
} // namespace interpolation

#endif