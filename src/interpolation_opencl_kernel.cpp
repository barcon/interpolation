#include "interpolation_opencl_kernel.hpp"

namespace interpolation
{
    namespace opencl
    {
        static ProgramPtr program{ nullptr };
        static ContextPtr context{ nullptr };
        static PlatformPtr platform{ nullptr };

        void InitPlatform()
        {
            platform = club::CreatePlatform();

            if (!platform)
            {
				throw std::runtime_error("Platform could not be created");
            }
        }
        void InitContext(const club::PlatformIndex& platformIndex, const club::DeviceIndices& deviceIndices)
        {
            context = club::CreateContext(platform, platformIndex, deviceIndices);

            if (!context)
            {
                throw std::runtime_error(utils::string::Format("Context ({}) could not be created", platformIndex));
            }
        }
        void InitProgram(const String& fileName)
        {
            program = club::CreateProgramFromFile(context, fileName);

            if (!program)
            {
                throw std::runtime_error(utils::string::Format("Program could not be created {}", fileName));
            }
        }
        void Initialize(const String& fileName, const club::PlatformIndex& platformIndex, const club::DeviceIndices& deviceIndices)
        {
            InitPlatform();
            InitContext(platformIndex, deviceIndices);
            InitProgram(fileName);
        }

        ContextPtr GetContext()
        {
            return context;
        }

        KernelInterpolationRBFPtr CreateKernelInterpolationRBF()
        {
            auto res = KernelInterpolationRBF::Create();

            return res;
        }
        KernelInterpolationRBFPtr KernelInterpolationRBF::Create()
        {
            class MakeSharedEnabler : public KernelInterpolationRBF
            {
            };

            auto res = std::make_shared<MakeSharedEnabler>();

            return res;
        }
        KernelInterpolationRBF::KernelInterpolationRBF()
        {
            kFunctionTPS1_ = club::CreateKernel(program, kFunctionTPS1, 1);
            kFunctionMQB1_ = club::CreateKernel(program, kFunctionMQB1, 1);
            kFunctionGAU1_ = club::CreateKernel(program, kFunctionGAU1, 1);
            kFunctionINQ1_ = club::CreateKernel(program, kFunctionINQ1, 1);
            
            kFunctionTPS2_ = club::CreateKernel(program, kFunctionTPS2, 1);
            kFunctionMQB2_ = club::CreateKernel(program, kFunctionMQB2, 1);
            kFunctionGAU2_ = club::CreateKernel(program, kFunctionGAU2, 1);
            kFunctionINQ2_ = club::CreateKernel(program, kFunctionINQ2, 1);
            
            kDotProduct_ = club::CreateKernel(program, kDotProduct, 1);
        }
 
    }
} /* namespace eilig */