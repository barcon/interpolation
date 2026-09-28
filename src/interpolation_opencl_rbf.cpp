#include "interpolation_opencl_rbf.hpp"
#include "interpolation_rbf.hpp"

#include <cmath>
#include <algorithm>

namespace interpolation
{
	namespace opencl
	{
		InterpolationRBFPtr CreateInterpolationRBF(IBasisPtr basis, Type function, Scalar shape)
		{
			auto res = InterpolationRBF::Create(basis, function, shape);

			return res;
		}
		InterpolationRBFPtr InterpolationRBF::Create(IBasisPtr basis, Type function, Scalar shape)
		{
			class MakeSharedEnabler : public InterpolationRBF
			{
			};

			auto res = std::make_shared<MakeSharedEnabler>();

			res->SetBasis(basis);
			res->SetFunction(function, shape);

			return res;
		}
		
		InterpolationRBF::InterpolationRBF()
		{
			kernel_ = KernelInterpolationRBF::Create();
		}
		Type InterpolationRBF::GetType() const
		{
			return type_;
		}
		Matrix InterpolationRBF::GetValue(const Vector& point) const
		{
			auto numberNodes = nodes_.size();
			auto numberDof = nodes_[0]->GetNumberDof();
			auto numberPoints = 1;
			auto numberRows = nodes_[0]->GetValue().GetRows();
			auto numberCols = nodes_[0]->GetValue().GetCols();

			Matrix res(numberRows, numberCols, eilig::matrix_zeros);
			
			eilig::opencl::Vector rbf(numberNodes);

			if (alphaGPU_ == nullptr)
			{
				throw std::runtime_error("Radial basis functions not initialized");
			}

			if (point.GetRows() != GetNumberCoordinates())
			{
				throw std::invalid_argument("Point dimensions do not match basis coordinates.");
			}

			auto pointsGPU = CreateBufferPoints({point});
			auto rbfGPU = (this->*function2_)(numberNodes, numberPoints, pointsGPU, nodesGPU_);
			auto resultGPU = DotProduct(numberNodes, numberPoints, numberDof, rbfGPU, alphaGPU_);

			return CreateResult(numberPoints, numberDof, resultGPU)[0];
		}
		Matrices InterpolationRBF::GetValue(const Vectors& points) const
		{
			auto numberNodes = nodes_.size();
			auto numberDof = nodes_[0]->GetNumberDof();
			auto numberPoints = points.size();
			auto numberRows = nodes_[0]->GetValue().GetRows();
			auto numberCols = nodes_[0]->GetValue().GetCols();

			Matrices res;

			if (alphaGPU_ == nullptr)
			{
				throw std::runtime_error("Radial basis functions not initialized");
			}

			auto pointsGPU = CreateBufferPoints(points);
			auto rbfGPU = (this->*function2_)(numberNodes, numberPoints, pointsGPU, nodesGPU_);
			auto resultGPU = DotProduct(numberNodes, numberPoints, numberDof, rbfGPU, alphaGPU_);			

			return CreateResult(numberPoints, numberDof, resultGPU);
		}

		NumberCoordinates InterpolationRBF::GetNumberCoordinates() const
		{
			return basis_->GetNumberCoordinates();
		}

		void InterpolationRBF::SetNodes(const Nodes& nodes)
		{
			Vectors nodalPoints;
			Status status{ eilig::EILIG_NOT_CONVERGED };

			if (nodes.size() == 0)
			{
				throw std::invalid_argument("Nodes cannot be empty.");
			}

			auto numberDof = nodes[0]->GetNumberDof();
			auto numberNodes = nodes.size();

			for (auto& node : nodes)
			{
				if (node == nullptr)
				{
					throw std::invalid_argument("Nodes cannot contain null pointers.");
				}

				if (node->GetNumberCoordinates() != GetNumberCoordinates())
				{
					throw std::invalid_argument("Node coordinates do not match basis coordinates.");
				}

				if (node->GetNumberDof() != numberDof)
				{
					throw std::invalid_argument("All nodes must have the same number of degrees of freedom.");
				}
			}

			std::transform(nodes.begin(), nodes.end(), std::back_inserter(nodalPoints), [](INodePtr node) -> Vector {
				return node->GetPoint();
				});

			nodes_ = nodes;
			nodesGPU_ = CreateBufferPoints(nodalPoints);
			alphaGPU_ = club::CreateBuffer(GetContext(), sizeof(Scalar) * numberNodes * numberDof);

			logger::Info(headerInterpolation, "Assembly of linear system started.");

			eilig::opencl::Matrix A = (this->*function1_)(numberNodes, nodesGPU_);
			eilig::opencl::Vector y(numberNodes);			

			logger::Info(headerInterpolation, "Assembly of linear system completed.");
			logger::Info(headerInterpolation, utils::string::Format("Solving InterpolationRBF linear system: {} equations", numberNodes));
		
			auto aux = Vector(numberNodes);
			for (DofIndex dof = 0; dof < numberDof; ++dof)
			{
				for (Index i = 0; i < numberNodes; ++i)
				{
					aux(i) = nodes_[i]->GetValue(dof);
				}

				auto d = eilig::opencl::Vector(aux);

				status = eilig::IterativeBiCGStab(A, y, d, callbackIterative);
				if (status != eilig::EILIG_SUCCESS)
				{
					logger::Error(headerInterpolation, utils::string::Format("Interpolation failed to initialize: {}", eilig::messages.at(eilig::EILIG_NOT_CONVERGED)));
					return;
				}

				alphaGPU_->Write(GetContext()->GetQueues()[deviceIndex_], sizeof(Scalar) * dof * numberNodes, sizeof(Scalar) * numberNodes, y.GetDataGPU()->Get(), CL_TRUE);
			}
		}
		void InterpolationRBF::SetBasis(IBasisPtr basis)
		{
			if (basis == nullptr)
			{
				throw std::invalid_argument("Basis cannot be null.");
			}

			if (basis->GetNumberCoordinates() != 3)
			{
				throw std::invalid_argument("Basis must have 3 coordinates.");
			}

			basis_ = basis;
		}
		void InterpolationRBF::SetFunction(Type function, Scalar shape)
		{
			switch (function)
			{
			case interpolation_function_INQ:
				function1_ = &InterpolationRBF::FunctionINQ1;
				function2_ = &InterpolationRBF::FunctionINQ2;
				break;
			case interpolation_function_GAU:
				function1_ = &InterpolationRBF::FunctionGAU1;
				function2_ = &InterpolationRBF::FunctionGAU2;
				break;
			case interpolation_function_MQB:
				function1_ = &InterpolationRBF::FunctionMQB1;
				function2_ = &InterpolationRBF::FunctionMQB2;
				break;
			case interpolation_function_TPS:
			default:
				function1_ = &InterpolationRBF::FunctionTPS1;
				function2_ = &InterpolationRBF::FunctionTPS2;
				break;
			}

			shape_ = shape;
		}

		BufferPtr InterpolationRBF::DotProduct(NumberNodes numberNodes, Number numberPoints, NumberDof numberDof, BufferPtr rbfGPU, BufferPtr alphaGPU) const
		{
			club::Events events(1);
			club::Error error;
			Index globalSize[1];

			auto res = club::CreateBuffer(GetContext(), sizeof(Scalar) * numberPoints * numberDof);

			const auto& dimension = kernel_->kDotProduct_->GetDim();
			const auto& localSize = GetContext()->GetLocalSize(deviceIndex_, dimension);

			globalSize[0] = GlobalSize(numberPoints, localSize[0]);

			kernel_->kDotProduct_->SetArg(0, sizeof(Index), &numberNodes);
			kernel_->kDotProduct_->SetArg(1, sizeof(Index), &numberPoints);
			kernel_->kDotProduct_->SetArg(2, sizeof(Index), &numberDof);
			kernel_->kDotProduct_->SetArg(3, sizeof(cl_mem), &rbfGPU->Get());
			kernel_->kDotProduct_->SetArg(4, sizeof(cl_mem), &alphaGPU->Get());
			kernel_->kDotProduct_->SetArg(5, sizeof(cl_mem), &res->Get());
			
			error = clEnqueueNDRangeKernel(GetContext()->GetQueues()[deviceIndex_],
				kernel_->kDotProduct_->GetKernel(),
				kernel_->kDotProduct_->GetDim(), NULL, globalSize,
				&localSize[0], 0, NULL, &events[0]);

			if (error != CL_SUCCESS)
			{
				logger::Error(headerInterpolation, utils::string::Format("Enqueueing kernel: {}", club::messages.at(error)));
			}

			clWaitForEvents(static_cast<cl_uint>(events.size()), &events[0]);

			return res;
		}
		BufferPtr InterpolationRBF::CreateBufferPoints(const Vectors& points) const
		{
			auto numberPoints = points.size();
			auto numberCoordinates = GetNumberCoordinates();

			Scalars buffer(numberPoints * numberCoordinates);

			for (Index i = 0; i < numberPoints; ++i)
			{
				if (points[i].GetRows() != numberCoordinates)
				{
					throw std::invalid_argument("Point dimensions do not match basis coordinates.");
				}

				for (Index j = 0; j < numberCoordinates; ++j)
				{
					buffer[i * numberCoordinates + j] = points[i](j);
				}
			}

			auto bufferGPU = club::CreateBuffer(GetContext(), sizeof(Scalar) * buffer.size());
			bufferGPU->Write(GetContext()->GetQueues()[deviceIndex_], 0, sizeof(Scalar) * buffer.size(), buffer.data(), CL_TRUE);

			return bufferGPU;
		}
		Matrices InterpolationRBF::CreateResult(Number numberPoints, NumberDof numberDof, BufferPtr resultGPU) const
		{
			Matrices res;
			Scalars buffer(numberPoints * numberDof);

			resultGPU->Read(GetContext()->GetQueues()[deviceIndex_], 0, sizeof(Scalar) * buffer.size(), &buffer[0], CL_TRUE);

			//for()

			return res;
		}
		
		eilig::opencl::Matrix InterpolationRBF::FunctionTPS1(NumberNodes numberNodes, BufferPtr nodes) const
		{
			club::Events events(1);
			club::Error error;
			Index globalSize[1];

			auto numberRows = numberNodes;
			auto numberCols = numberNodes;

			eilig::opencl::Matrix res(numberRows, numberCols);

			const auto& dimension = kernel_->kFunctionTPS1_->GetDim();
			const auto& localSize = GetContext()->GetLocalSize(deviceIndex_, dimension);

			globalSize[0] = GlobalSize(numberRows, localSize[0]);

			kernel_->kFunctionTPS1_->SetArg(0, sizeof(Index), &numberRows);
			kernel_->kFunctionTPS1_->SetArg(1, sizeof(Index), &numberCols);
			kernel_->kFunctionTPS1_->SetArg(2, sizeof(cl_mem), &nodes->Get());
			kernel_->kFunctionTPS1_->SetArg(3, sizeof(cl_mem), &res.GetDataGPU()->Get());

			error = clEnqueueNDRangeKernel(GetContext()->GetQueues()[deviceIndex_],
				kernel_->kFunctionTPS1_->GetKernel(),
				kernel_->kFunctionTPS1_->GetDim(), NULL, globalSize,
				&localSize[0], 0, NULL, &events[0]);

			if (error != CL_SUCCESS)
			{
				logger::Error(headerInterpolation, utils::string::Format("Enqueueing kernel: {}", club::messages.at(error)));
			}
			
			clWaitForEvents(static_cast<cl_uint>(events.size()), &events[0]);

			return res;
		}
		eilig::opencl::Matrix InterpolationRBF::FunctionMQB1(NumberNodes numberNodes, BufferPtr nodes) const
		{
			club::Events events(1);
			club::Error error;
			Index globalSize[1];

			auto numberRows = numberNodes;
			auto numberCols = numberNodes;

			eilig::opencl::Matrix res(numberRows, numberCols);

			const auto& dimension = kernel_->kFunctionMQB1_->GetDim();
			const auto& localSize = GetContext()->GetLocalSize(deviceIndex_, dimension);

			globalSize[0] = GlobalSize(numberRows, localSize[0]);

			kernel_->kFunctionMQB1_->SetArg(0, sizeof(Index), &numberRows);
			kernel_->kFunctionMQB1_->SetArg(1, sizeof(Index), &numberCols);
			kernel_->kFunctionMQB1_->SetArg(2, sizeof(cl_mem), &nodes->Get());
			kernel_->kFunctionMQB1_->SetArg(3, sizeof(cl_mem), &res.GetDataGPU()->Get());
			kernel_->kFunctionMQB1_->SetArg(4, sizeof(Scalar), &shape_);
			
			error = clEnqueueNDRangeKernel(GetContext()->GetQueues()[deviceIndex_],
				kernel_->kFunctionMQB1_->GetKernel(),
				kernel_->kFunctionMQB1_->GetDim(), NULL, globalSize,
				&localSize[0], 0, NULL, &events[0]);

			if (error != CL_SUCCESS)
			{
				logger::Error(headerInterpolation, utils::string::Format("Enqueueing kernel: {}", club::messages.at(error)));
			}
			
			clWaitForEvents(static_cast<cl_uint>(events.size()), &events[0]);

			return res;
		}
		eilig::opencl::Matrix InterpolationRBF::FunctionGAU1(NumberNodes numberNodes, BufferPtr nodes) const
		{
			club::Events events(1);
			club::Error error;
			Index globalSize[1];

			auto numberRows = numberNodes;
			auto numberCols = numberNodes;

			eilig::opencl::Matrix res(numberRows, numberCols);

			const auto& dimension = kernel_->kFunctionGAU1_->GetDim();
			const auto& localSize = GetContext()->GetLocalSize(deviceIndex_, dimension);

			globalSize[0] = GlobalSize(numberRows, localSize[0]);

			kernel_->kFunctionGAU1_->SetArg(0, sizeof(Index), &numberRows);
			kernel_->kFunctionGAU1_->SetArg(1, sizeof(Index), &numberCols);
			kernel_->kFunctionGAU1_->SetArg(2, sizeof(cl_mem), &nodes->Get());
			kernel_->kFunctionGAU1_->SetArg(3, sizeof(cl_mem), &res.GetDataGPU()->Get());
			kernel_->kFunctionGAU1_->SetArg(4, sizeof(Scalar), &shape_);
			
			error = clEnqueueNDRangeKernel(GetContext()->GetQueues()[deviceIndex_],
				kernel_->kFunctionGAU1_->GetKernel(),
				kernel_->kFunctionGAU1_->GetDim(), NULL, globalSize,
				&localSize[0], 0, NULL, &events[0]);

			if (error != CL_SUCCESS)
			{
				logger::Error(headerInterpolation, utils::string::Format("Enqueueing kernel: {}", club::messages.at(error)));
			}

			clWaitForEvents(static_cast<cl_uint>(events.size()), &events[0]);

			return res;
		}
		eilig::opencl::Matrix InterpolationRBF::FunctionINQ1(NumberNodes numberNodes, BufferPtr nodes) const
		{
			club::Events events(1);
			club::Error error;
			Index globalSize[1];

			auto numberRows = numberNodes;
			auto numberCols = numberNodes;

			eilig::opencl::Matrix res(numberRows, numberCols);

			const auto& dimension = kernel_->kFunctionINQ1_->GetDim();
			const auto& localSize = GetContext()->GetLocalSize(deviceIndex_, dimension);

			globalSize[0] = GlobalSize(numberRows, localSize[0]);

			kernel_->kFunctionINQ1_->SetArg(0, sizeof(Index), &numberRows);
			kernel_->kFunctionINQ1_->SetArg(1, sizeof(Index), &numberCols);
			kernel_->kFunctionINQ1_->SetArg(2, sizeof(cl_mem), &nodes->Get());
			kernel_->kFunctionINQ1_->SetArg(3, sizeof(cl_mem), &res.GetDataGPU()->Get());
			kernel_->kFunctionINQ1_->SetArg(4, sizeof(Scalar), &shape_);

			error = clEnqueueNDRangeKernel(GetContext()->GetQueues()[deviceIndex_],
				kernel_->kFunctionINQ1_->GetKernel(),
				kernel_->kFunctionINQ1_->GetDim(), NULL, globalSize,
				&localSize[0], 0, NULL, &events[0]);

			if (error != CL_SUCCESS)
			{
				logger::Error(headerInterpolation, utils::string::Format("Enqueueing kernel: {}", club::messages.at(error)));
			}

			clWaitForEvents(static_cast<cl_uint>(events.size()), &events[0]);

			return res;
		}

		BufferPtr InterpolationRBF::FunctionTPS2(NumberNodes numberNodes, Number numberPoints, BufferPtr points, BufferPtr nodes) const
		{
			club::Events events(1);
			club::Error error;
			Index globalSize[1];

			auto res = club::CreateBuffer(GetContext(), sizeof(Scalar) * numberPoints * numberNodes);

			const auto& dimension = kernel_->kFunctionTPS2_->GetDim();
			const auto& localSize = GetContext()->GetLocalSize(deviceIndex_, dimension);

			globalSize[0] = GlobalSize(numberPoints, localSize[0]);

			kernel_->kFunctionTPS2_->SetArg(0, sizeof(Index), &numberPoints);
			kernel_->kFunctionTPS2_->SetArg(1, sizeof(Index), &numberNodes);
			kernel_->kFunctionTPS2_->SetArg(2, sizeof(cl_mem), &points->Get());
			kernel_->kFunctionTPS2_->SetArg(3, sizeof(cl_mem), &nodes->Get());
			kernel_->kFunctionTPS2_->SetArg(4, sizeof(cl_mem), &res->Get());

			error = clEnqueueNDRangeKernel(GetContext()->GetQueues()[deviceIndex_],
				kernel_->kFunctionTPS2_->GetKernel(),
				kernel_->kFunctionTPS2_->GetDim(), NULL, globalSize,
				&localSize[0], 0, NULL, &events[0]);

			if (error != CL_SUCCESS)
			{
				logger::Error(headerInterpolation, utils::string::Format("Enqueueing kernel: {}", club::messages.at(error)));
			}

			clWaitForEvents(static_cast<cl_uint>(events.size()), &events[0]);

			return res;
		}
		BufferPtr InterpolationRBF::FunctionMQB2(NumberNodes numberNodes, Number numberPoints, BufferPtr points, BufferPtr nodes) const
		{
			club::Events events(1);
			club::Error error;
			Index globalSize[1];

			auto res = club::CreateBuffer(GetContext(), sizeof(Scalar) * numberPoints * numberNodes);

			const auto& dimension = kernel_->kFunctionMQB2_->GetDim();
			const auto& localSize = GetContext()->GetLocalSize(deviceIndex_, dimension);

			globalSize[0] = GlobalSize(numberPoints, localSize[0]);

			kernel_->kFunctionMQB2_->SetArg(0, sizeof(Index), &numberPoints);
			kernel_->kFunctionMQB2_->SetArg(1, sizeof(Index), &numberNodes);
			kernel_->kFunctionMQB2_->SetArg(2, sizeof(cl_mem), &points->Get());
			kernel_->kFunctionMQB2_->SetArg(3, sizeof(cl_mem), &nodes->Get());
			kernel_->kFunctionMQB2_->SetArg(4, sizeof(cl_mem), &res->Get());
			kernel_->kFunctionMQB2_->SetArg(5, sizeof(Scalar), &shape_);
			
			error = clEnqueueNDRangeKernel(GetContext()->GetQueues()[deviceIndex_],
				kernel_->kFunctionMQB2_->GetKernel(),
				kernel_->kFunctionMQB2_->GetDim(), NULL, globalSize,
				&localSize[0], 0, NULL, &events[0]);

			if (error != CL_SUCCESS)
			{
				logger::Error(headerInterpolation, utils::string::Format("Enqueueing kernel: {}", club::messages.at(error)));
			}

			clWaitForEvents(static_cast<cl_uint>(events.size()), &events[0]);

			return res;
		}
		BufferPtr InterpolationRBF::FunctionGAU2(NumberNodes numberNodes, Number numberPoints, BufferPtr points, BufferPtr nodes) const
		{
			club::Events events(1);
			club::Error error;
			Index globalSize[1];

			auto res = club::CreateBuffer(GetContext(), sizeof(Scalar) * numberPoints * numberNodes);

			const auto& dimension = kernel_->kFunctionGAU2_->GetDim();
			const auto& localSize = GetContext()->GetLocalSize(deviceIndex_, dimension);

			globalSize[0] = GlobalSize(numberPoints, localSize[0]);

			kernel_->kFunctionGAU2_->SetArg(0, sizeof(Index), &numberPoints);
			kernel_->kFunctionGAU2_->SetArg(1, sizeof(Index), &numberNodes);
			kernel_->kFunctionGAU2_->SetArg(2, sizeof(cl_mem), &points->Get());
			kernel_->kFunctionGAU2_->SetArg(3, sizeof(cl_mem), &nodes->Get());
			kernel_->kFunctionGAU2_->SetArg(4, sizeof(cl_mem), &res->Get());
			kernel_->kFunctionGAU2_->SetArg(5, sizeof(Scalar), &shape_);

			error = clEnqueueNDRangeKernel(GetContext()->GetQueues()[deviceIndex_],
				kernel_->kFunctionGAU2_->GetKernel(),
				kernel_->kFunctionGAU2_->GetDim(), NULL, globalSize,
				&localSize[0], 0, NULL, &events[0]);

			if (error != CL_SUCCESS)
			{
				logger::Error(headerInterpolation, utils::string::Format("Enqueueing kernel: {}", club::messages.at(error)));
			}

			clWaitForEvents(static_cast<cl_uint>(events.size()), &events[0]);

			return res;
		}
		BufferPtr InterpolationRBF::FunctionINQ2(NumberNodes numberNodes, Number numberPoints, BufferPtr points, BufferPtr nodes) const
		{
			club::Events events(1);
			club::Error error;
			Index globalSize[1];

			auto res = club::CreateBuffer(GetContext(), sizeof(Scalar) * numberPoints * numberNodes);

			const auto& dimension = kernel_->kFunctionINQ2_->GetDim();
			const auto& localSize = GetContext()->GetLocalSize(deviceIndex_, dimension);

			globalSize[0] = GlobalSize(numberPoints, localSize[0]);

			kernel_->kFunctionINQ2_->SetArg(0, sizeof(Index), &numberPoints);
			kernel_->kFunctionINQ2_->SetArg(1, sizeof(Index), &numberNodes);
			kernel_->kFunctionINQ2_->SetArg(2, sizeof(cl_mem), &points->Get());
			kernel_->kFunctionINQ2_->SetArg(3, sizeof(cl_mem), &nodes->Get());
			kernel_->kFunctionINQ2_->SetArg(4, sizeof(cl_mem), &res->Get());
			kernel_->kFunctionINQ2_->SetArg(5, sizeof(Scalar), &shape_);

			error = clEnqueueNDRangeKernel(GetContext()->GetQueues()[deviceIndex_],
				kernel_->kFunctionINQ2_->GetKernel(),
				kernel_->kFunctionINQ2_->GetDim(), NULL, globalSize,
				&localSize[0], 0, NULL, &events[0]);

			if (error != CL_SUCCESS)
			{
				logger::Error(headerInterpolation, utils::string::Format("Enqueueing kernel: {}", club::messages.at(error)));
			}

			clWaitForEvents(static_cast<cl_uint>(events.size()), &events[0]);

			return res;
		}
	}
} //namespace interpolation