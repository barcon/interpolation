#include "interpolation_rbf.hpp"
#include <cmath>

namespace interpolation
{
	Status CallbackIterative(Index iteration, Scalar residual)
	{
		Scalar tolerance = 1.0e-5;
		if (std::isnan(residual))
		{
			return eilig::EILIG_NOT_CONVERGED;
		}

		if (residual < tolerance)
		{
			return eilig::EILIG_SUCCESS;
		}

		return eilig::EILIG_CONTINUE;
	}

	InterpolationRBFPtr CreateInterpolationRBF()
	{
		auto res = InterpolationRBF::Create();

		return res;
	}
	InterpolationRBF::InterpolationRBF()
	{
		basis_ = basis::CreateBasisCartesian();
	}
	InterpolationRBFPtr InterpolationRBF::Create()
	{
		class MakeSharedEnabler : public InterpolationRBF
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

		return res;

	}
	InterpolationRBFPtr InterpolationRBF::GetPtr()
	{
		return this->shared_from_this();
	}
	ConstInterpolationRBFPtr InterpolationRBF::GetPtr() const
	{
		return const_cast<InterpolationRBF*>(this)->GetPtr();
	}
	Type InterpolationRBF::GetType() const
	{
		return type_;
	}
	Matrix InterpolationRBF::GetValue(Scalar x, Scalar y, Scalar z) const
	{
		auto output = node::CreateNode();

		output->SetPoint(x, y, z);

		GetValue(output);

		return output->GetValue();
	}
	void InterpolationRBF::GetValue(INodePtr output) const
	{
		Matrix vec(rows1_, 1, eilig::matrix_zeros);
		Matrix res(rows2_, cols2_, eilig::matrix_zeros);
		Scalar aux{ 0.0 };
		unsigned int counter{ 0 };

		if (alpha_.size() == 0)
		{
			logger::Error(headerInterpolation, "Radial basis functions not initialized");
			return;
		}

		auto& point1 = output->GetPoint();

		for (std::size_t i = 0; i < nn_; ++i)
		{
			auto& point2 = nodes_[i]->GetPoint();

			vec(i, 0) = (this->*function_)(point1, point2);
		}

		vec(nn_, 0) = 1.0;
		for (std::size_t j = 0; j < dim_; ++j)
		{
			vec(nn_ + j + 1, 0) = point1(j);
		}

		for (std::size_t u = 0; u < rows2_; ++u)
		{
			for (std::size_t v = 0; v < cols2_; ++v)
			{
				for (std::size_t i = 0; i < rows1_; ++i)
				{
					res(u, v) += vec(i, 0) * alpha_[counter](i);
				}
				++counter;
			}
		}

		output->SetValue(res);
	}
	void InterpolationRBF::SetNodes(const Nodes& nodes)
	{
		if (nodes.size() == 0)
		{
			logger::Error(headerInterpolation, "Incompatible number of nodes");
			return;
		}

		if (basis_ == nullptr)
		{
			logger::Error(headerInterpolation, "Basis pointer is nullptr. Set a basis first");
			return;

		}

		logger::Info(headerInterpolation, "Solving InterpolationRBF linear system");

		nodes_ = nodes;
		dim_ = basis_->GetDim();
		nn_ = nodes_.size();
		rows1_ = dim_ + nn_ + 1;
		cols1_ = dim_ + nn_ + 1;
		rows2_ = nodes_[0]->GetValue().GetRows();
		cols2_ = nodes_[0]->GetValue().GetCols();

		Matrix A(rows1_, cols1_, eilig::matrix_zeros);
		Vector d(rows1_, 0.0);
		Vector y(rows1_, 0.0);
		Status status{ eilig::EILIG_SUCCESS };

		for (std::size_t i = 0; i < nn_; ++i)
		{
			A(i, nn_) = 1.0;
			A(nn_, i) = 1.0;

			auto& point1 = nodes_[i]->GetPoint();

			for (std::size_t j = 0; j < nn_; ++j)
			{
				auto& point2 = nodes_[j]->GetPoint();

				A(i, j) = (this->*function_)(point1, point2);
			}

			for (std::size_t j = 0; j < dim_; ++j)
			{
				A(i, nn_ + j + 1) = point1(j);
				A(nn_ + j + 1, i) = point1(j);
			}
		}

		for (std::size_t u = 0; u < rows2_; ++u)
		{
			for (std::size_t v = 0; v < cols2_; ++v)
			{
				for (std::size_t i = 0; i < nn_; ++i)
				{
					d(i) = nodes_[i]->GetValue()(u, v);
				}

				status = eilig::IterativeBiCGStab(A, y, d);
				if (status != eilig::EILIG_SUCCESS)
				{
					logger::Error(headerInterpolation, "Interpolation failed to initialize : " + eilig::messages.at(eilig::EILIG_NOT_CONVERGED));
					return;
				}

				alpha_.push_back(y);
			}
		}
	}
	void InterpolationRBF::SetBasis(IBasisPtr basis)
	{
		basis_ = basis;
	}
	void InterpolationRBF::SetShape(Scalar shape)
	{
		a_ = shape;
	}
	void InterpolationRBF::SetFunction(Type function)
	{
		switch (function)
		{
		case interpolation_function_MQB:
			function_ = &InterpolationRBF::FunctionMQB;
			break;
		case interpolation_function_TPS:
		default:
			function_ = &InterpolationRBF::FunctionTPS;
			break;
		}
	}
	Scalar InterpolationRBF::FunctionTPS(const Vector& point1, const Vector& point2) const
	{
		Scalar res{ 0.0 };
		Scalar x = basis_->Distance(point1, point2);

		if (utils::math::IsAlmostEqual(x, 0.0, 2))
		{
			res = 0.0;
		}
		else
		{
			res = x * x * std::log(x);
		}

		return res;
	}
	Scalar InterpolationRBF::FunctionMQB(const Vector& point1, const Vector& point2) const
	{
		Scalar res{ 0.0 };
		Scalar x = basis_->Distance(point1, point2);

		res = std::sqrt(a_ * a_ + x * x);

		return res;
	}
} //namespace interpolation