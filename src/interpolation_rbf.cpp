#include "interpolation_rbf.hpp"
#include <cmath>

namespace interpolation
{
	InterpolationRBFPtr CreateInterpolationRBF()
	{
		auto res = InterpolationRBF::Create();

		return res;
	}
	InterpolationRBF::InterpolationRBF()
	{
		basis_ = basis::CreateBasisCartesian(0);

		polynomials_.push_back(&InterpolationRBF::FunctionP0C);
		polynomials_.push_back(&InterpolationRBF::FunctionP1X);
		polynomials_.push_back(&InterpolationRBF::FunctionP1Y);
		polynomials_.push_back(&InterpolationRBF::FunctionP1Z);
	}
	InterpolationRBFPtr InterpolationRBF::Create()
	{
		class MakeSharedEnabler : public InterpolationRBF
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

		return res;

	}
	Type InterpolationRBF::GetType() const
	{
		return type_;
	}
	Matrix InterpolationRBF::GetValue(const Vector& point) const
	{
		Number order = polynomials_.size();
		Vector rbf(numberNodes_ + order, 0.0);
		Matrix output(nodes_[0]->GetValue().GetRows(), nodes_[0]->GetValue().GetCols(), eilig::matrix_zeros);

		if (alpha_.size() == 0)
		{
			throw std::runtime_error("Radial basis functions not initialized");
		}

		for (Index i = 0; i < numberNodes_; ++i)
		{
			rbf(i) = (this->*function_)(point, nodes_[i]->GetPoint());
		}

		for (Index k = 0; k < order; ++k)
		{
			rbf(numberNodes_ + k) = (this->*polynomials_[k])(point);
		}

		for (DofIndex dof = 0; dof < numberDof_; ++dof)
		{
			output(dof) = eilig::Dot(rbf, alpha_[dof]);
		}

		return output;
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

		alpha_.clear();
		nodes_ = nodes;
		numberNodes_ = nodes_.size();
		numberDof_ = nodes_[0]->GetNumberDof();

		Number order = polynomials_.size();
		Matrix A(numberNodes_ + order, numberNodes_ + order, eilig::matrix_zeros);
		Vector d(numberNodes_ + order, 0.0);
		Vector y(numberNodes_ + order, 0.0);
		Status status{ eilig::EILIG_SUCCESS };

		for (Index i = 0; i < numberNodes_; ++i)
		{
			const auto& point1 = nodes_[i]->GetPoint();

			for (Index j = 0; j < numberNodes_; ++j)
			{
				const auto& point2 = nodes_[j]->GetPoint();

				A(i, j) = (this->*function_)(point2, point1);

				for (Index k = 0; k < order; ++k)
				{

					A(i, k + numberNodes_) = (this->*polynomials_[k])(point2);
					A(k + numberNodes_, i) = (this->*polynomials_[k])(point2);
				}
			}
		}

		for (DofIndex dof = 0; dof < numberDof_; ++dof)
		{
			for (Index i = 0; i < numberNodes_; ++i)
			{
				d(i) = nodes_[i]->GetValue(dof);
			}

			status = eilig::IterativeBiCGStab(A, y, d);
			if (status != eilig::EILIG_SUCCESS)
			{
				logger::Error(headerInterpolation, utils::string::Format("Interpolation failed to initialize: {}", eilig::messages.at(eilig::EILIG_NOT_CONVERGED)));
				return;
			}

			alpha_.push_back(y);
		}
	}
	void InterpolationRBF::SetBasis(IBasisPtr basis)
	{
		basis_ = basis;
	}
	void InterpolationRBF::SetFunction(Type function, Scalar shape)
	{
		switch (function)
		{
		case interpolation_function_INQ:
			function_ = &InterpolationRBF::FunctionINQ;
			break;
		case interpolation_function_GAU:
			function_ = &InterpolationRBF::FunctionGAU;
			break;
		case interpolation_function_MQB:
			function_ = &InterpolationRBF::FunctionMQB;
			break;
		case interpolation_function_TPS:
		default:
			function_ = &InterpolationRBF::FunctionTPS;
			break;
		}

		shape_ = shape;
	}
	Scalar InterpolationRBF::FunctionTPS(const Vector& point1, const Vector& point2) const
	{
		Scalar r = basis_->Distance(point1, point2);

		return utils::math::IsAlmostEqual(r, 0.0, 2) ? 0.0 : r * r * std::log(r);
	}
	Scalar InterpolationRBF::FunctionMQB(const Vector& point1, const Vector& point2) const
	{
		Scalar r = basis_->Distance(point1, point2);

		return std::sqrt(shape_ * shape_ + r * r);
	}
	Scalar InterpolationRBF::FunctionGAU(const Vector& point1, const Vector& point2) const
	{
		Scalar r = basis_->Distance(point1, point2);

		return std::exp(-(shape_ * r) * (shape_ * r));

	}
	Scalar InterpolationRBF::FunctionINQ(const Vector& point1, const Vector& point2) const
	{
		Scalar r = basis_->Distance(point1, point2);

		return 1.0 / (r * r + shape_ * shape_);

	}
	
	Scalar InterpolationRBF::FunctionP0C(const Vector& point) const
	{
		return 1.0;
	}
	Scalar InterpolationRBF::FunctionP1X(const Vector& point) const
	{
		return point(0);
	}
	Scalar InterpolationRBF::FunctionP1Y(const Vector& point) const
	{
		return point(1);
	}
	Scalar InterpolationRBF::FunctionP1Z(const Vector& point) const
	{
		return point(2);
	}
} //namespace interpolation