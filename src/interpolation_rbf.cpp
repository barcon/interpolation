#include "interpolation_rbf.hpp"
#include <cmath>

namespace interpolation
{
	CallbackIterative callbackIterative = [](Index iteration, Scalar residual) -> long long int
		{
			Scalar tolerance{ 1e-2 };

			if (std::isnan(residual))
			{
				return eilig::EILIG_NOT_CONVERGED;
			}

			if (residual < tolerance)
			{
				return eilig::EILIG_SUCCESS;
			}

			logger::Info(headerInterpolation, utils::string::Format("Iteration: {}, Residual: {}", iteration, residual));

			return eilig::EILIG_CONTINUE;
		};

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
	Type InterpolationRBF::GetType() const
	{
		return type_;
	}
	Matrix InterpolationRBF::GetValue(const Vector& point) const
	{
		auto numberNodes = nodes_.size();
		auto numberDof = nodes_[0]->GetNumberDof();

		Vector rbf(numberNodes, 0.0);
		Matrix output(nodes_[0]->GetValue().GetRows(), nodes_[0]->GetValue().GetCols(), eilig::matrix_zeros);

		if (alpha_.size() == 0)
		{
			throw std::runtime_error("Radial basis functions not initialized");
		}

		for (Index i = 0; i < numberNodes; ++i)
		{
			rbf(i) = (this->*function_)(point, nodes_[i]->GetPoint());
		}

		for (DofIndex dof = 0; dof < numberDof; ++dof)
		{
			output(dof) = eilig::Dot(rbf, alpha_[dof]);
		}

		return output;
	}
	NumberCoordinates InterpolationRBF::GetNumberCoordinates() const
	{
		return basis_->GetNumberCoordinates();
	}
	void InterpolationRBF::SetNodes(const Nodes& nodes)
	{
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

		nodes_ = nodes;
		alpha_.clear();

		logger::Info(headerInterpolation, utils::string::Format("Solving InterpolationRBF linear system: {} equations", numberNodes));

		Matrix A(numberNodes, numberNodes, eilig::matrix_zeros);
		Vector d(numberNodes, 0.0);
		Vector y(numberNodes, 0.0);
		Status status{ eilig::EILIG_NOT_CONVERGED };

		for (Index i = 0; i < numberNodes; ++i)
		{
			const auto& point1 = nodes_[i]->GetPoint();

			for (Index j = 0; j < numberNodes; ++j)
			{
				const auto& point2 = nodes_[j]->GetPoint();

				A(i, j) = (this->*function_)(point2, point1);

			}
		}
		logger::Info(headerInterpolation, "Assembly linear system completed.");

		for (DofIndex dof = 0; dof < numberDof; ++dof)
		{
			for (Index i = 0; i < numberNodes; ++i)
			{
				d(i) = nodes_[i]->GetValue(dof);
			}

			status = eilig::IterativeBiCGStab(A, y, d, callbackIterative);
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
		if (basis == nullptr)
		{
			throw std::invalid_argument("Basis cannot be null.");
		}

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
} //namespace interpolation

/*
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


for (Index k = 0; k < order; ++k)
{

	A(i, k + numberNodes) = (this->*polynomials_[k])(point2);
	A(k + numberNodes, i) = (this->*polynomials_[k])(point2);
}


switch (res->GetNumberCoordinates())
{
case 1:
	res->polynomials_.push_back(&InterpolationRBF::FunctionP0C);
	res->polynomials_.push_back(&InterpolationRBF::FunctionP1X);
	break;
case 2:
	res->polynomials_.push_back(&InterpolationRBF::FunctionP0C);
	res->polynomials_.push_back(&InterpolationRBF::FunctionP1X);
	res->polynomials_.push_back(&InterpolationRBF::FunctionP1Y);
	break;
case 3:
	res->polynomials_.push_back(&InterpolationRBF::FunctionP0C);
	res->polynomials_.push_back(&InterpolationRBF::FunctionP1X);
	res->polynomials_.push_back(&InterpolationRBF::FunctionP1Y);
	res->polynomials_.push_back(&InterpolationRBF::FunctionP1Z);
	break;
default:
	throw std::invalid_argument("Unsupported number of coordinates for RBF interpolation.");
}*/