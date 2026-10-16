#include "interpolation_rbf.hpp"
#include <cmath>

namespace interpolation
{
	InterpolationRBFPtr CreateInterpolationRBF(IBasisPtr basis, FunctionType functionType, Scalar shape)
	{
		auto res = InterpolationRBF::Create(basis, functionType, shape);

		return res;
	}
	InterpolationRBFPtr InterpolationRBF::Create(IBasisPtr basis, FunctionType functionType, Scalar shape)
	{
		class MakeSharedEnabler : public InterpolationRBF
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

		res->SetBasis(basis);
		res->SetFunction(functionType, shape);

		return res;
	}

	InterpolationType InterpolationRBF::GetType() const
	{
		return type_;
	}
	Matrix InterpolationRBF::GetValue(const Vector& point) const
	{
		auto numberNodes = nodes_.size();
		auto numberDof = nodes_[0]->GetNumberDof();

		Vector rbf(numberNodes, 0.0);
		Matrix res(nodes_[0]->GetValue().GetRows(), nodes_[0]->GetValue().GetCols(), MatrixType::zeros);

		if (point.GetRows() != GetNumberCoordinates())
		{
			throw std::invalid_argument("Point dimensions do not match basis coordinates.");
		}

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
			res(dof) = eilig::Dot(rbf, alpha_[dof]);
		}

		return res;
	}
	Matrices InterpolationRBF::GetValue(const Vectors& points) const
	{
		auto numberNodes = nodes_.size();
		auto numberDof = nodes_[0]->GetNumberDof();

		Matrices res;
		Matrix aux(nodes_[0]->GetValue().GetRows(), nodes_[0]->GetValue().GetCols(), MatrixType::zeros);
		Vector rbf(numberNodes, 0.0);

		if (alpha_.size() == 0)
		{
			throw std::runtime_error("Radial basis functions not initialized");
		}

		for (const auto& point : points)
		{
			if (point.GetRows() != GetNumberCoordinates())
			{
				throw std::invalid_argument("Point dimensions do not match basis coordinates.");
			}

			for (Index i = 0; i < numberNodes; ++i)
			{
				rbf(i) = (this->*function_)(point, nodes_[i]->GetPoint());
			}

			for (DofIndex dof = 0; dof < numberDof; ++dof)
			{
				aux(dof) = eilig::Dot(rbf, alpha_[dof]);
			}

			res.push_back(aux);
		}

		return res;
	}
	NumberCoordinates InterpolationRBF::GetNumberCoordinates() const
	{
		return basis_->GetNumberCoordinates();
	}

	void InterpolationRBF::SetNodes(const Nodes& nodes)
	{
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

		nodes_ = nodes;
		alpha_.clear();

		logger::Info(headerInterpolation, "Assembly of linear system started.");

		Matrix A(numberNodes, numberNodes);
		Vector y(numberNodes);
		Vector d(numberNodes);;

		for (Index i = 0; i < numberNodes; ++i)
		{
			const auto& point1 = nodes_[i]->GetPoint();

			for (Index j = i; j < numberNodes; ++j)
			{
				const auto& point2 = nodes_[j]->GetPoint();

				A(i, j) = (this->*function_)(point2, point1);
				A(j, i) = A(i, j);
			}
		}
		logger::Info(headerInterpolation, "Assembly of linear system completed.");
		logger::Info(headerInterpolation, utils::string::Format("Solving InterpolationRBF linear system: {} equations", numberNodes));

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
	void InterpolationRBF::SetFunction(FunctionType functionType, Scalar shape)
	{
		switch (functionType)
		{
		case FunctionType::INQ:
			function_ = &InterpolationRBF::FunctionINQ;
			break;
		case FunctionType::GAU:
			function_ = &InterpolationRBF::FunctionGAU;
			break;
		case FunctionType::MQB:
			function_ = &InterpolationRBF::FunctionMQB;
			break;
		case FunctionType::TPS:
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