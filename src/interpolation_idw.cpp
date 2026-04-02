#include "interpolation_idw.hpp"

namespace interpolation
{
	InterpolationIDWPtr CreateInterpolationIDW(IBasisPtr basis, const Nodes& nodes)
	{
		auto res = InterpolationIDW::Create(basis);

		res->SetNodes(nodes);

		return res;
	}
	InterpolationIDWPtr InterpolationIDW::Create(IBasisPtr basis)
	{
		class MakeSharedEnabler : public InterpolationIDW
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

		res->SetBasis(basis);

		return res;
	}
	Type InterpolationIDW::GetType() const
	{
		return type_;
	}
	Matrix InterpolationIDW::GetValue(const Vector& point) const
	{
		Scalar sum = 0.;
		Scalars distance(numberNodes_);
		Scalars weight(numberNodes_);
		Nodes found;
		Matrix res;
		INodePtr output = node::CreateNode(0, point);

		if (tree_ == nullptr)
		{
			throw std::runtime_error("Search tree not initialized");
		}

		if (point.GetRows() != GetNumberCoordinates())
		{
			throw std::invalid_argument("Point dimensions do not match basis coordinates.");
		}

		found = tree_->SearchNumber(output, numberNodes_);

		if (found.size() != numberNodes_)
		{
			throw std::runtime_error("Min. number of nodes could not be found for the interpolation.");
		}

		for (Index i = 0; i < numberNodes_; ++i)
		{
			distance[i] = basis_->Distance(found[i]->GetPoint(), output->GetPoint());

			if (utils::math::IsAlmostEqual(distance[i], 0., 2))
			{
				res = found[i]->GetValue();
				return res;
			}

			weight[i] = std::pow(distance[i], -shape_);
			sum += weight[i];
		}

		res = weight[0] * found[0]->GetValue();
		for (Index i = 1; i < numberNodes_; ++i)
		{
			res = res + weight[i] * found[i]->GetValue();
		}

		return (1.0 / sum) * res;
	}
	NumberCoordinates InterpolationIDW::GetNumberCoordinates() const
	{
		return basis_->GetNumberCoordinates();
	}
	void InterpolationIDW::SetNodes(const Nodes & nodes)
	{
		if (nodes.size() == 0)
		{
			throw std::invalid_argument("Nodes cannot be empty.");
		}

		if (nodes.size() < numberNodes_)
		{
			throw std::invalid_argument("Number of nodes must be at least equal to the specified number of nodes for interpolation: " + std::to_string(numberNodes_));
		}

		auto numberDof = nodes[0]->GetNumberDof();

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

		tree_ = kdtree::CreateTree();
		tree_->SetBasis(basis_);
		tree_->MakeTree(nodes);
	}
	void InterpolationIDW::SetBasis(IBasisPtr basis)
	{
		if (basis == nullptr)
		{
			throw std::invalid_argument("Basis cannot be null.");
		}

		basis_ = basis;
	}
	void InterpolationIDW::SetShape(Scalar shape)
	{
		shape_ = shape;
	}
} //namespace interpolation