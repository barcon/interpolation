#include "interpolation_nn.hpp"

namespace interpolation
{
	InterpolationNNPtr CreateInterpolationNN(IBasisPtr basis, const Nodes& nodes)
	{
		auto res = InterpolationNN::Create(basis);

		res->SetNodes(nodes);

		return res;
	}
	InterpolationNNPtr InterpolationNN::Create(IBasisPtr basis)
	{
		class MakeSharedEnabler : public InterpolationNN
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

		res->SetBasis(basis);

		return res;
	}
	Type InterpolationNN::GetType() const
	{
		return type_;
	}
	Matrix InterpolationNN::GetValue(const Vector& point) const
	{
		INodePtr found;
		INodePtr output = node::CreateNode(0, point);

		if (tree_ == nullptr)
		{
			throw std::runtime_error("Search tree not initialized");
		}

		if (point.GetRows() != GetNumberCoordinates())
		{
			throw std::invalid_argument("Point dimensions do not match basis coordinates.");
		}

		found = tree_->Search(output);

		if (found == nullptr)
		{
			throw std::runtime_error("Min. number of nodes could not be found for the interpolation.");
		}

		return found->GetValue();
	}
	NumberCoordinates InterpolationNN::GetNumberCoordinates() const
	{
		return basis_->GetNumberCoordinates();
	}
	void InterpolationNN::SetNodes(const Nodes& nodes)
	{
		if (nodes.size() == 0)
		{
			throw std::invalid_argument("Nodes cannot be empty.");
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
	void InterpolationNN::SetBasis(IBasisPtr basis)
	{
		if (basis == nullptr)
		{
			throw std::invalid_argument("Basis cannot be null.");
		}

		basis_ = basis;
	}
} //namespace interpolation