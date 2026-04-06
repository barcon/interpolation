#include "interpolation_fixed.hpp"

namespace interpolation
{
	InterpolationFixedPtr CreateInterpolationFixed(IBasisPtr basis, const Nodes& nodes, Index nodeIndex)
	{
		auto res = InterpolationFixed::Create(basis);

		res->SetNodes(nodes);
		res->SetNodeIndex(nodeIndex);

		return res;
	}
	InterpolationFixedPtr InterpolationFixed::Create(IBasisPtr basis)
	{
		class MakeSharedEnabler : public InterpolationFixed
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

		res->SetBasis(basis);

		return res;
	}
	Type InterpolationFixed::GetType() const
	{
		return type_;
	}
	Matrix InterpolationFixed::GetValue(const Vector& point) const
	{
		if(point.GetRows() != GetNumberCoordinates())
		{
			throw std::invalid_argument("Point dimensions do not match basis coordinates.");
		}

		return nodes_[nodeIndex_]->GetValue();
	}
	NumberCoordinates InterpolationFixed::GetNumberCoordinates() const
	{
		return basis_->GetNumberCoordinates();
	}
	void InterpolationFixed::SetNodes(const Nodes& nodes)
	{
		if (nodes.size() == 0)
		{
			throw std::invalid_argument("Nodes cannot be empty.");
		}

		auto numberDof = nodes[0]->GetNumberDof();

		for(auto& node: nodes)
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

		SetNodeIndex(0);
	}
	void InterpolationFixed::SetNodeIndex(Index nodeIndex)
	{
		if(nodeIndex >= nodes_.size())
		{
			throw std::out_of_range("Node index is out of range.");
		}

		nodeIndex_ = nodeIndex;
	}
	void InterpolationFixed::SetBasis(IBasisPtr basis)
	{
		if (basis == nullptr)
		{
			throw std::invalid_argument("Basis cannot be null.");
		}

		basis_ = basis;
	}
} //namespace interpolation