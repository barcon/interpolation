#include "interpolation_nn.hpp"

namespace interpolation
{
	InterpolationNNPtr CreateInterpolationNN()
	{
		auto res = InterpolationNN::Create();

		return res;
	}
	InterpolationNN::InterpolationNN()
	{
		basis_ = basis::CreateBasisCartesian(0);
	}
	InterpolationNNPtr InterpolationNN::Create()
	{
		class MakeSharedEnabler : public InterpolationNN
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

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

		found = tree_->Search(output);

		if (found == nullptr)
		{
			throw std::runtime_error("Min. number of nodes could not be found for the interpolation.");
		}

		return found->GetValue();
	}
	void InterpolationNN::SetNodes(const Nodes& nodes)
	{
		if (nodes.size() == 0)
		{
			logger::Error(headerInterpolation, "Incompatible size: number of nodes");
			return;
		}

		if (basis_ == nullptr)
		{
			logger::Error(headerInterpolation, "Basis pointer is nullptr. Set a basis first");
			return;
		}

		tree_ = kdtree::CreateTree();
		tree_->SetBasis(basis_);
		tree_->MakeTree(nodes);
	}
	void InterpolationNN::SetBasis(IBasisPtr basis)
	{
		basis_ = basis;
	}
} //namespace interpolation