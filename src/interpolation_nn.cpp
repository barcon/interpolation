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
		basis_ = basis::CreateBasisCartesian();
	}
	InterpolationNNPtr InterpolationNN::Create()
	{
		class MakeSharedEnabler : public InterpolationNN
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

		return res;

	}
	InterpolationNNPtr InterpolationNN::GetPtr()
	{
		return this->shared_from_this();
	}
	ConstInterpolationNNPtr InterpolationNN::GetPtr() const
	{
		return const_cast<InterpolationNN*>(this)->GetPtr();
	}
	Type InterpolationNN::GetType() const
	{
		return type_;
	}
	Matrix InterpolationNN::GetValue(Scalar x, Scalar y, Scalar z) const
	{
		auto output = node::CreateNode();

		output->SetPoint(x, y, z);

		GetValue(output);

		return output->GetValue();
	}
	void InterpolationNN::GetValue(INodePtr output) const
	{
		INodePtr found;

		if (tree_ == nullptr)
		{
			logger::Error(headerInterpolation, "Search tree not initialized");
			return;
		}

		found = tree_->Search(output);

		if (found == nullptr)
		{
			logger::Error(headerInterpolation, "Min. number of nodes could not be found for the interpolation.");
			return;
		}

		output->SetValue(found->GetValue());
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