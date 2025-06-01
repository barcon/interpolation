#include "interpolation_fixed.hpp"

namespace interpolation
{
	InterpolationFixedPtr CreateInterpolationFixed()
	{
		auto res = InterpolationFixed::Create();

		return res;
	}
	InterpolationFixed::InterpolationFixed()
	{
		basis_ = basis::CreateBasisCartesian();
	}
	InterpolationFixedPtr InterpolationFixed::Create()
	{
		class MakeSharedEnabler : public InterpolationFixed
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

		return res;

	}
	InterpolationFixedPtr InterpolationFixed::GetPtr()
	{
		return this->shared_from_this();
	}
	ConstInterpolationFixedPtr InterpolationFixed::GetPtr() const
	{
		return const_cast<InterpolationFixed*>(this)->GetPtr();
	}
	Type InterpolationFixed::GetType() const
	{
		return type_;
	}
	Matrix InterpolationFixed::GetValue(Scalar x, Scalar y, Scalar z) const
	{
		auto output = nodes::CreateNode();

		output->SetPoint(x, y, z);

		GetValue(output);

		return output->GetValue();
	}
	void InterpolationFixed::GetValue(INodePtr output) const
	{
		output->SetValue(nodes_[n_]->GetValue());
	}
	void InterpolationFixed::SetNodes(const Nodes& nodes)
	{
		nodes_ = nodes;
	}
	void InterpolationFixed::SetNodeIndex(Index index)
	{
		n_ = index;
	}
	void InterpolationFixed::SetBasis(IBasisPtr basis)
	{
		basis_ = basis;
	}
} //namespace interpolation