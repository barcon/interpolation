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
		basis_ = basis::CreateBasisCartesian(0);
	}
	InterpolationFixedPtr InterpolationFixed::Create()
	{
		class MakeSharedEnabler : public InterpolationFixed
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

		return res;

	}
	Type InterpolationFixed::GetType() const
	{
		return type_;
	}
	Matrix InterpolationFixed::GetValue(const Vector& point) const
	{
		return nodes_[n_]->GetValue();
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