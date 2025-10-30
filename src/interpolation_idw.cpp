#include "interpolation_idw.hpp"

namespace interpolation
{
	InterpolationIDWPtr CreateInterpolationIDW()
	{
		auto res = InterpolationIDW::Create();

		return res;
	}
	InterpolationIDW::InterpolationIDW()
	{
		basis_ = basis::CreateBasisCartesian();
	}
	InterpolationIDWPtr InterpolationIDW::Create()
	{
		class MakeSharedEnabler : public InterpolationIDW
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

		return res;

	}
	InterpolationIDWPtr InterpolationIDW::GetPtr()
	{
		return this->shared_from_this();
	}
	ConstInterpolationIDWPtr InterpolationIDW::GetPtr() const
	{
		return const_cast<InterpolationIDW*>(this)->GetPtr();
	}
	Type InterpolationIDW::GetType() const
	{
		return type_;
	}
	Matrix InterpolationIDW::GetValue(Scalar x, Scalar y, Scalar z) const
	{
		auto output = node::CreateNode();

		output->SetPoint(x, y, z);

		GetValue(output);

		return output->GetValue();
	}
	void InterpolationIDW::GetValue(INodePtr output) const
	{
		Scalar sum = 0.;
		Scalars distance(n_);
		Scalars weight(n_);
		Nodes found;
		Matrix res;

		if (tree_ == nullptr)
		{
			logger::Error(headerInterpolation, "Search tree not initialized");
			return;
		}

		found = tree_->SearchNumber(output, n_);

		if (found.size() != n_)
		{
			logger::Error(headerInterpolation, "Min. number of nodes could not be found for the interpolation.");
			return;
		}

		for (unsigned int i = 0; i < n_; ++i)
		{
			distance[i] = basis_->Distance(found[i]->GetPoint(), output->GetPoint());

			if (utils::math::IsAlmostEqual(distance[i], 0., 2))
			{
				res = found[i]->GetValue();
				output->SetValue(res);
				return;
			}

			weight[i] = std::pow(distance[i], -p_);
			sum += weight[i];
		}

		res = weight[0] * found[0]->GetValue();
		for (unsigned int i = 1; i < n_; ++i)
		{
			res = res + weight[i] * found[i]->GetValue();
		}

		sum = 1.0 / sum;
		res = res * sum;
		output->SetValue(res);
	}
	void InterpolationIDW::SetNodes(const Nodes & nodes)
	{
		tree_ = kdtree::CreateTree();
		tree_->SetBasis(basis_);
		tree_->MakeTree(nodes);
	}
	void InterpolationIDW::SetBasis(IBasisPtr basis)
	{
		basis_ = basis;
	}
	void InterpolationIDW::SetShape(Scalar shape)
	{
		p_ = shape;
	}
	void InterpolationIDW::SetNumberNodes(NumberNodes numberNodes)
	{
		n_ = numberNodes;
	}
} //namespace interpolation