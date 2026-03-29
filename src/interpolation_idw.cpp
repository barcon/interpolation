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
		basis_ = basis::CreateBasisCartesian(0);
	}
	InterpolationIDWPtr InterpolationIDW::Create()
	{
		class MakeSharedEnabler : public InterpolationIDW
		{
		};

		auto res = std::make_shared<MakeSharedEnabler>();

		return res;

	}
	Type InterpolationIDW::GetType() const
	{
		return type_;
	}
	Matrix InterpolationIDW::GetValue(const Vector& point) const
	{
		Scalar sum = 0.;
		Scalars distance(n_);
		Scalars weight(n_);
		Nodes found;
		Matrix res;
		INodePtr output = node::CreateNode(0, point);

		if (tree_ == nullptr)
		{
			throw std::runtime_error("Search tree not initialized");
		}

		found = tree_->SearchNumber(output, n_);

		if (found.size() != n_)
		{
			throw std::runtime_error("Min. number of nodes could not be found for the interpolation.");
		}

		for (Index i = 0; i < n_; ++i)
		{
			distance[i] = basis_->Distance(found[i]->GetPoint(), output->GetPoint());

			if (utils::math::IsAlmostEqual(distance[i], 0., 2))
			{
				res = found[i]->GetValue();
				return res;
			}

			weight[i] = std::pow(distance[i], -p_);
			sum += weight[i];
		}

		res = weight[0] * found[0]->GetValue();
		for (Index i = 1; i < n_; ++i)
		{
			res = res + weight[i] * found[i]->GetValue();
		}

		sum = 1.0 / sum;
		res = res * sum;

		return res;
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