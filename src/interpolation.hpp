#ifndef INTERPOLATION_HPP_
#define INTERPOLATION_HPP_

#include "interpolation_types.hpp"
#include "interpolation_fixed.hpp"
#include "interpolation_idw.hpp"
#include "interpolation_rbf.hpp"
#include "interpolation_nn.hpp"

#ifdef EILIG_ENABLE_OPENCL
#include "interpolation_opencl_kernel.hpp"
#include "interpolation_opencl_rbf.hpp"
#endif // EILIG_ENABLE_OPENCL


#endif /* INTERPOLATION_HPP_*/