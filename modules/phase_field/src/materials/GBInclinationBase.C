#include "GBInclinationBase.h"

InputParameters
GBInclinationBase::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Base material to determine inclination dependent properties for AGG. "
                             "Builds up to the inclination and its derivatives wrt $\nabla\eta$.");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  MooseEnum gb_id_method("graintracker=0 ffc=1", "graintracker");
  params.addParam<MooseEnum>(
      "gb_id_method", gb_id_method, "Which GB OP identification/selection option.");
  params.addParam<UserObjectName>("grain_tracker",
                                  "The GrainTracker UserObject to get values from.");
  params.addParam<UserObjectName>("ffc", "The FFC UserObject to get values from.");
  MooseEnum angular_func("atan_2D=0 atan_3D=1", "atan_2D");
  params.addParam<MooseEnum>("angular_func",
                             angular_func,
                             "Which angular distance function to use. "
                             "atan_2D: oriented angle atan2(y,x) in [0,2pi); "
                             "atan_3D: oriented azimuth around +x in xy plane (atan2(y,x)) in "
                             "[0,2pi) and polar from +/-z (acos(z)) in [0,pi/2]- NOT BUILT YET.");
  params.addParam<Real>("alpha_tol", 0, "alpha tolerance");
  params.addParam<Real>("hgbalpha_tol", 0, "hgb*alpha tolerance");
  params.addParam<MaterialPropertyName>("hgb", "hgb", "Name of gb switching function.");
  params.addParam<bool>("limit_umag",
                        true,
                        "Limit the $u=\nabla\eta_i - \nabla\eta_j$ based on i and a tolerances, "
                        "else skip the qp/element.");
  return params;
}

GBInclinationBase::GBInclinationBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_name(_op_num),
    // Angular distance to the x axis
    _theta_ij(declareProperty<std::vector<Real>>("theta_ij")),
    _dtheta_dgradeta(declareProperty<std::vector<RealGradient>>("dtheta_dgradeta")),
    _d2theta_dgradeta2(declareProperty<std::vector<RealTensorValue>>("d2theta_dgradeta2")),
    // Polar angle (from z axis)- 90 degrees in 2D with global reference
    _polar_ij(declareProperty<std::vector<Real>>("polar_ij")),
    _dpolar_dgradeta(declareProperty<std::vector<RealGradient>>("dpolar_dgradeta")),
    _d2polar_dgradeta2(declareProperty<std::vector<RealTensorValue>>("d2polar_dgradeta2")),
    _ij_i(declareProperty<std::vector<unsigned int>>("ij_i")),
    _ij_j(declareProperty<std::vector<unsigned int>>("ij_j")),
    _op_is_present(declareProperty<std::vector<unsigned char>>("op_is_present")),
    _op_has_active_pair(declareProperty<std::vector<unsigned char>>("op_has_active_pair")),
    _ug_i(declareProperty<std::vector<unsigned int>>("ug_i")),
    _ug_j(declareProperty<std::vector<unsigned int>>("ug_j")),
    // Grain Tracker/FFC for GB identification
    _gb_case(getParam<MooseEnum>("gb_id_method")),
    _grain_tracker(isParamValid("grain_tracker") ? &getUserObject<GrainTracker>("grain_tracker")
                                                 : nullptr),
    _ffc_tracker(isParamValid("ffc") ? &getUserObject<FeatureFloodCount>("ffc") : nullptr),
    // Thresholds and settings
    _angular_func(getParam<MooseEnum>("angular_func")),
    _alphatol(getParam<Real>("alpha_tol")),
    _hgbatol(getParam<Real>("hgbalpha_tol")),
    _hgb(getMaterialProperty<Real>(getParam<MaterialPropertyName>("hgb"))),
    _no_ij_pairs(declareProperty<bool>("no_ij_pairs")),
    _gtnum(declareProperty<Real>("gt_num")),
    _limit_umag(getParam<bool>("limit_umag"))

{
  if (_op_num < 2)
    mooseError("Inclination properties requires op_num >= 2");

  _vals.resize(_op_num);
  _grad_vals.resize(_op_num);

  for (unsigned int i = 0; i < _op_num; ++i)
  {
    _vals[i] = &coupledValue("v", i);
    _vals_name[i] = coupledName("v", i);
    _grad_vals[i] = &coupledGradient("v", i);
  }
}

void
GBInclinationBase::computeQpProperties()
{
  auto & theta = _theta_ij[_qp];
  auto & dtheta = _dtheta_dgradeta[_qp];
  auto & d2theta = _d2theta_dgradeta2[_qp];

  auto & polar = _polar_ij[_qp];
  auto & dpolar = _dpolar_dgradeta[_qp];
  auto & d2polar = _d2polar_dgradeta2[_qp];

  auto & ij_i = _ij_i[_qp];
  auto & ij_j = _ij_j[_qp];
  auto & ug_i = _ug_i[_qp];
  auto & ug_j = _ug_j[_qp];

  auto & op_is_present = _op_is_present[_qp];
  auto & op_has_active_pair = _op_has_active_pair[_qp];

  theta.clear();
  dtheta.clear();
  d2theta.clear();
  polar.clear();
  dpolar.clear();
  d2polar.clear();
  ij_i.clear();
  ij_j.clear();
  ug_i.clear();
  ug_j.clear();

  op_is_present.assign(_op_num, 0);
  op_has_active_pair.assign(_op_num, 0);

  _gb_pairs.clear();

  _no_ij_pairs[_qp] = true;

  // Find out the number of boundary unique_id and save them

  switch (_gb_case)
  {
    case 0:
    {
      // Grain Tracker Case
      if (!_grain_tracker)
        mooseError(name(), ": gb_id_method=graintracker requires 'grain_tracker' user object.");
      const auto & op_to_grains = (*_grain_tracker).getVarToFeatureVector(_current_elem->id());
      for (auto i : index_range(op_to_grains))
      {
        if (op_to_grains[i] == FeatureFloodCount::invalid_id)
          continue;

        // _gb_ij_list.push_back(i);
        _gb_pairs.push_back({i, op_to_grains[i]});
        op_is_present[i] = 1;
      }
      break;
    }

    case 1:
    {
      // FFC Case
      if (!_ffc_tracker)
        mooseError(name(), ": gb_id_method=ffc requires 'ffc' user object.");
      const auto & op_to_grains = (*_ffc_tracker).getVarToFeatureVector(_current_elem->id());
      for (auto i : index_range(op_to_grains))
      {
        if (op_to_grains[i] == FeatureFloodCount::invalid_id)
          continue;

        // _gb_ij_list.push_back(i);
        _gb_pairs.push_back({i, op_to_grains[i]});
        op_is_present[i] = 1;
      }
      break;
    }

    default:
      mooseError("Unknown gb_id_method = ", _gb_case);
      break;
  }

  // Sort by first value:
  std::sort(_gb_pairs.begin(), _gb_pairs.end()); // sorts by .first by default
  const auto n_active_ops = _gb_pairs.size();
  const auto n_active_pairs = n_active_ops * (n_active_ops - 1) / 2;

  theta.reserve(n_active_pairs);
  dtheta.reserve(n_active_pairs);
  d2theta.reserve(n_active_pairs);
  polar.reserve(n_active_pairs);
  dpolar.reserve(n_active_pairs);
  d2polar.reserve(n_active_pairs);
  ij_i.reserve(n_active_pairs);
  ij_j.reserve(n_active_pairs);
  ug_i.reserve(n_active_pairs);
  ug_j.reserve(n_active_pairs);

  // Grain Tracker check for number of grains
  _gtnum[_qp] = _gb_pairs.size();

  switch (_gb_pairs.size())
  {
    case 0:
      // Zero OP - Skip
      break;
    case 1:
      // One OP - Skip
      break;

    default:
      // do all ij pairs (i<j) if more than 2 vars/features
      for (std::size_t idx1 = 0; idx1 < _gb_pairs.size(); ++idx1)
        for (std::size_t idx2 = idx1 + 1; idx2 < _gb_pairs.size(); ++idx2)
        {
          unsigned int i = _gb_pairs[idx1].first;
          unsigned int j = _gb_pairs[idx2].first;
          unsigned int ugi = _gb_pairs[idx1].second;
          unsigned int ugj = _gb_pairs[idx2].second;
          // const unsigned k = GBPairPacking::pack_upper(i, j);
          // if j-i points outward from lower grain number
          // if i-j like in paper points inward on lower grain number
          // Can also invert it by doing ngb = -ngb
          RealGradient ngb = (*_grad_vals[i])[_qp] - (*_grad_vals[j])[_qp];
          Real a_dist = 0.0;
          RealGradient dtheta_dgetai(0.0, 0.0, 0.0);
          RealTensorValue d2theta_dgetai2(0.0);

          if (ngb.norm() > 1.0e-10)
          {
            RealGradient uxyz(0.0, 0.0, 0.0);
            Real alpha = 0.0;
            uxyz = ngb;
            // uxyz = ngb / ngb.norm();
            ngb /= ngb.norm(); // normalize ngb now that we have uxyz for the un-normalized
            // alpha check
            const bool hov_enabled = (_hgbatol != 0.0);
            const bool inv_enabled = (_alphatol != 0.0);
            bool hov_trigger = false;
            bool inv_trigger = false;
            if (hov_enabled)
            {
              const Real hovtol = _hgb[_qp] / _hgbatol; // safe: _hgbatol != 0 here
              hov_trigger = (uxyz.norm() < hovtol);
              if (_limit_umag && hov_trigger)
              {
                uxyz *= hovtol / uxyz.norm();
                hov_trigger = false;
              }
            }
            if (inv_enabled)
            {
              const Real invtol = 1.0 / _alphatol; // safe: _alphatol != 0 here
              inv_trigger = (uxyz.norm() < invtol);
              if (_limit_umag && inv_trigger)
              {
                uxyz *= invtol / uxyz.norm();
                inv_trigger = false;
              }
            }
            // Combine logic: if both enabled, trip if EITHER trips
            bool alpha_skip = (hov_trigger || inv_trigger);

            // Now calculate the inclination using arc-trig functions
            if (!alpha_skip)
            {
              // moved here for checking if we arent skipping alpha
              alpha = 1 / uxyz.norm();
              // "atan_2D=0 atan_3D=1"
              switch (_angular_func)
              {
                case 0: // atan_2D
                {
                  // ATAN returning full 360 degrees instead of just 0-180 for +/-180
                  // xy plane componenet only
                  const Real x = ngb(0);
                  const Real y = ngb(1);

                  // Value ([-pi, pi] -> [0, 2pi))
                  Real angle = std::atan2(y, x);
                  if (angle < 0.0)
                    angle += 2.0 * libMesh::pi;
                  a_dist = angle;

                  // POLAR angle- in 2D hardcoded to 90 degrees from z axis
                  Real pang = libMesh::pi / 2;

                  // --- First derivatives wrt phi (= components of normalized ngb) ---
                  // angle = atan2(y,x)
                  // dθ/dx = -y/(x^2+y^2), dθ/dy = x/(x^2+y^2), dθ/dz = 0
                  const Real r2 = x * x + y * y;
                  const Real r2_eps = 1e-20; // robust near axis
                  const Real rdenom = std::max(r2, r2_eps);

                  RealGradient dthetadphi(0.0, 0.0, 0.0);
                  dthetadphi(0) = -y / rdenom; // dθ/dx
                  dthetadphi(1) = x / rdenom;  // dθ/dy
                  dthetadphi(2) = 0.0;         // dθ/dz

                  // --- Second derivatives (Hessian) wrt phi ---
                  // ∂²θ/∂x² =  2xy / r2^2
                  // ∂²θ/∂y² = -2xy / r2^2
                  // ∂²θ/∂x∂y = ∂²θ/∂y∂x = (y^2 - x^2) / r2^2
                  const Real inv_r2 = 1.0 / rdenom;
                  const Real inv_r2_2 = inv_r2 * inv_r2;

                  RealTensorValue d2thetadphi2(0.0);
                  d2thetadphi2(0, 0) = 2.0 * x * y * inv_r2_2;
                  d2thetadphi2(1, 1) = -2.0 * x * y * inv_r2_2;
                  const Real dxy = (y * y - x * x) * inv_r2_2;
                  d2thetadphi2(0, 1) = dxy;
                  d2thetadphi2(1, 0) = dxy;
                  // z-rows/cols remain zero

                  // --- Chain rule to grad(eta_i) and its second derivatives ---
                  RealTensorValue dij(1, 0, 0, 0, 1, 0, 0, 0, 1);
                  RankTwoTensor dphi_dgradetai;
                  RankThreeTensor d2phi_dgradetai2;
                  for (unsigned int p = 0; p < 3; ++p)
                  {
                    for (unsigned int q = 0; q < 3; ++q)
                    {
                      dphi_dgradetai(p, q) =
                          alpha * dij(p, q) - alpha * alpha * alpha * uxyz(p) * uxyz(q);
                      for (unsigned int r = 0; r < 3; ++r)
                      {
                        d2phi_dgradetai2(p, q, r) =
                            alpha * alpha * alpha *
                            (3 * alpha * alpha * uxyz(p) * uxyz(q) * uxyz(r) - uxyz(p) * dij(q, r) -
                             uxyz(q) * dij(p, r) - uxyz(r) * dij(p, q));
                      }
                    }
                    for (unsigned int q = 0; q < 3; ++q)
                    {
                      dtheta_dgetai(p) += dthetadphi(q) * dphi_dgradetai(q, p);
                      for (unsigned int r = 0; r < 3; ++r)
                        for (unsigned int s = 0; s < 3; ++s)
                        {
                          d2theta_dgetai2(q, r) +=
                              d2thetadphi2(p, s) * dphi_dgradetai(p, r) * dphi_dgradetai(s, q) +
                              dthetadphi(p) * d2phi_dgradetai2(p, q, r);
                        }
                    }
                  }
                  // // Finally, after all the checks and calcs we save this pair to ij -> k
                  ij_i.push_back(i);
                  ij_j.push_back(j);
                  ug_i.push_back(ugi);
                  ug_j.push_back(ugj);
                  // Kernel reduction precompute
                  op_has_active_pair[i] = 1;
                  op_has_active_pair[j] = 1;
                  // In plane angle
                  theta.push_back(a_dist);
                  dtheta.push_back(dtheta_dgetai);    // this is d/d \nabla \eta_i = -d/dj
                  d2theta.push_back(d2theta_dgetai2); // d2/d \nabla \eta_i^2 = d2/djj = -d2/dij
                  // Polar angle- hardcoded 90 degrees and no derivatives
                  polar.push_back(pang);
                  dpolar.push_back(RealGradient(0.0));
                  d2polar.push_back(RealTensorValue(0.0));
                  break;
                }

                default:
                  mooseError("Unknown angular_func = ", _angular_func);
                  break;
              }
            }
          }
        }
  }
  // Save a bool for if we need to treat this qp as isotropic in children materials
  // true means skip all the vectored calcs at this qp
  _no_ij_pairs[_qp] = _theta_ij[_qp].empty();
}
