#include "GBInclinationBase.h"

registerMooseObject("PhaseFieldApp", GBInclinationBase);

InputParameters
GBInclinationBase::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Base material to determine inclination dependent properties for AGG. "
                             "Builds up to the inclination and its derivatives wrt $\nabla\eta$.");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  MooseEnum gb_id_method("graintracker=0 ffc=1 all=2", "graintracker");
  params.addParam<MooseEnum>(
      "gb_id_method", gb_id_method, "Which GB OP identification/selection option.");
  params.addParam<UserObjectName>("grain_tracker",
                                  "The GrainTracker UserObject to get values from.");
  params.addParam<UserObjectName>("ffc", "The FFC UserObject to get values from.");
  MooseEnum angular_func("atan_2D=0 atan_3D=1 acos=2 atan_half=3", "atan_2D");
  params.addParam<MooseEnum>(
      "angular_func",
      angular_func,
      "Which angular distance function to use. "
      "atan_2D: oriented angle atan2(y,x) in [0,2pi); "
      "atan_3D: oriented azimuth around +x (atan2(z,y)) in [0,2pi) *UNTESTED*;"
      "acos: acos(x) in [0,pi]; "
      "atan_half: atan2(sqrt(y^2+z^2), x) in [0,pi].");
  params.addParam<Real>("intol", 100, "hgbalpha tolerance");
  params.addParam<Real>("altol", 100, "alpha tolerance");
  params.addParam<Real>("gt_tol", 0.001, "alpha tolerance");
  params.addParam<MaterialPropertyName>("hgb", "hgb", "Name of gb switching function.");
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
    _ij_i(declareProperty<std::vector<unsigned int>>("ij_i")),
    _ij_j(declareProperty<std::vector<unsigned int>>("ij_j")),
    // Grain Tracker/FFC for GB identification
    _gb_case(getParam<MooseEnum>("gb_id_method")),
    _grain_tracker(isParamValid("grain_tracker") ? &getUserObject<GrainTracker>("grain_tracker")
                                                 : nullptr),
    _ffc_tracker(isParamValid("ffc") ? &getUserObject<FeatureFloodCount>("ffc") : nullptr),
    _gt_tol(getParam<Real>("gt_tol")),
    _gtnum(declareProperty<Real>("gtnum")),
    _angular_func(getParam<MooseEnum>("angular_func")),
    _intol(getParam<Real>("intol")),
    _altol(getParam<Real>("altol")),
    _hgb(getMaterialProperty<Real>(getParam<MaterialPropertyName>("hgb"))),
    _no_ij_pairs(declareProperty<bool>("no_ij_pairs")),
    _testout3(declareProperty<Real>("testout3")),
    _aval(declareProperty<RealGradient>("a_val")),
    _ival(declareProperty<RealGradient>("i_val")),
    _acut(declareProperty<RealGradient>("a_cut")),
    _icut(declareProperty<RealGradient>("i_cut"))

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
  // Flatpack vector matrices definitions
  const unsigned num_pairs = GBPairPacking::count_upper(_op_num);
  auto & theta = _theta_ij[_qp];
  auto & dtheta = _dtheta_dgradeta[_qp];
  auto & d2theta = _d2theta_dgradeta2[_qp];
  // Size once per QP (no-op if already correct); then initialize values
  if (theta.size() != num_pairs)
    theta.resize(num_pairs);
  if (dtheta.size() != num_pairs)
    dtheta.resize(num_pairs);
  if (d2theta.size() != num_pairs)
    d2theta.resize(num_pairs);
  std::fill(theta.begin(), theta.end(), -1.0);
  std::fill(dtheta.begin(), dtheta.end(), RealGradient(0.0));
  std::fill(d2theta.begin(), d2theta.end(), RealTensorValue(0.0));

  auto & k2i = _ij_i[_qp];
  auto & k2j = _ij_j[_qp];
  if (k2i.size() != num_pairs)
    k2i.resize(num_pairs);
  if (k2j.size() != num_pairs)
    k2j.resize(num_pairs);
  // optional: initialize to a sentinel (here UINT_MAX)
  std::fill(k2i.begin(), k2i.end(), 4);
  std::fill(k2j.begin(), k2j.end(), 4);

  // Find out the number of boundary unique_id and save them
  _gb_ij_list.clear();
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

        _gb_ij_list.push_back(i);
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

        _gb_ij_list.push_back(i);
      }
      break;
    }

    case 2:
    {
      // Variable threshold cutoffs- BROKEN/FAILS
      for (unsigned int i = 0; i < _op_num; ++i)
      {
        if ((*_vals[i])[_qp] >= _gt_tol)
        {
          _gb_ij_list.push_back(i);
        }
      }
      break;
    }

    default:
      mooseError("Unknown gb_case = ", _gb_case);
      break;
  }

  _testout3[_qp] = 0.0;
  // _aval[_qp] = RealGradient(0.0);
  // _ival[_qp] = RealGradient(0.0);
  RealGradient aval(0.0);
  RealGradient ival(0.0);
  RealGradient acut(0.0);
  RealGradient icut(0.0);

  _gtnum[_qp] = _gb_ij_list.size();
  std::sort(_gb_ij_list.begin(), _gb_ij_list.end());

  switch (_gb_ij_list.size())
  {
    case 0:
      // _inclination_distance[_qp] = -1.0; // angular distance out of 0-2pi range (not counted)
      // _inclination[_qp] = 1.0; // f = 1 + cos() = 1
      break;
    case 1:
      // _inclination_distance[_qp] = -1.0; // angular distance out of 0-2pi range (not counted)
      // _inclination[_qp] = 1.0; // f = 1 + cos() = 1
      break;

    default:
      // do all ij pairs (i>j) if more than 2 vars/features
      for (std::size_t idx1 = 0; idx1 < _gb_ij_list.size(); ++idx1)
        for (std::size_t idx2 = idx1 + 1; idx2 < _gb_ij_list.size(); ++idx2)
        {
          unsigned int i = _gb_ij_list[idx1];
          unsigned int j = _gb_ij_list[idx2];
          const unsigned k = GBPairPacking::pack_upper(i, j);
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
            // DEBUG OUT
            if (k < 3)
            {
              aval(k) = _hgb[_qp] / uxyz.norm();
              ival(k) = 1 / uxyz.norm();
              // values
            }
            // alpha check
            const bool hov_enabled = (_altol != 0.0);
            const bool inv_enabled = (_intol != 0.0);
            bool hov_trigger = false;
            bool inv_trigger = false;
            if (hov_enabled)
            {
              const Real hovtol = _hgb[_qp] / _altol; // safe: _altol != 0 here
              hov_trigger = (uxyz.norm() < hovtol);
            }
            if (inv_enabled)
            {
              const Real invtol = 1.0 / _intol; // safe: _intol != 0 here
              inv_trigger = (uxyz.norm() < invtol);
            }
            // Combine logic: if both enabled, trip if EITHER trips
            bool alpha_skip = (hov_trigger || inv_trigger);
            // Check whats triggering where
            if (idx1 == 0 && idx2 == 1)
            {
              if (hov_trigger && inv_trigger)
                _testout3[_qp] = 3;
              else if (hov_trigger && !inv_trigger)
                _testout3[_qp] = 1;
              else if (!hov_trigger && inv_trigger)
                _testout3[_qp] = 2;
            }
            // DEBUG
            if (k < 3)
            {
              acut(k) = aval(k);
              icut(k) = ival(k);
              if (hov_trigger && inv_trigger)
              {
                acut(k) = -1;
                icut(k) = -1;
              }
              else if (hov_trigger && !inv_trigger)
                acut(k) = -1;
              else if (!hov_trigger && inv_trigger)
                icut(k) = -1;
            }

            // bool alpha_skip = false;
            // if (((uxyz.norm() < hovtol) && (_altol != 0.0)) ||
            //     ((uxyz.norm() < invtol) && (_intol != 0.0))) // hgb (altol)
            // {
            //   alpha = 0.0;
            //   alpha_skip = true;
            //   if (k == 1)
            //     _testout3[_qp] = 1;
            // }
            // else
            // {
            //   alpha = 1 / uxyz.norm();
            // }

            // Now calculate the inclination using arc-trig functions
            if (!alpha_skip)
            {
              // moved here for checking if we arent skipping alpha
              alpha = 1 / uxyz.norm();
              // "atan_2D=0 atan_3D=1 acos=2 atan_half=3"
              switch (_angular_func)
              {
                case 0: // atan_2D
                {
                  // ATAN returning full 360 degrees instead of just 0-180 for +/-180
                  // --- Oriented azimuth around +x, referenced from +y, increasing toward +z ---
                  const Real x = ngb(0);
                  const Real y = ngb(1);

                  // Value ([-pi, pi] -> [0, 2pi))
                  Real angle = std::atan2(y, x);
                  if (angle < 0.0)
                    angle += 2.0 * libMesh::pi;
                  a_dist = angle;

                  // --- First derivatives wrt phi (= components of normalized ngb) ---
                  // angle = atan2(z, y)
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
                  // COPIED FROM PREVIOUS PARTS
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
                  // Finally, after all the checks and calcs we save this pair to ij -> k
                  theta[k] = a_dist;
                  dtheta[k] = dtheta_dgetai;    // this is d/d \nabla \eta_i = -d/dj
                  d2theta[k] = d2theta_dgetai2; // this is d2/d \nabla \eta_i^2 = d2/djj = -d2/dij
                  k2i[k] = i;
                  k2j[k] = j;

                  // _inclination[_qp] = 1.0;
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
  _no_ij_pairs[_qp] =
      std::all_of(theta.begin(), theta.end(), [](const Real v) { return v == -1.0; });
  // DEBUG
  _aval[_qp] = aval;
  _ival[_qp] = ival;
  _acut[_qp] = acut;
  _icut[_qp] = icut;
}
