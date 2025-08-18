#include "GGInclinationVector.h"

registerMooseObject("PhaseFieldApp", GGInclinationVector);

InputParameters
GGInclinationVector::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Measures the inclination as a vector along the grain boundaries.");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  // params.addParam<std::string>(
  //     "grain_tracker", "", "The GrainTracker UserObject to get values from.");
  params.addParam<UserObjectName>("grain_tracker",
                                  "The GrainTracker UserObject to get values from.");
  params.addParam<UserObjectName>("ffc", "The FFC UserObject to get values from.");
  params.addParam<MaterialPropertyName>("hgb", "Name of GB switching function.");
  // params.addParam<std::string>(
  //     "hgb", "", "The GB switching function if manually defining it (GBID=SWITCH).");
  MooseEnum gbident("GRAINTRACKER=0 HGB=1 SWITCH=2 FFC=3", "GRAINTRACKER");
  params.addParam<MooseEnum>("gb_id_method", gbident, "Which approach to use to select the GB.");
  params.addParam<Real>(
      "hgb_threshold", 0.6, "Cutoff for GB switching function to select GB region.");
  return params;
}

GGInclinationVector::GGInclinationVector(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    // : Material(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _vals_name(_op_num),
    // Inclination vector for polar plots
    _inclination_vector(declareProperty<RealGradient>("inclination_vector")),
    _ang_dist(declareProperty<Real>("ang_dist")),
    // Grain Tracker/EBSD for GB identification
    _gb_id_type(getParam<MooseEnum>("gb_id_method")),
    _grain_tracker(isParamValid("grain_tracker") ? &getUserObject<GrainTracker>("grain_tracker")
                                                 : nullptr),
    _ffc_tracker(isParamValid("ffc") ? &getUserObject<FeatureFloodCount>("ffc") : nullptr),
    _hgb_external(isParamValid("hgb") ? &getMaterialProperty<Real>("hgb") : nullptr),
    _hgb_threshold(getParam<Real>("hgb_threshold"))
{
  if (_op_num == 0)
    mooseError("Model requires op_num > 0");

  // Check for GT user object if using that case
  // std::string gt_name = getParam<std::string>("grain_tracker");
  // mooseWarning("GTName = ", gt_name);
  // if ((_gb_id_type == 0) && (!gt_name.empty()))
  //   _grain_tracker = &getUserObject<GrainTracker>(gt_name);
  // else if ((_gb_id_type == 3) && (!gt_name.empty()))
  //   _ffc_tracker = &getUserObject<FeatureFloodCount>(gt_name);
  // else if ((_gb_id_type == 0) && (gt_name.empty()))
  //   mooseError(
  //       "GGInclinationVector requires a grain_tracker specified if gb_id_method =
  //       GRAINTRACKER.");
  // else if ((_gb_id_type == 3) && (gt_name.empty()))
  //   mooseError("GGInclinationVector requires a grain_tracker specified if gb_id_method = FFC.");

  // // Check for hgb input if using case 2 (SWITCH)
  // std::string hgb_name = getParam<std::string>("hgb");
  // if ((_gb_id_type == 2) && (!hgb_name.empty()))
  //   _hgb_external = &getMaterialProperty<Real>(hgb_name);
  // else if ((_gb_id_type == 2) && (hgb_name.empty()))
  //   mooseError("GGInclinationVector requires a hgb specified if gb_id_method = SWITCH.");

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
GGInclinationVector::computeQpProperties()
{
  _gb_ij_pairs.clear();
  _hgb_pairs.clear();
  _inc_vec_pairs.clear();

  switch (_gb_id_type)
  {
    case 0:
    {
      // Grain Tracker method
      const auto & op_to_grains = (*_grain_tracker).getVarToFeatureVector(_current_elem->id());
      for (auto i : index_range(op_to_grains))
      {
        if (op_to_grains[i] == FeatureFloodCount::invalid_id)
          continue;

        _gb_ij_pairs.push_back(i);
      }
      break;
    }
    case 1:
    {
      // HGB method
      // Check hgb threshold
      Real hbulk = 0.0;
      Real hgb_thresh = 0.0;
      for (unsigned int i = 0; i < _op_num; ++i)
        hbulk += (*_vals[i])[_qp] * (*_vals[i])[_qp];
      hgb_thresh = 4 * (1 - hbulk) * (1 - hbulk);
      if (hgb_thresh > _hgb_threshold)
      {
        _gb_ij_pairs.resize(_op_num);
        std::iota(_gb_ij_pairs.begin(), _gb_ij_pairs.end(), 0);
      }
      else
      {
        _gb_ij_pairs.clear();
      }
      break;
    }
    case 2:
    {
      // Input HGB method
      // Check hgb threshold
      if ((*_hgb_external)[_qp] > _hgb_threshold)
      {
        _gb_ij_pairs.resize(_op_num);
        std::iota(_gb_ij_pairs.begin(), _gb_ij_pairs.end(), 0);
      }
      else
      {
        _gb_ij_pairs.clear();
      }
      break;
    }
    case 3:
    {
      // FFC method
      const auto & op_to_grains = (*_ffc_tracker).getVarToFeatureVector(_current_elem->id());
      for (auto i : index_range(op_to_grains))
      {
        if (op_to_grains[i] == FeatureFloodCount::invalid_id)
          continue;

        _gb_ij_pairs.push_back(i);
      }
      break;
    }
    default:
      mooseError("Unknown gb_id_method = ", _gb_id_type);
  }

  // Make a copy and sort
  _gb_ij_sorted = _gb_ij_pairs;
  std::sort(_gb_ij_sorted.begin(), _gb_ij_sorted.end());

  switch (_gb_ij_pairs.size())
  {
    case 0:
      _inclination_vector[_qp] = RealGradient(0.0, 0.0, 0.0);
      _ang_dist[_qp] = -1;
      break;
    case 1:
      _inclination_vector[_qp] = RealGradient(0.0, 0.0, 0.0);
      _ang_dist[_qp] = -1;
      break;
    default:
      // do all ij pairs if more than 2 vars/features
      for (std::size_t idx1 = 0; idx1 < _gb_ij_sorted.size(); ++idx1)
        for (std::size_t idx2 = idx1 + 1; idx2 < _gb_ij_sorted.size(); ++idx2)
        {
          unsigned int i = _gb_ij_sorted[idx1];
          unsigned int j = _gb_ij_sorted[idx2];
          RealGradient ngb = (*_grad_vals[i])[_qp] - (*_grad_vals[j])[_qp];
          // _hgb_pairs.push_back((*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] *
          //                      (*_vals[j])[_qp]);
          if (ngb.norm() > 1.0e-10)
          {
            ngb /= ngb.norm(); // Really dont think this should be here?
            _inc_vec_pairs.push_back(ngb);
            _hgb_pairs.push_back((*_vals[i])[_qp] * (*_vals[i])[_qp] * (*_vals[j])[_qp] *
                                 (*_vals[j])[_qp]);
          }
          else
          {
            ngb = 0.0;
          }
        }

      // Now take the assumed multiple and convert
      if (_hgb_pairs.size() > 0)
      {
        RealGradient numer(0.0, 0.0, 0.0);
        Real denom = 0.0;
        Real ang_numer = 0.0;
        for (std::size_t n = 0; n < _hgb_pairs.size(); ++n)
        {
          RealGradient ngb_n = _inc_vec_pairs[n];
          numer += ngb_n * _hgb_pairs[n];
          denom += _hgb_pairs[n];
          Real R = std::sqrt((ngb_n(1) * ngb_n(1)) + (ngb_n(2) * ngb_n(2)));
          Real a_dist = std::atan2(R, ngb_n(0));
          ang_numer += a_dist * _hgb_pairs[n];
        }
        if (denom > 1.0e-10)
        {
          _inclination_vector[_qp] = numer / denom;
          _ang_dist[_qp] = ang_numer / denom;
        }
        else
        {
          _inclination_vector[_qp] = RealGradient(0.0, 0.0, 0.0);
          _ang_dist[_qp] = -1;
        }
      }
      else
      {
        _inclination_vector[_qp] = RealGradient(0.0, 0.0, 0.0);
        _ang_dist[_qp] = -1;
      }
  }
}
