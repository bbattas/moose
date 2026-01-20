//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElementalGammaMaterial.h"
#include "libmesh/quadrature.h"

registerMooseObject("PhaseFieldApp", ElementalGammaMaterial);

InputParameters
ElementalGammaMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Generates a diffusion tensor to distinguish between the bulk, grain "
                             "boundaries, and surfaces");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<MaterialPropertyName>(
      "gb_energy_iso_name", "sigma_iso", "Isotropic GB energy before inclination dependence.");
  params.addParam<MaterialPropertyName>(
      "kappa", "kappa", "Gradient energy constant kappa material name.");
  params.addParam<Real>(
      "free_energy_m", 1, "Free energy function constant m (or mu in PF kernels).");
  params.addParam<MaterialPropertyName>("L0", "L0", "Name of L prefactor/iso value.");
  params.addParam<MaterialPropertyName>(
      "L_qp", "L_qp", "Name of L with some qps dropped to be averaged.");
  params.addParam<bool>("aniso_L", false, "Is AC mobility L an inclination dependent variable.");
  params.addParam<bool>("well", false, "Using well function for inclination?");
  params.addParam<bool>(
      "skip", false, "Skip all the averaging calculations (for testing purposes).");
  return params;
}

ElementalGammaMaterial::ElementalGammaMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _gamma_in(getMaterialProperty<Real>("gamma_qp")),
    _gamma_out(declareProperty<Real>("gamma_asymm")),
    _no_ij_pairs(getMaterialProperty<bool>("no_ij_pairs")),
    _elem_no_ij(declareProperty<bool>("elem_no_ij")),
    _int_noij(declareProperty<Real>("int_noij")),
    _gbe_iso(getMaterialProperty<Real>("gb_energy_iso_name")),
    _kappa(getMaterialProperty<Real>(getParam<MaterialPropertyName>("kappa"))),
    _const_m(getParam<Real>("free_energy_m")),
    _aniso_L(getParam<bool>("aniso_L")),
    _L0(getMaterialProperty<Real>(getParam<MaterialPropertyName>("L0"))),
    _L_qp(_aniso_L ? &getMaterialProperty<Real>(getParam<MaterialPropertyName>("L_qp")) : nullptr),
    _L_elem(_aniso_L ? &declareProperty<Real>("L") : nullptr),
    _well(getParam<bool>("well")),
    _int_width_in(getMaterialProperty<Real>("int_width_qp")),
    _int_width_out(declareProperty<Real>("int_width")),
    _skip(getParam<bool>("skip"))
{
  // if (_op_num == 0)
  //   mooseError("Model requires op_num > 0");

  // _vals.resize(_op_num);
  // _grad_vals.resize(_op_num);
  // for (unsigned int i = 0; i < _op_num; ++i)
  // {
  //   _vals[i] = &coupledValue("v", i);
  //   _grad_vals[i] = &coupledGradient("v", i);
  //   _vals_name[i] = coupledName("v", i);
  //   if (!isCoupledConstant(_vals_name[i]))
  //     _dDdeta[i] = &declarePropertyDerivative<RealTensorValue>(_diffusivity_name, _vals_name[i]);
  // }
}

void
ElementalGammaMaterial::computeProperties()
{
  // If skipping, just copy input gamma/int_width to outputs on all qps
  if (_skip)
  {
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    {
      _gamma_out[_qp] = _gamma_in[_qp];
      _int_width_out[_qp] = _int_width_in[_qp];

      // Set flags to some value
      _elem_no_ij[_qp] = false;
      _int_noij[_qp] = 0;

      if (_aniso_L)
        (*_L_elem)[_qp] = (*_L_qp)[_qp];
    }

    return; // Skip all remaining computations
  }

  // Hard-coded coefficients (for poly gamma)
  constexpr Real a1 = -3.0944; // coefficient for g2^4
  constexpr Real a2 = -1.8169; // coefficient for g2^3
  constexpr Real a3 = 10.323;  // coefficient for g2^2
  constexpr Real a4 = -8.1819; // coefficient for g2
  constexpr Real a5 = 2.0033;  // constant term

  bool elem_skip = false;
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    const bool qp_noij = _no_ij_pairs[_qp];
    if (qp_noij)
    {
      elem_skip = true;
      break;
    }
  }

  bool gamma_change = false;
  Real min_gamma = 0.0;
  Real max_gamma = 0.0;
  Real min_iw = 0.0;
  Real max_iw = 0.0;
  if (_well)
  {
    // initialize min/max from first qp
    min_gamma = _gamma_in[0];
    max_gamma = _gamma_in[0];

    min_iw = _int_width_in[0];
    max_iw = _int_width_in[0];

    for (_qp = 1; _qp < _qrule->n_points(); ++_qp)
    {
      // check if any qp is different from the first one
      if (_gamma_in[_qp] != _gamma_in[0])
        gamma_change = true;

      // update min / max
      if (_gamma_in[_qp] < min_gamma)
      {
        min_gamma = _gamma_in[_qp];
        max_iw = _int_width_in[_qp];
      }
      if (_gamma_in[_qp] > max_gamma)
      {
        max_gamma = _gamma_in[_qp];
        min_iw = _int_width_in[_qp];
      }
    }
  }

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    _elem_no_ij[_qp] = elem_skip;
    if (elem_skip)
    {
      _int_noij[_qp] = 1;
      // Gamma
      Real g = _gbe_iso[_qp] / (std::sqrt(_kappa[_qp] * _const_m));
      Real g2 = g * g;
      Real pg = (((a1 * g2 + a2) * g2 + a3) * g2 + a4) * g2 + a5;
      _gamma_out[_qp] = 1 / pg;
      if (_aniso_L)
        (*_L_elem)[_qp] = _L0[_qp];
      // IW
      Real f0_int =
          (((((0.0788 * pg - 0.4955) * pg + 1.2244) * pg - 1.5281) * pg + 1.0686) * pg - 0.5563) *
              pg +
          0.2907;
      _int_width_out[_qp] = (std::sqrt(_kappa[_qp] / _const_m)) * (std::sqrt(1 / f0_int));
    }
    else if (gamma_change)
    {
      _int_noij[_qp] = 2;
      _gamma_out[_qp] = max_gamma;
      if (_aniso_L)
        (*_L_elem)[_qp] = (*_L_qp)[_qp];
      _int_width_out[_qp] = min_iw;
    }
    else
    {
      _int_noij[_qp] = 0;
      _gamma_out[_qp] = _gamma_in[_qp];
      if (_aniso_L)
        (*_L_elem)[_qp] = (*_L_qp)[_qp];
      _int_width_out[_qp] = _int_width_in[_qp];
    }
  }
}
