//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACInterfaceAnisoGamma.h"

registerMooseObject("PhaseFieldApp", ACInterfaceAnisoGamma);

InputParameters
ACInterfaceAnisoGamma::validParams()
{
  InputParameters params = JvarMapKernelInterface<Kernel>::validParams();
  params.addClassDescription("Gradient energy Allen-Cahn Kernel");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  params.addParam<MaterialPropertyName>(
      "gamma_name", "gamma_aniso", "The gamma used with the kernel");
  params.addParam<MaterialPropertyName>(
      "dgamma_dgradop_name",
      "The derivative of gamma with respect to the gradient of the variable being operated on");
  params.addParam<MaterialPropertyName>(
      "d2gamma_dgradop2_name",
      "The 2nd derivative of gamma with respect to the gradient of the variable being operated on");
  params.addParam<MaterialPropertyName>(
      "dL_dgradop_name",
      "The derivative of L with respect to the gradient of the variable being operated on");
  params.addParam<MaterialPropertyName>(
      "d2L_dgradop2_name",
      "The 2nd derivative of L with respect to the gradient of the variable being operated on");
  // params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "The kappa used with the
  // kernel");
  params.addParam<bool>("variable_L",
                        true,
                        "The mobility is a function of any MOOSE variable (if "
                        "this is set to false L must be constant over the "
                        "entire domain!)");
  params.addParam<bool>("skip_off", false, "Skip the off-diagonal part, for testing/debugging.");
  params.addParam<MaterialPropertyName>("mask_name",
                                        "Name of a MaterialProperty to use as a mask.  "
                                        "If empty, mask = 1.0 everywhere.");
  // params.addRequiredCoupledVar("v",
  //                              "Array of coupled order parameter names for OTHER order
  //                              parameters");
  // params.addRequiredCoupledVar("var_name_base", "Array of coupled variable names");
  params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");
  params.addRequiredCoupledVar(
      "coupled_variables", "Other grain order parameters whose values/gradients enter this kernel");
  return params;
}

ACInterfaceAnisoGamma::ACInterfaceAnisoGamma(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    // Second order values needed for $L=f(\nabla u)$
    _second_u(second()),
    _second_test(secondTest()),
    _second_phi(secondPhi()),
    _L(getMaterialProperty<Real>("mob_name")),
    _gamma(getMaterialProperty<Real>(getParam<MaterialPropertyName>("gamma_name"))),
    _dgammadgrad_op(getMaterialProperty<RealGradient>("dgamma_dgradop_name")),
    _d2gammadgrad_op2(getMaterialProperty<RealTensorValue>("d2gamma_dgradop2_name")),
    _skip_off(getParam<bool>("skip_off")),
    _mask(isParamValid("mask_name") ? &getMaterialProperty<Real>("mask_name") : nullptr),
    _mask_tf(isParamValid("mask_name")),
    // Grain OP list as per ACGrGrPoly
    _grain_ids(),
    // _grain_set(),
    // _op_num(coupledComponents("v")),
    // _vals(coupledValues("v")),
    // _vals_var(coupledIndices("v")),
    _var_name_base(getParam<std::string>("var_name_base")),
    _variable_L(getParam<bool>("variable_L")),
    _dLdop(_variable_L ? &getMaterialPropertyDerivative<Real>("mob_name", _var.name()) : nullptr),
    _d2Ldop2(_variable_L
                 ? &getMaterialPropertyDerivative<Real>("mob_name", _var.name(), _var.name())
                 : nullptr),
    _dLdgrad_op(_variable_L ? &getMaterialProperty<RealGradient>("dL_dgradop_name") : nullptr),
    _d2Ldgrad_op2(_variable_L ? &getMaterialProperty<RealTensorValue>("d2L_dgradop2_name")
                              : nullptr),
    _d2Ldopdgrad_op(
        _variable_L ? &getMaterialPropertyDerivative<RealGradient>("dL_dgradop_name", _var.name())
                    : nullptr),
    _dLdarg(_n_args),
    _d2Ldargdop(_n_args),
    _d2Ldarg2(_n_args),
    _dLdgradarg(_n_args),
    _d2Ldgradargdop(_n_args),
    _d2Ldargdgradop(_n_args),
    _d2Ldgradarg2(_n_args),
    _d2Ldgradargdarg(_n_args),
    // _dkappadarg(_n_args),
    // _arg(_n_args),
    _eta_vals(),
    _gradarg(_n_args),
    _second_arg(_n_args)
{
  // // Pull a mask value if there is one
  // const std::string & mask_name = getParam<std::string>("mask_name");
  // if (!mask_name.empty())
  // {
  //   // mooseWarning("Mask Name is ", mask_name);
  //   _mask = &getMaterialProperty<Real>(mask_name);
  //   mooseWarning("Mask name is ", mask_name);
  //   mooseWarning("Mask is ", _mask);
  // }
  // else
  // {
  //   mooseWarning("Mask name is empty!");
  // }

  // Get mobility and kappa derivatives and coupled variable gradients
  // mooseWarning("NAMEBASE: ", _var_name_base);
  // std::vector<unsigned int> _grain_ids;
  // mooseWarning("Coupled stuff: ", _n_args);
  for (unsigned int i = 0; i < _n_args; ++i)
  {
    // mooseWarning("Into initial");
    MooseVariable * ivar = _coupled_standard_moose_vars[i];
    const VariableName iname = ivar->name();
    // mooseWarning("coupled name: ", iname);
    if (iname == _var.name())
    {
      if (isCoupled("args"))
        paramError("args",
                   "The kernel variable should not be specified in the coupled `args` parameter.");
      else
        paramError("coupled_variables",
                   "The kernel variable should not be specified in the coupled `coupled_variables` "
                   "parameter.");
    }

    if (_variable_L)
    {
      _dLdarg[i] = &getMaterialPropertyDerivative<Real>("mob_name", i);
      _d2Ldargdop[i] = &getMaterialPropertyDerivative<Real>("mob_name", iname, _var.name());
      _d2Ldargdgradop[i] = &getMaterialPropertyDerivative<RealGradient>("dL_dgradop_name", iname);
      _d2Ldarg2[i].resize(_n_args);
      for (unsigned int j = 0; j < _n_args; ++j)
        _d2Ldarg2[i][j] = &getMaterialPropertyDerivative<Real>("mob_name", i, j);

      _gradarg[i] = &(ivar->gradSln());
      _second_arg[i] = &(ivar->secondSln());
    }

    // ID the number at the end of the args to determine which dLdgrad to use
    // Check if coupled var starts with the var_name_base (grain op)
    if (MooseUtils::beginsWith(iname, _var_name_base) && iname != _var.name())
    {
      _grain_ids.push_back(i);
      _eta_vals.push_back(&coupledValue("coupled_variables", i));
      if (_variable_L)
      {
        // Pull the number at the end of the grain op
        const std::size_t pos = iname.find_last_not_of("0123456789");
        if (pos == std::string::npos || pos + 1 == iname.size())
          mooseError("Variable '",
                     iname,
                     "' does not end in a numeric suffix, despite starting with var_name_base.");
        const std::string digits = iname.substr(pos + 1);
        // mooseWarning("digits = ", digits);
        _dLdgradarg[i] = &getMaterialProperty<RealGradient>("dLdgrad_eta_" + digits);
        _d2Ldgradarg2[i] = &getMaterialProperty<RealTensorValue>("d2Ldgrad_eta2_" + digits);
        _d2Ldgradargdarg[i].resize(_n_args);
        // _d2Ldgradargdarg[i] =
        //     &getMaterialPropertyDerivative<RealGradient>("dLdgrad_eta_" + digits, iname);
        _d2Ldgradargdop[i] =
            &getMaterialPropertyDerivative<RealGradient>("dLdgrad_eta_" + digits, _var.name());
        for (unsigned int j = 0; j < _n_args; ++j)
          _d2Ldgradargdarg[i][j] = &getMaterialPropertyDerivative<RealGradient>(
              "dLdgrad_eta_" + digits, _coupled_standard_moose_vars[j]->name());
      }
    }
  }
}

void
ACInterfaceAnisoGamma::initialSetup()
{
  validateCoupling<Real>("mob_name");
}

RealGradient
ACInterfaceAnisoGamma::gradL() // Includes grad op dependence
{
  RealGradient g = _grad_u[_qp] * (*_dLdop)[_qp];
  g += _second_u[_qp] * (*_dLdgrad_op)[_qp];
  for (unsigned int i = 0; i < _n_args; ++i)
  {
    g += (*_gradarg[i])[_qp] * (*_dLdarg[i])[_qp];
    g += (*_second_arg[i])[_qp] * (*_dLdgradarg[i])[_qp];
  }
  return g;
}

RealGradient
ACInterfaceAnisoGamma::nablaLPsi() // RH $L \nabla \psi$
{
  // sum is the product rule gradient \f$ \nabla (L\psi) \f$
  RealGradient sum = _L[_qp] * _grad_test[_i][_qp];

  if (_variable_L)
    sum += gradL() * _test[_i][_qp];

  return sum;
}

Real
ACInterfaceAnisoGamma::sumSqEtaj()
{
  // Sum all other (grain?) order parameters
  Real SumOPj = 0.0;
  // for (unsigned int j : _grain_ids)
  for (unsigned int k = 0; k < _eta_vals.size(); ++k)
  {
    const Real eta_j = (*_eta_vals[k])[_qp];
    SumOPj += eta_j * eta_j;
    // SumOPj += (*_coupled_standard_moose_vars[j]).sln()[_qp] *
    //           (*_coupled_standard_moose_vars[j]).sln()[_qp];
  }
  return SumOPj;
}

Real
ACInterfaceAnisoGamma::computeQpResidual()
{
  if (_mask_tf)
    return (*_mask)[_qp] * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() * nablaLPsi();
  else
    return _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() * nablaLPsi();
}

Real
ACInterfaceAnisoGamma::computeQpJacobian()
{
  // dsum is the derivative of R wrt u
  Real ddir = 2 * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() * _phi[_j][_qp] * nablaLPsi();

  // dind is the derivative of R wrt \nabla u
  Real dind =
      _u[_qp] * _u[_qp] * _d2gammadgrad_op2[_qp] * sumSqEtaj() * _grad_phi[_j][_qp] * nablaLPsi();

  if (_variable_L)
  {
    // Grad L partials
    static const RealTensorValue I(1, 0, 0, 0, 1, 0, 0, 0, 1);
    // The direct u pieces
    RealGradient dgradLdu =
        (*_d2Ldop2)[_qp] * _grad_u[_qp] + (*_d2Ldopdgrad_op)[_qp] * _second_u[_qp];
    RealTensorValue dgradLdgradu = libMesh::outer_product((*_d2Ldopdgrad_op)[_qp], _grad_u[_qp]) +
                                   I * (*_dLdop)[_qp] + (*_d2Ldgrad_op2)[_qp] * _second_u[_qp];
    // The cross terms with eta/gradeta dependence in grad L
    for (unsigned int i = 0; i < _n_args; ++i)
    {
      dgradLdu += (*_d2Ldargdop[i])[_qp] * (*_gradarg[i])[_qp] +
                  (*_d2Ldgradargdop[i])[_qp] * (*_second_arg[i])[_qp];
      dgradLdgradu += libMesh::outer_product((*_d2Ldargdgradop[i])[_qp], (*_gradarg[i])[_qp]) -
                      (*_d2Ldgrad_op2)[_qp] * (*_second_arg[i])[_qp];
    }

    // Direct L dependence
    ddir += (*_dLdop)[_qp] * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
            _phi[_j][_qp] * _grad_test[_i][_qp];
    // Direct grad L dependence
    ddir += dgradLdu * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() * _phi[_j][_qp] *
            _test[_i][_qp];
    // Indirect L dependence (of grad u)
    dind += (*_dLdgrad_op)[_qp] * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
            _grad_phi[_j][_qp] * _grad_test[_i][_qp];
    // Indirect grad L dependence (of grad u)
    dind += dgradLdgradu * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
            _grad_phi[_j][_qp] * _test[_i][_qp];
  }

  if (_mask_tf)
    return (*_mask)[_qp] * (ddir + dind);
  else
    return (ddir + dind);
}

Real
ACInterfaceAnisoGamma::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_skip_off)
  {
    return 0.0;
  }
  else
  {
    // get the coupled variable jvar is referring to
    const unsigned int cvar = mapJvarToCvar(jvar);

    // Check if the offdiag wrt variable is one of the grain ops
    bool is_grain = std::find(_grain_ids.begin(), _grain_ids.end(), cvar) != _grain_ids.end();
    if (is_grain)
    {
      // Direct dependence on arg
      Real ddir = 2 * (*_coupled_standard_moose_vars[cvar]).sln()[_qp] * _u[_qp] * _u[_qp] *
                  _dgammadgrad_op[_qp] * _phi[_j][_qp] * nablaLPsi();
      // Indirect dependence on grad_arg
      Real dind = _u[_qp] * _u[_qp] * (-_d2gammadgrad_op2[_qp]) * sumSqEtaj() * _grad_phi[_j][_qp] *
                  nablaLPsi();

      if (_variable_L)
      {
        // Grad L partials
        static const RealTensorValue I(1, 0, 0, 0, 1, 0, 0, 0, 1);
        // The direct u pieces
        RealGradient dgradLdarg = (*_d2Ldargdop[cvar])[_qp] * _grad_u[_qp] +
                                  (*_d2Ldargdgradop[cvar])[_qp] * _second_u[_qp];
        RealTensorValue dgradLdgradarg =
            libMesh::outer_product((*_d2Ldgradargdop[cvar])[_qp], _grad_u[_qp]) -
            (*_d2Ldgrad_op2)[_qp] * _second_u[_qp] + I * (*_dLdarg[cvar])[_qp] -
            (*_d2Ldgradarg2[cvar])[_qp] * (*_second_arg[cvar])[_qp];
        // Could have more of the dgradLdgradarg into cross terms if those were defined that way
        // The cross terms with eta/gradeta dependence in grad L
        for (unsigned int i = 0; i < _n_args; ++i)
        {
          dgradLdarg += (*_d2Ldarg2[cvar][i])[_qp] * (*_gradarg[i])[_qp] +
                        (*_d2Ldgradargdarg[i][cvar])[_qp] * (*_second_arg[i])[_qp];
          dgradLdgradarg +=
              libMesh::outer_product((*_d2Ldgradargdarg[cvar][i])[_qp], (*_gradarg[i])[_qp]);
        }
        // Now combine the offdiag grad L pieces
        // Direct L dependence
        ddir += (*_dLdarg[cvar])[_qp] * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
                _phi[_j][_qp] * _grad_test[_i][_qp];
        // Direct grad L dependence
        ddir += dgradLdarg * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
                _phi[_j][_qp] * _test[_i][_qp];
        // Indirect L dependence (of grad arg)
        dind += (*_dLdgradarg[cvar])[_qp] * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
                _grad_phi[_j][_qp] * _grad_test[_i][_qp];
        // Indirect grad L dependence (of grad arg)
        dind += dgradLdgradarg * _u[_qp] * _u[_qp] * _dgammadgrad_op[_qp] * sumSqEtaj() *
                _grad_phi[_j][_qp] * _test[_i][_qp];
      }

      // Output the grain_op based offdiagonal
      if (_mask_tf)
        return (*_mask)[_qp] * (ddir + dind);
      else
        return (ddir + dind);
    }

    // Non-grain OP offdiag
    return 0.0;
  }
}
