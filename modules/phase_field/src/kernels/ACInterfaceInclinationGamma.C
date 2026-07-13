#include "ACInterfaceInclinationGamma.h"

registerMooseObject("PhaseFieldApp", ACInterfaceInclinationGamma);

InputParameters
ACInterfaceInclinationGamma::validParams()
{
  InputParameters params = JvarMapKernelInterface<Kernel>::validParams();
  params.addClassDescription("Gradient energy Allen-Cahn Kernel for dgamma/dgradeta dependence.");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  // params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "The kappa used with the
  // kernel");
  params.addParam<bool>("variable_L",
                        false,
                        "The mobility is a function of any MOOSE variable (if "
                        "this is set to false L must be constant over the "
                        "entire domain!)");
  params.addParam<bool>("skip_off", false, "Skip the off-diagonal part, for testing/debugging.");
  params.addParam<bool>("debug_kernel",
                        false,
                        "If true, print extra mapping info and sanity warnings at construction.");
  params.addParam<MaterialPropertyName>("mask_name",
                                        "Name of a MaterialProperty to use as a mask.  "
                                        "If empty, mask = 1.0 everywhere.");
  params.addRequiredParam<std::string>("var_name_base", "specifies the base name of the variables");
  params.addRequiredParam<unsigned int>("op_num", "Total number of grain order parameters.");
  return params;
}

ACInterfaceInclinationGamma::ACInterfaceInclinationGamma(const InputParameters & parameters)
  : DerivativeMaterialInterface<JvarMapKernelInterface<Kernel>>(parameters),
    // Include the second order parts for if _variable_L
    _variable_L(getParam<bool>("variable_L")),
    _second_u(_variable_L ? &second() : nullptr),
    _second_test(_variable_L ? &secondTest() : nullptr),
    _second_phi(_variable_L ? &secondPhi() : nullptr),
    _op_num(getParam<unsigned int>("op_num")),
    // IJ pair information
    _var_name_base(getParam<std::string>("var_name_base")),
    _no_ij_pairs(getMaterialProperty<bool>("no_ij_pairs")), // elem_no_ij
    _op_has_active_pair(getMaterialProperty<std::vector<unsigned char>>("op_has_active_pair")),
    _ij_i(getMaterialProperty<std::vector<unsigned int>>("ij_i")),
    _ij_j(getMaterialProperty<std::vector<unsigned int>>("ij_j")),
    // mu or const_m
    _mu(getMaterialProperty<Real>("mu")),
    // Gamma
    _dgammaij_dgradeta(getMaterialProperty<std::vector<RealGradient>>("dgamma_dgradeta")),
    _d2gammaij_dgradeta2(getMaterialProperty<std::vector<RealTensorValue>>("d2gamma_dgradeta2")),
    // AC Mobility
    _L(getMaterialProperty<Real>("mob_name")),
    _skip_off(getParam<bool>("skip_off")),
    _debug_kernel(getParam<bool>("debug_kernel")),
    _mask(isParamValid("mask_name") ? &getMaterialProperty<Real>("mask_name") : nullptr),
    _mask_tf(isParamValid("mask_name")),
    // Variable L derivatives for u and gradeta
    _grain_slots(_n_args),
    _dLdu(_variable_L ? &getMaterialPropertyDerivative<Real>("mob_name", _var.name()) : nullptr),
    _d2Ldu2(_variable_L ? &getMaterialPropertyDerivative<Real>("mob_name", _var.name(), _var.name())
                        : nullptr),
    _dL_dgradu(),
    _dL_dgradeta(_n_args),
    _d2L_dgradudu(),
    _d2L_dgradudarg(_n_args),
    _d2L_dgradudgradeta(),
    _d2L_dgradetadarg(_n_args),
    _d2L_dgradetadu(_n_args),
    // for n_args
    _dLdarg(_n_args),
    _d2Ldargdu(_n_args),
    _d2Ldarg2(_n_args),
    // _d2Ldargdgrad_eta(_n_args),
    _gradarg(_n_args),
    _second_arg(_n_args)

{
  // Parse my OP index from my variable name
  _my_op = OpIndexUtils::parseOpIndex(_var.name(), _var_name_base);
  // _grain_slots.clear();
  // _grain_slots.reserve(_n_args);

  if (_variable_L)
  {
    _dL_dgradu = &getMaterialProperty<RealGradient>("dL_dgradeta_" + Moose::stringify(_my_op));
    _d2L_dgradudgradeta = &getMaterialProperty<std::vector<RealTensorValue>>(
        "d2L_dgradeta2_" + Moose::stringify(_my_op));
    _d2L_dgradudu = &getMaterialPropertyDerivative<RealGradient>(
        "dL_dgradeta_" + Moose::stringify(_my_op), _var.name());
  }

  // Pre-size using op_num — no need for max_op_seen
  _eta_by_op.assign(_op_num, nullptr);
  _op_to_jvar.assign(_op_num, -1);

  // Make sure my own slot is populated even if I'm not listed among coupled vars
  _eta_by_op[_my_op] = &_u; // use the Kernel's own VariableValue
  _op_to_jvar[_my_op] = _var.number();

  // Fill entries for coupled variables
  for (unsigned i = 0; i < _n_args; ++i)
  {
    MooseVariable * ivar = _coupled_standard_moose_vars[i];
    const VariableName iname = ivar->name();
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

    _grain_slots[i] = false; // default to not a grain

    if (MooseUtils::beginsWith(iname, _var_name_base) && (iname != _var.name()))
    {
      const unsigned k = OpIndexUtils::parseOpIndex(ivar->name(), _var_name_base);
      if (k == _my_op)
        mooseError("Arg variable ",
                   iname,
                   " has OP index ",
                   k,
                   " but the variable for this kernel is ",
                   _var.name());
      // Validate against op_num
      if (k >= _op_num)
        mooseError(
            "Coupled variable '", iname, "' has OP index ", k, " but op_num = ", _op_num, ".");
      // _grain_val_by_k.emplace(k, &coupledValue("coupled_variables", i));
      // _grain_val_by_argi.emplace(i, &coupledValue("coupled_variables", i));
      _eta_by_op[k] = &(ivar->sln());
      _op_to_jvar[k] = static_cast<int>(ivar->number());
      // _grain_slots.push_back(i);
      _grain_slots[i] = true; // this arg index is a grain
      if (_variable_L)
      {
        _dL_dgradeta[i] = &getMaterialProperty<RealGradient>("dL_dgradeta_" + Moose::stringify(k));
        _d2L_dgradetadu[i] = &getMaterialPropertyDerivative<RealGradient>(
            "dL_dgradeta_" + Moose::stringify(k), iname);
        _d2L_dgradetadarg[i].resize(_n_args);
        for (unsigned j = 0; j < _n_args; ++j)
        {
          _d2L_dgradetadarg[i][j] = &getMaterialPropertyDerivative<RealGradient>(
              "dL_dgradeta_" + Moose::stringify(k), _coupled_standard_moose_vars[j]->name());
        }
        // _gradarg[k] = &(ivar->gradSln());
        // _second_arg[k] = &(ivar->secondSln());
      }
    }

    // Moelans L
    if (_variable_L)
    {
      _gradarg[i] = &(ivar->gradSln());
      _second_arg[i] = &(ivar->secondSln());
      // L dependencies
      _dLdarg[i] = &getMaterialPropertyDerivative<Real>("mob_name", i);
      _d2Ldargdu[i] = &getMaterialPropertyDerivative<Real>("mob_name", iname, _var.name());
      _d2L_dgradudarg[i] = &getMaterialPropertyDerivative<RealGradient>(
          "dL_dgradeta_" + Moose::stringify(_my_op), iname);
      // _d2Ldargdgrad_eta[i] =
      //     &getMaterialPropertyDerivative<std::vector<RealGradient>>("dL_dgradeta", iname);
      _d2Ldarg2[i].resize(_n_args);
      for (unsigned int j = 0; j < _n_args; ++j)
        _d2Ldarg2[i][j] = &getMaterialPropertyDerivative<Real>("mob_name", i, j);
    }
  }

  for (unsigned k = 0; k < _op_num; ++k)
    if (!_eta_by_op[k])
      mooseError("ACInterfaceInclinationGamma for variable ",
                 _var.name(),
                 " is missing coupled value for OP index ",
                 k,
                 ". All grain OPs must be coupled so compact active pairs can be evaluated.");

  // debugging dump
  if (_debug_kernel && this->processor_id() == 0)
  {
    std::ostringstream oss;
    oss << "ACInterfaceInclinationGamma mapping (base='" << _var_name_base
        << "', op_num=" << _op_num << "):\n";
    oss << "  my_op=" << _my_op << "  var='" << _var.name() << "'  jvar=" << _var.number() << "\n";
    for (unsigned k = 0; k < _op_num; ++k)
    {
      oss << "  OP " << k << ": "
          << "val_ptr=" << (void *)_eta_by_op[k] << "  jvar=" << _op_to_jvar[k]
          << "  cvar=" << mapJvarToCvar(_op_to_jvar[k]);
      // Try to recover a name for clarity (best-effort)
      bool named = false;
      for (unsigned i = 0; i < _n_args && !named; ++i)
        if (OpIndexUtils::parseOpIndex(_coupled_standard_moose_vars[i]->name(), _var_name_base) ==
            k)
        {
          oss << "  name='" << _coupled_standard_moose_vars[i]->name() << "'";
          named = true;
        }
      if (!named && k == _my_op)
        oss << "  name='" << _var.name() << "'";
      oss << "\n";
    }
    for (unsigned i = 0; i < _n_args; ++i)
    {
      oss << "arg_i=" << i << "  grain_val_i=" << _grain_slots[i] << "\n";
    }
    mooseWarning(oss.str()); // use mooseInfo if you prefer non-yellow output
  }
}

void
ACInterfaceInclinationGamma::initialSetup()
{
  validateCoupling<Real>("mob_name");
}

RealGradient
ACInterfaceInclinationGamma::gradL()
{
  mooseAssert(_variable_L, "gradL() should only be called when variable_L = true.");
  mooseAssert(_second_u, "second_u pointer is null in gradL().");

  RealGradient g = _grad_u[_qp] * (*_dLdu)[_qp];
  g += (*_second_u)[_qp] * (*_dL_dgradu)[_qp];

  for (unsigned int i = 0; i < _n_args; ++i)
  {
    g += (*_gradarg[i])[_qp] * (*_dLdarg[i])[_qp];

    if (_grain_slots[i])
      g += (*_second_arg[i])[_qp] * (*_dL_dgradeta[i])[_qp];
  }

  return g;
}

RealGradient
ACInterfaceInclinationGamma::nablaLPsi() // RH $L \nabla \psi$
{
  // sum is the product rule gradient \f$ \nabla (L\psi) \f$
  RealGradient sum = _L[_qp] * _grad_test[_i][_qp];

  if (_variable_L)
    sum += gradL() * _test[_i][_qp];

  return sum;
}

Real
ACInterfaceInclinationGamma::computeQpResidual()
{
  // Exit shortcuts
  if (_no_ij_pairs[_qp])
    return 0.0;

  const auto & op_has_active_pair = _op_has_active_pair[_qp];

  mooseAssert(_my_op < op_has_active_pair.size(),
              "my_op index is out of range for op_has_active_pair.");

  if (!op_has_active_pair[_my_op])
    return 0.0;

  // const RealGradient nabla_Lpsi = nablaLPsi(); // L * grad(test) for now

  const auto & ii = _ij_i[_qp];
  const auto & jj = _ij_j[_qp];
  const auto & dgamma = _dgammaij_dgradeta[_qp];

  mooseAssert(ii.size() == jj.size() && ii.size() == dgamma.size(),
              "ij_i/ij_j/dgamma length mismatch at qp");

  Real sum = 0.0;
  // Added to try to speed up
  const RealGradient nabla_Lpsi = _variable_L ? nablaLPsi() : _L[_qp] * _grad_test[_i][_qp];

  for (std::size_t k = 0; k < ii.size(); ++k)
  {
    const unsigned i = ii[k];
    const unsigned j = jj[k];

    // skip if my OP isn't involved
    if (_my_op != i && _my_op != j)
      continue;

    // Determine which is the "other" OP in the pair and get its value
    const unsigned other_op = (_my_op == i) ? j : i;
    const Real eta_other = eta_at(other_op);
    // const Real eta_other = get_val_by_k(other_op);
    // const Real eta_fac = _u[_qp] * _u[_qp] * (eta_other * eta_other); // η_u^2 * η_other^2

    // sign from symmetry: dγ/d∇η_i = - dγ/d∇η_j
    const Real sgn = (_my_op == i) ? 1.0 : -1.0;

    // dot with ∇(L ψ)
    sum += sgn * (dgamma[k] * nabla_Lpsi) * _u[_qp] * _u[_qp] * (eta_other * eta_other);
  }

  if (_mask_tf)
    return (*_mask)[_qp] * _mu[_qp] * sum;
  else
    return _mu[_qp] * sum;
}

Real
ACInterfaceInclinationGamma::computeQpJacobian()
{
  // Exit shortcuts
  if (_no_ij_pairs[_qp])
    return 0.0;

  const auto & op_has_active_pair = _op_has_active_pair[_qp];

  mooseAssert(_my_op < op_has_active_pair.size(),
              "my_op index is out of range for op_has_active_pair.");

  if (!op_has_active_pair[_my_op])
    return 0.0;

  const auto & ii = _ij_i[_qp];
  const auto & jj = _ij_j[_qp];
  const auto & dgamma = _dgammaij_dgradeta[_qp];
  const auto & d2gamma = _d2gammaij_dgradeta2[_qp];

  mooseAssert(ii.size() == jj.size(), "ij_i/ij_j length mismatch at qp");
  mooseAssert(ii.size() == dgamma.size(), "ij_i/dgamma length mismatch at qp");
  mooseAssert(ii.size() == d2gamma.size(), "ij_i/d2gamma length mismatch at qp");

  Real ddir = 0.0; // derivative wrt u
  Real dind = 0.0; // derivative wrt grad_u

  // Added to try to speed up
  const RealGradient nabla_Lpsi = _variable_L ? nablaLPsi() : _L[_qp] * _grad_test[_i][_qp];

  for (std::size_t k = 0; k < ii.size(); ++k)
  {
    const unsigned i = ii[k], j = jj[k];
    // if (_my_op != i && _my_op != j)
    //   continue;

    // // Determine which is the "other" OP in the pair and get its value
    // const Real sgn = (_my_op == i) ? 1.0 : -1.0;
    // const unsigned other_op = (_my_op == i) ? j : i;
    // Alternate faster?:
    const bool my_is_i = (_my_op == i);
    if (!my_is_i && _my_op != j)
      continue;

    const Real sgn = my_is_i ? 1.0 : -1.0;
    const unsigned other_op = my_is_i ? j : i;
    const Real eta_other = eta_at(other_op);
    // const Real eta_other = get_val_by_k(other_op);

    ddir += 2 * _u[_qp] * sgn * dgamma[k] * eta_other * eta_other * _phi[_j][_qp] * nabla_Lpsi;
    dind +=
        _u[_qp] * _u[_qp] * d2gamma[k] * eta_other * eta_other * _grad_phi[_j][_qp] * nabla_Lpsi;
    // d2gamma is ii or jj here so both +

    // Variable L Moelans
    if (_variable_L)
    {
      // Grad L partials
      static const RealTensorValue I(1, 0, 0, 0, 1, 0, 0, 0, 1);
      // The direct u pieces
      RealGradient dgradLdu =
          (*_d2Ldu2)[_qp] * _grad_u[_qp] + (*_d2L_dgradudu)[_qp] * (*_second_u)[_qp];
      RealTensorValue dgradLdgradu = libMesh::outer_product((*_d2L_dgradudu)[_qp], _grad_u[_qp]) +
                                     I * (*_dLdu)[_qp] +
                                     (*_d2L_dgradudgradeta)[_qp][_my_op] * (*_second_u)[_qp];
      // The cross terms with eta/gradeta dependence in grad L
      for (unsigned int i = 0; i < _n_args; ++i)
      {
        dgradLdu += (*_d2Ldargdu[i])[_qp] * (*_gradarg[i])[_qp];
        if (_grain_slots[i])
          dgradLdu += (*_d2L_dgradetadu[i])[_qp] * (*_second_arg[i])[_qp];
        dgradLdgradu += libMesh::outer_product((*_d2L_dgradudarg[i])[_qp], (*_gradarg[i])[_qp]) -
                        (*_d2L_dgradudgradeta)[_qp][_my_op] * (*_second_arg[i])[_qp];
      }

      // Direct L dependence
      ddir += (*_dLdu)[_qp] * _u[_qp] * _u[_qp] * sgn * dgamma[k] * eta_other * eta_other *
              _phi[_j][_qp] * _grad_test[_i][_qp];
      // Direct grad L dependence
      ddir += dgradLdu * _u[_qp] * _u[_qp] * sgn * dgamma[k] * eta_other * eta_other *
              _phi[_j][_qp] * _test[_i][_qp];
      // Indirect L dependence (of grad u)
      dind += (*_dL_dgradu)[_qp] * _u[_qp] * _u[_qp] * sgn * dgamma[k] * eta_other * eta_other *
              _grad_phi[_j][_qp] * _grad_test[_i][_qp];
      // Indirect grad L dependence (of grad u)
      dind += dgradLdgradu * _u[_qp] * _u[_qp] * sgn * dgamma[k] * eta_other * eta_other *
              _grad_phi[_j][_qp] * _test[_i][_qp];
    }
  }

  if (_mask_tf)
    return (*_mask)[_qp] * _mu[_qp] * (ddir + dind);
  else
    return _mu[_qp] * (ddir + dind);
}

Real
ACInterfaceInclinationGamma::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_skip_off || _no_ij_pairs[_qp])
    return 0.0;

  int v_int = -1;
  for (unsigned op = 0; op < _op_to_jvar.size(); ++op)
    if (_op_to_jvar[op] == static_cast<int>(jvar))
    {
      v_int = static_cast<int>(op);
      break;
    }

  if (v_int < 0)
    return 0.0;

  const unsigned int v = static_cast<unsigned int>(v_int);

  if (v == _my_op)
    return 0.0;

  const auto & op_has_active_pair = _op_has_active_pair[_qp];

  if (_my_op >= op_has_active_pair.size() || !op_has_active_pair[_my_op])
    return 0.0;

  if (v >= op_has_active_pair.size() || !op_has_active_pair[v])
    return 0.0;

  const auto & ii = _ij_i[_qp];
  const auto & jj = _ij_j[_qp];
  const auto & dgamma = _dgammaij_dgradeta[_qp];
  const auto & d2gamma = _d2gammaij_dgradeta2[_qp];

  mooseAssert(ii.size() == jj.size(), "ij_i/ij_j length mismatch at qp");
  mooseAssert(ii.size() == dgamma.size(), "ij_i/dgamma length mismatch at qp");
  mooseAssert(ii.size() == d2gamma.size(), "ij_i/d2gamma length mismatch at qp");

  const Real u = _u[_qp];
  const Real u2 = u * u;
  const Real phi = _phi[_j][_qp];

  const RealGradient & grad_phi = _grad_phi[_j][_qp];
  const RealGradient & grad_test = _grad_test[_i][_qp];

  const RealGradient nabla_Lpsi = _variable_L ? nablaLPsi() : _L[_qp] * grad_test;

  Real ddir = 0.0;
  Real dind = 0.0;

  for (std::size_t k = 0; k < ii.size(); ++k)
  {
    const unsigned int i = ii[k];
    const unsigned int j = jj[k];

    const bool my_is_i = (_my_op == i);
    const bool my_is_j = (_my_op == j);

    if (!my_is_i && !my_is_j)
      continue;

    const unsigned int other_op = my_is_i ? j : i;

    if (other_op != v)
      continue;

    const Real sgn = my_is_i ? 1.0 : -1.0;
    const Real eta_other = eta_at(other_op);
    const Real eta_other2 = eta_other * eta_other;

    ddir += 2.0 * eta_other * u2 * sgn * phi * (dgamma[k] * nabla_Lpsi);

    dind += u2 * eta_other2 * ((-d2gamma[k] * grad_phi) * nabla_Lpsi);

    if (_variable_L)
    {
      const unsigned int cvar = mapJvarToCvar(jvar);

      static const RealTensorValue I(1, 0, 0, 0, 1, 0, 0, 0, 1);

      RealGradient dgradLdarg = (*_d2Ldargdu[cvar])[_qp] * _grad_u[_qp] +
                                (*_d2L_dgradudarg[cvar])[_qp] * (*_second_u)[_qp];

      RealTensorValue dgradLdgradarg(0.0);

      if (_grain_slots[cvar])
      {
        dgradLdgradarg += libMesh::outer_product((*_d2L_dgradetadu[cvar])[_qp], _grad_u[_qp]) +
                          (*_d2L_dgradudgradeta)[_qp][other_op] * (*_second_u)[_qp] +
                          I * (*_dLdarg[cvar])[_qp];
      }

      for (unsigned int a = 0; a < _n_args; ++a)
      {
        dgradLdarg += (*_d2Ldarg2[cvar][a])[_qp] * (*_gradarg[a])[_qp];

        if (_grain_slots[a])
          dgradLdarg += (*_d2L_dgradetadarg[a][cvar])[_qp] * (*_second_arg[a])[_qp];

        if (_grain_slots[cvar])
          dgradLdgradarg +=
              libMesh::outer_product((*_d2L_dgradetadarg[cvar][a])[_qp], (*_gradarg[a])[_qp]);
      }

      const Real eta_factor = eta_other2;
      const Real common = u2 * sgn * eta_factor;

      ddir += (*_dLdarg[cvar])[_qp] * common * (dgamma[k] * grad_test) * phi;

      ddir += (dgradLdarg * dgamma[k]) * common * phi * _test[_i][_qp];

      if (_grain_slots[cvar])
        dind += ((*_dL_dgradeta[cvar])[_qp] * dgamma[k]) * common * (grad_phi * grad_test);

      dind += ((dgradLdgradarg * dgamma[k]) * grad_phi) * common * _test[_i][_qp];
    }

    // only one possible matching compact pair
    break;
  }

  const Real prefactor = (_mask_tf ? (*_mask)[_qp] : 1.0) * _mu[_qp];

  return prefactor * (ddir + dind);
}
