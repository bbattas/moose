#include "ACInterfaceInclinationGamma.h"

registerMooseObject("PhaseFieldApp", ACInterfaceInclinationGamma);

InputParameters
ACInterfaceInclinationGamma::validParams()
{
  InputParameters params = JvarMapKernelInterface<Kernel>::validParams();
  params.addClassDescription("Gradient energy Allen-Cahn Kernel");
  params.addParam<MaterialPropertyName>("mob_name", "L", "The mobility used with the kernel");
  // params.addParam<MaterialPropertyName>("kappa_name", "kappa_op", "The kappa used with the
  // kernel");
  params.addParam<bool>("variable_L",
                        true,
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
    _second_u(second()),
    _second_test(secondTest()),
    _second_phi(secondPhi()),
    _op_num(getParam<unsigned int>("op_num")),
    // IJ pair information
    _var_name_base(getParam<std::string>("var_name_base")),
    _no_ij_pairs(getMaterialProperty<bool>("elem_no_ij")), // no_ij_pairs
    _ij_i(getMaterialProperty<std::vector<unsigned int>>("ij_i")),
    _ij_j(getMaterialProperty<std::vector<unsigned int>>("ij_j")),
    // mu or const_m
    _mu(getMaterialProperty<Real>("mu")),
    // Gamma
    _dgammaij_dgradeta(getMaterialProperty<std::vector<RealGradient>>("dgamma_dgradeta")),
    _d2gammaij_dgradeta2(getMaterialProperty<std::vector<RealTensorValue>>("d2gamma_dgradeta2")),
    // AC Mobility
    _L(getMaterialProperty<Real>("mob_name")),
    _variable_L(getParam<bool>("variable_L")),
    _skip_off(getParam<bool>("skip_off")),
    _debug_kernel(getParam<bool>("debug_kernel")),
    _mask(isParamValid("mask_name") ? &getMaterialProperty<Real>("mask_name") : nullptr),
    _mask_tf(isParamValid("mask_name")),
    _gradarg(_n_args),
    // other approach
    // _grain_args(),
    // _grain_idx(),
    // _grain_vals(),
    // _grain_k_test(),
    _grain_val_by_k(),
    _grain_val_by_argi()
{
  mooseWarning("Coupled args in agg kernel: ", _n_args);
  // Parse my OP index from my variable name
  _my_op = OpIndexUtils::parseOpIndex(_var.name(), _var_name_base);

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
      _grain_val_by_k.emplace(k, &coupledValue("coupled_variables", i));
      _grain_val_by_argi.emplace(i, &coupledValue("coupled_variables", i));
      //
      _eta_by_op[k] = &(ivar->sln());
      _op_to_jvar[k] = static_cast<int>(ivar->number());
    }

    // _eta_by_op[k] = &(ivar->sln());
    // _op_to_jvar[k] = static_cast<int>(ivar->number());
    _gradarg[i] = &(ivar->gradSln()); // your existing gradient cache
  }
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
          << "val_ptr=" << (void *)_eta_by_op[k] << "  jvar=" << _op_to_jvar[k];
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
    for (unsigned i = 0; i < _op_num; ++i)
    {
      oss << "k=" << i << "  grain_val_k=" << get_val_by_k(i) << "\n";
    }
    for (unsigned i = 0; i < _n_args; ++i)
    {
      oss << "arg_i=" << i << "  grain_val_i=" << get_val_by_argi(i) << "\n";
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
ACInterfaceInclinationGamma::nablaLPsi() // RH $L \nabla \psi$
{
  // sum is the product rule gradient \f$ \nabla (L\psi) \f$
  RealGradient sum = _L[_qp] * _grad_test[_i][_qp];

  // if (_variable_L)
  //   sum += gradL() * _test[_i][_qp];

  return sum;
}

Real
ACInterfaceInclinationGamma::computeQpResidual()
{
  if (_no_ij_pairs[_qp])
    return 0.0;

  // const RealGradient nabla_Lpsi = nablaLPsi(); // L * grad(test) for now

  const auto & ii = _ij_i[_qp];
  const auto & jj = _ij_j[_qp];
  const auto & dgamma = _dgammaij_dgradeta[_qp];

  mooseAssert(ii.size() == jj.size() && ii.size() == dgamma.size(),
              "ij_i/ij_j/dgamma length mismatch at qp");

  Real sum = 0.0;

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
    sum += sgn * (dgamma[k] * nablaLPsi()) * _u[_qp] * _u[_qp] * (eta_other * eta_other);
  }

  if (_mask_tf)
    return (*_mask)[_qp] * _mu[_qp] * sum;
  else
    return _mu[_qp] * sum;
}

Real
ACInterfaceInclinationGamma::computeQpJacobian()
{
  if (_no_ij_pairs[_qp])
    return 0.0;

  const auto & ii = _ij_i[_qp];
  const auto & jj = _ij_j[_qp];
  const auto & dgamma = _dgammaij_dgradeta[_qp];
  const auto & d2gamma = _d2gammaij_dgradeta2[_qp];

  Real ddir = 0.0; // derivative wrt u
  Real dind = 0.0; // derivative wrt grad_u

  for (std::size_t k = 0; k < ii.size(); ++k)
  {
    const unsigned i = ii[k], j = jj[k];
    if (_my_op != i && _my_op != j)
      continue;

    // const Real eta_i = eta_at(i);
    // const Real eta_j = eta_at(j);
    // const Real eta_fac = eta_i * eta_i * eta_j * eta_j;
    // const Real sgn = (_my_op == i) ? 1.0 : -1.0;

    // // (1) factor derivative part: (dγ·∇(Lψ)) * d/d(u)[η_i^2 η_j^2] * φ_j
    // const Real other = (_my_op == i) ? eta_j : eta_i;
    // const Real d_fac_du = 2.0 * eta_at(_my_op) * other * other;
    // J += sgn * (dgamma[k] * nabla_Lpsi) * d_fac_du * _phi[_j][_qp];

    // // (2) dγ term derivative: ( (T · ∇φ_j) · ∇(Lψ) ) * η_i^2 η_j^2
    // // Here d2gamma is ∂²γ/∂(∇η_i)². Using your symmetry, this is also
    // // the diagonal second derivative for u; the sign 'sgn' already picks i vs j.
    // const RealGradient T_times_gradphi = d2gamma[k] * _grad_phi[_j][_qp];
    // J += sgn * (T_times_gradphi * nabla_Lpsi) * eta_fac;

    // OR my way kinda
    // Determine which is the "other" OP in the pair and get its value
    const Real sgn = (_my_op == i) ? 1.0 : -1.0;
    const unsigned other_op = (_my_op == i) ? j : i;
    const Real eta_other = eta_at(other_op);
    // const Real eta_other = get_val_by_k(other_op);

    ddir += 2 * _u[_qp] * sgn * dgamma[k] * eta_other * eta_other * _phi[_j][_qp] * nablaLPsi();
    dind +=
        _u[_qp] * _u[_qp] * d2gamma[k] * eta_other * eta_other * _grad_phi[_j][_qp] * nablaLPsi();
  } // d2gamma is ii or jj here so both +

  if (_variable_L)
  {
    // add L dependencies
  }

  if (_mask_tf)
    return (*_mask)[_qp] * _mu[_qp] * (ddir + dind);
  else
    return _mu[_qp] * (ddir + dind);
}

Real
ACInterfaceInclinationGamma::computeQpOffDiagJacobian(unsigned int jvar)
{
  if ((_skip_off) || (_no_ij_pairs[_qp]))
  {
    return 0.0;
  }
  else
  {
    // Which OP index does this jvar correspond to?
    // If not a known OP, bail out early.
    int v = -1;
    for (unsigned op = 0; op < _op_to_jvar.size(); ++op)
      if (_op_to_jvar[op] == static_cast<int>(jvar))
      {
        v = static_cast<int>(op);
        break;
      }
    if (v < 0)
      return 0.0;

    // shorthand
    const auto & ii = _ij_i[_qp];
    const auto & jj = _ij_j[_qp];
    const auto & dgamma = _dgammaij_dgradeta[_qp];
    const auto & d2gamma = _d2gammaij_dgradeta2[_qp];

    // Real J = 0.0;
    Real ddir = 0.0; // derivative wrt u
    Real dind = 0.0; // derivative wrt grad_u

    for (std::size_t k = 0; k < ii.size(); ++k)
    {
      const unsigned i = ii[k], j = jj[k];

      // We need pairs that include both u and v
      const bool contains_u = (_my_op == i) || (_my_op == j);
      const bool contains_v = (static_cast<unsigned>(v) == i) || (static_cast<unsigned>(v) == j);
      if (!(contains_u && contains_v))
        continue;

      // const Real eta_i = eta_at(i), eta_j = eta_at(j);
      // const Real eta_fac = eta_i * eta_i * eta_j * eta_j;

      // const bool u_is_i = (_my_op == i);
      // const unsigned other_of_u = u_is_i ? j : i;

      // sign for dγ wrt ∇η_u:
      // const Real sgn_u = u_is_i ? 1.0 : -1.0;

      const Real sgn = (_my_op == i) ? 1.0 : -1.0;
      const unsigned other_op = (_my_op == i) ? j : i;
      const Real eta_other = eta_at(other_op);

      // (1) factor derivative w.r.t. eta_v
      // d/d(eta_v) [eta_i^2 eta_j^2] = 2*eta_v*(other)^2
      // const Real other = (static_cast<unsigned>(v) == i) ? eta_j : eta_i;
      // const Real d_fac_dv = 2.0 * eta_at(static_cast<unsigned>(v)) * other * other;
      // J += sgn_u * (dgamma[k] * nablaLPsi()) * d_fac_dv * _phi[_j][_qp];

      // (2) mixed second derivative piece:
      // d/d(∇eta_v) [∂γ/∂∇eta_u] : ∇φ_j
      // Using your stated symmetry: mixed = - (diagonal wrt u)
      // const RealGradient mixed_times_gradphi = -(d2gamma[k] * _grad_phi[_j][_qp]);
      // J += sgn_u * (mixed_times_gradphi * nablaLPsi()) * eta_fac;

      // Direct dependence on arg
      ddir += 2 * eta_other * _u[_qp] * _u[_qp] * sgn * dgamma[k] * _phi[_j][_qp] * nablaLPsi();
      // Indirect dependence on grad_arg
      dind += _u[_qp] * _u[_qp] * (-d2gamma[k] * _grad_phi[_j][_qp]) * eta_other * eta_other *
              nablaLPsi();
    }

    if (_mask_tf)
      return (*_mask)[_qp] * _mu[_qp] * (ddir + dind);
    else
      return _mu[_qp] * (ddir + dind);
  }
}
