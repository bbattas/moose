/****************************************************************/
/*                  DO NOT MODIFY THIS HEADER                   */
/*                           Marmot                             */
/*                                                              */
/*            (c) 2017 Battelle Energy Alliance, LLC            */
/*                     ALL RIGHTS RESERVED                      */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*             Under Contract No. DE-AC07-05ID14517             */
/*             With the U. S. Department of Energy              */
/*                                                              */
/*             See COPYRIGHT for full restrictions              */
/****************************************************************/

#include "GBSegregation.h"
#include "libmesh/quadrature.h"

registerMooseObject("MarmotApp", GBSegregation);

template<>
InputParameters validParams<GBSegregation>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Grain boundary concentration segregation");
  params.addRequiredCoupledVarWithAutoBuild("etas", "var_name_base", "op_num", "Array of order parameters that define the solid phase");
  params.addParam<Real>("c_prefactor", 1.0, "Arrhenius equation prefactor for the concentration");
  params.addRequiredParam<Real>("formation_energy", "Formation energy for defect");
  params.addRequiredParam<Real>("GB_energy", "Reduction of formation energy on grain boundaries");
  params.addParam<Real>("max_bulk_concentration", 0.01, "Maximum allowable value of bulk defect concentration");
  params.addParam<Real>("max_GB_concentration", 0.1, "Maximum allowable value of GB defect concentration");
  params.addRequiredCoupledVar("Temperature", "Temperature variable in Kelvin");
  params.addParam<std::string>("c_name", "c_eq", "Name of concentration being calculated");
  return params;
}

GBSegregation::GBSegregation(const InputParameters & parameters) : DerivativeMaterialInterface<Material>(parameters),
  _A(getParam<Real>("c_prefactor")),
  _Ef(getParam<Real>("formation_energy")),
  _Egb(getParam<Real>("GB_energy")),
  _max_bulk(getParam<Real>("max_bulk_concentration")),
  _max_gb(getParam<Real>("max_GB_concentration")),
  _T(coupledValue("Temperature")),
  _nop(coupledComponents("etas")),
  _ops(_nop),
  _op_names(_nop),
  _ceq_name(getParam<std::string>("c_name")),
  _ceq(declareProperty<Real>(_ceq_name)),
  _dceqdopa(_nop),
  _d2ceqdopa2(_nop),
  _d2ceqdopadopb(_nop)
{
  for (auto i = 0; i < _nop; ++i)
  {
    _ops[i] = &coupledValue("etas", i);
    _op_names[i] = getVar("etas", i)->name();
    _dceqdopa[i] = &declarePropertyDerivative<Real>(_ceq_name, _op_names[i]);
    _d2ceqdopa2[i] = &declarePropertyDerivative<Real>(_ceq_name, _op_names[i], _op_names[i]);
    _d2ceqdopadopb[i].resize(i);
    for (auto j = 0; j < i; ++j)
      _d2ceqdopadopb[i][j] = &declarePropertyDerivative<Real>(_ceq_name, _op_names[i], _op_names[j]);
  }
}

void
GBSegregation::computeQpProperties()
{
  Real kB = 8.617343e-5;

  Real sum_op = 0.0;
  for (auto i = 0; i < _nop; ++i)
    sum_op += (*_ops[i])[_qp] * (*_ops[i])[_qp];

  Real cbeq = _A * std::exp(-_Ef / (kB * _T[_qp]));
  if (cbeq > _max_bulk)
    cbeq = _max_bulk;
  Real cgbeq = _A * std::exp(-(_Ef - _Egb) / (kB * _T[_qp]));
  if (cgbeq > _max_gb)
    cgbeq = _max_gb;

  _ceq[_qp] = cbeq + (cgbeq - cbeq) * 4.0 * (1.0 - sum_op) * (1.0 - sum_op);
  for (auto i = 0; i < _nop; ++i)
  {
    (*_dceqdopa[i])[_qp] = -16.0 * (*_ops[i])[_qp] * (cgbeq - cbeq) * (1.0 - sum_op);
    (*_d2ceqdopa2[i])[_qp] = -16.0 * (cgbeq - cbeq) * (1.0 - sum_op - 2.0 * (*_ops[i])[_qp] * (*_ops[i])[_qp]);

    for (auto j = 0; j < i; ++j)
      (*_d2ceqdopadopb[i][j])[_qp] = 32.0 * (*_ops[i])[_qp] * (*_ops[j])[_qp] * (cgbeq - cbeq);
  }
}
