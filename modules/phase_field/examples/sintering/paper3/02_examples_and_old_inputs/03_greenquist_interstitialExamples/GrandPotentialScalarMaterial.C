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

#include "GrandPotentialScalarMaterial.h"
#include "libmesh/quadrature.h"

registerMooseObject("MarmotApp", GrandPotentialScalarMaterial);

template <>
InputParameters
validParams<GrandPotentialScalarMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Generates non-uniform mobility scalar with no directional dependence");
  params.addRequiredCoupledVarWithAutoBuild("etas", "var_name_base", "op_num", "Array of coupled order parameters");
  params.addRequiredCoupledVar("void_variable", "Name of variable that defines the void phase");
  params.addRequiredParam<NonlinearVariableName>("chemical_potential", "Name of chemical potential variable for this species");
  params.addRequiredCoupledVar("Temperature", "Temperature variable name");
  params.addRequiredParam<Real>("D0", "Diffusivity prefactor");
  params.addRequiredParam<Real>("migration_energy", "Migration energy for species");
  params.addParam<Real>("bulk_index", 1.0, "Diffusivity multiplier weight for the bulk region");
  params.addParam<Real>("GB_index", 1.0, "Diffusivity multiplier weight for the GB region");
  params.addParam<Real>("surface_index", 1.0, "Diffusivity multiplier weight for the surface region");
  params.addRequiredParam<MaterialPropertyName>("susceptibility", "Susceptibility material property name");
  params.addParam<std::string>("D_name", "Name of diffusivity scalar to output");
  params.addParam<std::string>("M_name", "Name of mobility scalar to output");
  return params;
}

GrandPotentialScalarMaterial::GrandPotentialScalarMaterial(const InputParameters & parameters) : DerivativeMaterialInterface<Material>(parameters),
  _nop(coupledComponents("etas")),
  _ops(_nop),
  _op_names(_nop),
  _phi(coupledValue("void_variable")),
  _phi_name(getVar("void_variable", 0)->name()),
  _w_name(getParam<NonlinearVariableName>("chemical_potential")),
  _T(coupledValue("Temperature")),
  _D0(getParam<Real>("D0")),
  _Em(getParam<Real>("migration_energy")),
  _wB(getParam<Real>("bulk_index")),
  _wGB(getParam<Real>("GB_index")),
  _wS(getParam<Real>("surface_index")),
  _chi_name(getParam<MaterialPropertyName>("susceptibility")),
  _chi(getMaterialProperty<Real>(_chi_name)),
  _dchidw(getMaterialPropertyDerivativeByName<Real>(_chi_name, _w_name)),
  _dchidphi(getMaterialPropertyDerivativeByName<Real>(_chi_name, _phi_name)),
  _d2chidw2(getMaterialPropertyDerivativeByName<Real>(_chi_name, _w_name, _w_name)),
  _d2chidphidw(getMaterialPropertyDerivativeByName<Real>(_chi_name, _phi_name, _w_name)),
  _d2chidphi2(getMaterialPropertyDerivativeByName<Real>(_chi_name, _phi_name, _phi_name)),
  _Dname(getParam<std::string>("D_name")),
  _D(declareProperty<Real>(_Dname)),
  _dDdphi(declarePropertyDerivative<Real>(_Dname, _phi_name)),
  _dDdeta(_nop),
  _d2Ddphi2(declarePropertyDerivative<Real>(_Dname, _phi_name, _phi_name)),
  _d2Ddeta2(_nop),
  _Mname(getParam<std::string>("M_name")),
  _M(declareProperty<Real>(_Mname)),
  _dMdw(declarePropertyDerivative<Real>(_Mname, _w_name)),
  _dMdphi(declarePropertyDerivative<Real>(_Mname, _phi_name)),
  _dMdeta(_nop),
  _d2Mdw2(declarePropertyDerivative<Real>(_Mname, _w_name, _w_name)),
  _d2Mdphidw(declarePropertyDerivative<Real>(_Mname, _phi_name, _w_name)),
  _d2Mdetadw(_nop),
  _d2Mdphi2(declarePropertyDerivative<Real>(_Mname, _phi_name, _phi_name)),
  _d2Mdetadphi(_nop),
  _d2Mdeta2(_nop)
{
  for (unsigned int i = 0; i < _nop; ++i)
  {
    _ops[i] = &coupledValue("etas", i);
    _op_names[i] = getVar("etas", i)->name();

    _dDdeta[i] = &declarePropertyDerivative<Real>(_Dname, _op_names[i]);
    _d2Ddeta2[i].resize(i);
    _dMdeta[i] = &declarePropertyDerivative<Real>(_Mname, _op_names[i]);
    _d2Mdetadw[i] = &declarePropertyDerivative<Real>(_Mname, _op_names[i], _w_name);
    _d2Mdetadphi[i] = &declarePropertyDerivative<Real>(_Mname, _op_names[i], _phi_name);
    _d2Mdeta2[i].resize(i);

    for (unsigned int j = 0; j < i; ++j)
    {
      _d2Ddeta2[i][j] = &declarePropertyDerivative<Real>(_Dname, _op_names[i], _op_names[j]);
      _d2Mdeta2[i][j] = &declarePropertyDerivative<Real>(_Mname, _op_names[i], _op_names[j]);
    }
  }
}

void
GrandPotentialScalarMaterial::computeQpProperties()
{
  Real kB = 8.617343e-5; // eV/K
  Real D_bulk = _D0 * std::exp(-_Em / kB / _T[_qp]);

  Real gb_region = 0;
  for (unsigned int i = 0; i < _nop; ++i)
    for (unsigned int j = 0; j < i; ++j)
      gb_region += 2 * (*_ops[i])[_qp] * (*_ops[j])[_qp];

  Real surf_region = _phi[_qp] * _phi[_qp] * (1.0 - _phi[_qp]) * (1.0 - _phi[_qp]);

  _D[_qp] = D_bulk * (_wB + _wGB * gb_region + _wS * surf_region);
  _dDdphi[_qp] = D_bulk * _wS * 2.0 * _phi[_qp] * (1.0 - 3.0 * _phi[_qp] + 2.0 * _phi[_qp] * _phi[_qp]);
  _d2Ddphi2[_qp] = 2.0 * D_bulk * _wS * (1.0 - 6.0 * _phi[_qp] + 6.0 * _phi[_qp] * _phi[_qp]);

  _M[_qp] = _chi[_qp] * _D[_qp];
  _dMdw[_qp] = _dchidw[_qp] * _D[_qp];
  _dMdphi[_qp] = _dchidphi[_qp] * _D[_qp] + _chi[_qp] * _dDdphi[_qp];
  _d2Mdw2[_qp] = _d2chidw2[_qp] * _D[_qp];
  _d2Mdphidw[_qp] = _d2chidphidw[_qp] * _D[_qp] + _dchidw[_qp] * _dDdphi[_qp];
  _d2Mdphi2[_qp] = _d2chidphi2[_qp] * _D[_qp] + 2.0 * _dchidphi[_qp] * _dDdphi[_qp] + _chi[_qp] * _d2Ddphi2[_qp];

  for (unsigned int i = 0; i < _nop; ++i)
  {
    Real sum_not_i = 0.0;
    for (unsigned int j = 0; j < _nop; ++j)
      if (j != i)
        sum_not_i += (*_ops[j])[_qp];

    (*_dDdeta[i])[_qp] = D_bulk * _wGB * sum_not_i;
    (*_dMdeta[i])[_qp] = _chi[_qp] * (*_dDdeta[i])[_qp];
    (*_d2Mdetadw[i])[_qp] = _dchidw[_qp] * (*_dDdeta[i])[_qp];
    (*_d2Mdetadphi[i])[_qp] = _dchidphi[_qp] * (*_dDdeta[i])[_qp];

    for (unsigned int j = 0; j < i; ++j)
    {
      (*_d2Ddeta2[i][j])[_qp] = D_bulk * _wGB;
      (*_d2Mdeta2[i][j])[_qp] = _chi[_qp] * (*_d2Ddeta2[i][j])[_qp];
    }
  }
}
