//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GrandPotentialMultiSinteringMaterial.h"
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

registerMooseObject("PhaseFieldApp", GrandPotentialMultiSinteringMaterial);

InputParameters
GrandPotentialMultiSinteringMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Material to calculate switching and thermodynamic properties for the "
      "grand potential sintering model when using both vacancies and interstitials, so two mu "
      "values. Hard coded to only work for 2 mu and the parabolic solution.");
  params.addRequiredCoupledVarWithAutoBuild(
      "etas", "var_name_base", "op_num", "Array of order parameters that describe solid phase");
  params.addRequiredCoupledVar(
      "chemical_potentials",
      "The name of the vacancy and interstitial chemical potential variables");
  // params.addRequiredCoupledVar("chemical_potential_int",
  //                              "The name of the interstitial chemical potential variable");
  params.addRequiredCoupledVar("void_op", "The name of the void phase order parameter");
  params.addRequiredCoupledVar("Temperature", "Name of the temperature variable with units of K");
  // Parabolic coefficients
  params.addParam<MaterialPropertyName>("vac_solid_energy_coefficient",
                                        1.0,
                                        "Parabolic vacancy solid energy coefficient "
                                        "(energy/volume). Only used for parabolic energy.");
  params.addRequiredParam<MaterialPropertyName>(
      "vac_void_energy_coefficient", "Parabolic vacancy void energy coefficient (energy/volume)");
  params.addParam<MaterialPropertyName>("int_solid_energy_coefficient",
                                        1.0,
                                        "Parabolic interstitial solid energy coefficient "
                                        "(energy/volume). Only used for parabolic energy.");
  params.addRequiredParam<MaterialPropertyName>(
      "int_void_energy_coefficient",
      "Parabolic interstitial void energy coefficient (energy/volume)");
  // GB/Surface energy
  params.addRequiredParam<MaterialPropertyName>(
      "surface_energy", "Material of surface energy in units of problem (energy/area)");
  params.addRequiredParam<MaterialPropertyName>(
      "grainboundary_energy",
      "Material of grain boundary energy in units of problem (energy/area)");
  params.addParam<Real>("int_width", 1, "Interface width in units of problem (length)");
  params.addParam<Real>("surface_switch_value",
                        0.3,
                        "Value between 0 and 1 that determines when the interface begins to switch "
                        "from surface to GB. Small values give less error while large values "
                        "converge better.");
  // Defect equilibrium concentrations
  params.addRequiredParam<MaterialPropertyName>(
      "equilibrium_vacancy_concentration",
      "Name of material that determines the equilibrium vacancy concentration in the solid phase");
  params.addRequiredParam<MaterialPropertyName>("equilibrium_interstitial_concentration",
                                                "Name of material that determines the equilibrium "
                                                "interstitial concentration in the solid phase");
  params.addParam<Real>("atomic_volume", 0.04092, "Atomic volume of material");
  MooseEnum solid_energy_model("PARABOLIC DILUTE IDEAL", "PARABOLIC");
  params.addParam<MooseEnum>("solid_energy_model",
                             solid_energy_model,
                             "Type of energy function to use for the solid phase.");
  params.addParam<bool>(
      "mass_conservation", false, "imposing strict mass conservation formulation");
  return params;
}
// Will be using u for vacancies and i for interstitials since v and s are already void and solid
GrandPotentialMultiSinteringMaterial::GrandPotentialMultiSinteringMaterial(
    const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _neta(coupledComponents("etas")),
    _eta(_neta),
    _eta_name(_neta),
    // chemical potentials
    _wu(coupledValue("chemical_potential_vac")),
    _wu_name(coupledName("chemical_potential_vac", 0)),
    _wi(coupledValue("chemical_potential_int")),
    _wi_name(coupledName("chemical_potential_int", 0)),
    _phi(coupledValue("void_op")),
    _phi_name(coupledName("void_op", 0)),
    // Defect concentrations
    _csu_eq_name(getParam<MaterialPropertyName>("equilibrium_vacancy_concentration")),
    _csu_eq(getMaterialProperty<Real>(_csu_eq_name)),
    _dcsu_eq(_neta),
    _d2csu_eq(_neta),
    _csi_eq_name(getParam<MaterialPropertyName>("equilibrium_interstitial_concentration")),
    _csi_eq(getMaterialProperty<Real>(_csi_eq_name)),
    _dcsi_eq(_neta),
    _d2csi_eq(_neta),
    _T(coupledValue("Temperature")),
    // parabolic coeffs
    _kvu(getMaterialProperty<Real>("vac_void_energy_coefficient")),
    _ksu(getMaterialProperty<Real>("vac_solid_energy_coefficient")),
    _kvi(getMaterialProperty<Real>("int_void_energy_coefficient")),
    _ksi(getMaterialProperty<Real>("int_solid_energy_coefficient")),
    _hv(declareProperty<Real>("hv")),
    _dhv(declarePropertyDerivative<Real>("hv", _phi_name)),
    _d2hv(declarePropertyDerivative<Real>("hv", _phi_name, _phi_name)),
    _hs(declareProperty<Real>("hs")),
    _dhs(declarePropertyDerivative<Real>("hs", _phi_name)),
    _d2hs(declarePropertyDerivative<Real>("hs", _phi_name, _phi_name)),
    // Chi
    _chiu(declareProperty<Real>("chiu")),
    _dchiudphi(declarePropertyDerivative<Real>("chiu", _phi_name)),
    _dchiudw(declarePropertyDerivative<Real>("chiu", _wu_name)),
    _d2chiudphi2(declarePropertyDerivative<Real>("chiu", _phi_name, _phi_name)),
    _d2chiudw2(declarePropertyDerivative<Real>("chiu", _wu_name, _wu_name)),
    _d2chiudphidw(declarePropertyDerivative<Real>("chiu", _phi_name, _wu_name)),
    _chii(declareProperty<Real>("chii")),
    _dchiidphi(declarePropertyDerivative<Real>("chii", _phi_name)),
    _dchiidw(declarePropertyDerivative<Real>("chii", _wi_name)),
    _d2chiidphi2(declarePropertyDerivative<Real>("chii", _phi_name, _phi_name)),
    _d2chiidw2(declarePropertyDerivative<Real>("chii", _wi_name, _wi_name)),
    _d2chiidphidw(declarePropertyDerivative<Real>("chii", _phi_name, _wi_name)),
    // Defect densities
    // Void phase
    _rhovu(declareProperty<Real>("rhovu")),
    _drhovudw(declarePropertyDerivative<Real>("rhovu", _wu_name)),
    _rhovi(declareProperty<Real>("rhovi")),
    _drhovidw(declarePropertyDerivative<Real>("rhovi", _wi_name)),
    // Solid phase
    _rhosu(declareProperty<Real>("rhosu")),
    _drhosudw(declarePropertyDerivative<Real>("rhosu", _wu_name)),
    _d2rhosudw2(declarePropertyDerivative<Real>("rhosu", _wu_name, _wu_name)),
    _drhosu(_neta),
    _d2rhosudwdeta(_neta),
    _d2rhosu(_neta),
    _rhosi(declareProperty<Real>("rhosi")),
    _drhosidw(declarePropertyDerivative<Real>("rhosi", _wi_name)),
    _d2rhosidw2(declarePropertyDerivative<Real>("rhosi", _wi_name, _wi_name)),
    _drhosi(_neta),
    _d2rhosidwdeta(_neta),
    _d2rhosi(_neta),
    // Omega stay just v/s but need derivates wrt both chempots (omegas_dwdeta too)
    _omegav(declareProperty<Real>("omegav")),
    _domegavdwu(declarePropertyDerivative<Real>("omegav", _wu_name)),
    _d2omegavdwu2(declarePropertyDerivative<Real>("omegav", _wu_name, _wu_name)),
    _domegavdwi(declarePropertyDerivative<Real>("omegav", _wi_name)),
    _d2omegavdwi2(declarePropertyDerivative<Real>("omegav", _wi_name, _wi_name)),
    _omegas(declareProperty<Real>("omegas")),
    _domegasdwu(declarePropertyDerivative<Real>("omegas", _wu_name)),
    _d2omegasdwu2(declarePropertyDerivative<Real>("omegas", _wu_name, _wu_name)),
    _domegasdwi(declarePropertyDerivative<Real>("omegas", _wi_name)),
    _d2omegasdwi2(declarePropertyDerivative<Real>("omegas", _wi_name, _wi_name)),
    _domegasdeta(_neta),
    _d2omegasdwudeta(_neta),
    _d2omegasdwideta(_neta),
    _d2omegasdetadeta(_neta),
    // Thermodynamic parameters
    _mu(declareProperty<Real>("mu")),
    _dmu(declarePropertyDerivative<Real>("mu", _phi_name)),
    _d2mu(declarePropertyDerivative<Real>("mu", _phi_name, _phi_name)),
    _kappa(declareProperty<Real>("kappa")),
    _dkappa(declarePropertyDerivative<Real>("kappa", _phi_name)),
    _d2kappa(declarePropertyDerivative<Real>("kappa", _phi_name, _phi_name)),
    _gamma(declareProperty<Real>("gamma")),
    // mass conservation
    _hv_c_min(declareProperty<Real>("hv_c_min")),
    _hs_c_min(declareProperty<Real>("hs_c_min")),
    _hv_over_kVa(declareProperty<Real>("hv_over_kVa")),
    _hs_over_kVa(declareProperty<Real>("hs_over_kVa")),

    // _sigma_s(getParam<Real>("surface_energy")),
    // _sigma_gb(getParam<Real>("grainboundary_energy")),
    _sigma_s(getMaterialProperty<Real>("surface_energy")),
    _sigma_gb(getMaterialProperty<Real>("grainboundary_energy")),
    _int_width(getParam<Real>("int_width")),
    _switch(getParam<Real>("surface_switch_value")),
    _Va(getParam<Real>("atomic_volume")),
    _solid_energy(getParam<MooseEnum>("solid_energy_model")),
    // _mu_s(6.0 * _sigma_s / _int_width),
    // _mu_gb(6.0 * _sigma_gb / _int_width),
    // _kappa_s(0.75 * _sigma_s * _int_width),
    // _kappa_gb(0.75 * _sigma_gb * _int_width),
    _kB(8.617343e-5), // eV/K
    _mass_conservation(getParam<bool>("mass_conservation"))
{
  if ((_switch > 1.0) || (_switch < 0.0))
    mooseError(
        "GrandPotentialMultiSinteringMaterial: surface_switch_value should be between 0 and 1");

  if (_mass_conservation && _solid_energy > 0)
    mooseError("GrandPotentialMultiSinteringMaterial: strict mass conservation is currently only "
               "applicable to parabolic free energy");

  if (_mass_conservation)
    mooseError(
        "GrandPotentialMultiSinteringMaterial: strict mass conservation is currently disabled");

  for (unsigned int i = 0; i < _neta; ++i)
  {
    _eta[i] = &coupledValue("etas", i);
    _eta_name[i] = coupledName("etas", i);
    // vac and int cs deta
    _dcsu_eq[i] = &getMaterialPropertyDerivativeByName<Real>(_csu_eq_name, _eta_name[i]);
    _d2csu_eq[i].resize(_neta);
    _dcsi_eq[i] = &getMaterialPropertyDerivativeByName<Real>(_csi_eq_name, _eta_name[i]);
    _d2csi_eq[i].resize(_neta);
    // vac and int rhos deta
    _drhosu[i] = &declarePropertyDerivative<Real>("rhosu", _eta_name[i]);
    _d2rhosu[i].resize(_neta);
    _d2rhosudwdeta[i] = &declarePropertyDerivative<Real>("rhosu", _wu_name, _eta_name[i]);
    _drhosi[i] = &declarePropertyDerivative<Real>("rhosi", _eta_name[i]);
    _d2rhosi[i].resize(_neta);
    _d2rhosidwdeta[i] = &declarePropertyDerivative<Real>("rhosi", _wi_name, _eta_name[i]);
    // omegas deta
    _domegasdeta[i] = &declarePropertyDerivative<Real>("omegas", _eta_name[i]);
    _d2omegasdwudeta[i] = &declarePropertyDerivative<Real>("omegas", _wu_name, _eta_name[i]);
    _d2omegasdwideta[i] = &declarePropertyDerivative<Real>("omegas", _wi_name, _eta_name[i]);
    _d2omegasdetadeta[i].resize(_neta);

    for (unsigned int j = 0; j <= i; ++j)
    {
      _d2csu_eq[j][i] =
          &getMaterialPropertyDerivativeByName<Real>(_csu_eq_name, _eta_name[j], _eta_name[i]);
      _d2csi_eq[j][i] =
          &getMaterialPropertyDerivativeByName<Real>(_csi_eq_name, _eta_name[j], _eta_name[i]);
      _d2rhosu[j][i] = &declarePropertyDerivative<Real>("rhosu", _eta_name[j], _eta_name[i]);
      _d2rhosi[j][i] = &declarePropertyDerivative<Real>("rhosi", _eta_name[j], _eta_name[i]);
      _d2omegasdetadeta[j][i] =
          &declarePropertyDerivative<Real>("omegas", _eta_name[j], _eta_name[i]);
    }
  }
}

void
GrandPotentialMultiSinteringMaterial::computeQpProperties()
{
  // Calculate phase switching functions
  _hv[_qp] = 0.0;
  _dhv[_qp] = 0.0;
  _d2hv[_qp] = 0.0;

  if (_phi[_qp] >= 1.0)
    _hv[_qp] = 1.0;
  else if (_phi[_qp] > 0.0)
  {
    _hv[_qp] = _phi[_qp] * _phi[_qp] * _phi[_qp] * (10.0 + _phi[_qp] * (-15.0 + _phi[_qp] * 6.0));
    _dhv[_qp] = 30.0 * _phi[_qp] * _phi[_qp] * (_phi[_qp] - 1.0) * (_phi[_qp] - 1.0);
    _d2hv[_qp] = 60.0 * _phi[_qp] * (2.0 * _phi[_qp] - 1.0) * (_phi[_qp] - 1.0);
  }

  _hs[_qp] = 1.0 - _hv[_qp];
  _dhs[_qp] = -_dhv[_qp];
  _d2hs[_qp] = -_d2hv[_qp];

  // Calculate interface switching function
  Real phi = _phi[_qp] / _switch;
  Real f = 0.0;
  Real df = 0.0;
  Real d2f = 0.0;
  if (phi >= 1.0)
    f = 1.0;
  else if (phi > 0.0)
  {
    f = phi * phi * phi * (10.0 + phi * (-15.0 + phi * 6.0));
    df = 30.0 / _switch * phi * phi * (phi - 1.0) * (phi - 1.0);
    d2f = 60.0 * phi / (_switch * _switch) * (2.0 * phi - 1.0) * (phi - 1.0);
  }

  // Equilibrium void phase defect concentrations
  Real cvu_eq = 1.0;
  Real cvi_eq = 0.0;

  // Calculate the void phase density and potentials
  _rhovu[_qp] = _wu[_qp] / (_Va * _Va * _kvu[_qp]) + cvu_eq / _Va;
  _drhovudw[_qp] = 1.0 / (_Va * _Va * _kvu[_qp]);
  _rhovi[_qp] = _wi[_qp] / (_Va * _Va * _kvi[_qp]) + cvi_eq / _Va;
  _drhovidw[_qp] = 1.0 / (_Va * _Va * _kvi[_qp]);

  _omegav[_qp] = -0.5 * _wu[_qp] * _wu[_qp] / (_Va * _Va * _kvu[_qp]) - _wu[_qp] * cvu_eq / _Va +
                 -0.5 * _wi[_qp] * _wi[_qp] / (_Va * _Va * _kvi[_qp]) - _wi[_qp] * cvi_eq / _Va;
  _domegavdwu[_qp] = -_rhovu[_qp];
  _d2omegavdwu2[_qp] = -_drhovudw[_qp];
  _domegavdwi[_qp] = -_rhovi[_qp];
  _d2omegavdwi2[_qp] = -_drhovidw[_qp];

  // Calculate solid phase density and potential
  Real d3rhosudw3 = 0;
  Real d3rhosidw3 = 0;
  switch (_solid_energy)
  {
    case 0: // PARABOLIC
    {
      _rhosu[_qp] = _wu[_qp] / (_Va * _Va * _ksu[_qp]) + _csu_eq[_qp] / _Va;
      _drhosudw[_qp] = 1.0 / (_Va * _Va * _ksu[_qp]);
      _d2rhosudw2[_qp] = 0.0;
      d3rhosudw3 = 0.0;
      _rhosi[_qp] = _wi[_qp] / (_Va * _Va * _ksi[_qp]) + _csi_eq[_qp] / _Va;
      _drhosidw[_qp] = 1.0 / (_Va * _Va * _ksi[_qp]);
      _d2rhosidw2[_qp] = 0.0;
      d3rhosidw3 = 0.0;

      _omegas[_qp] =
          -0.5 * _wu[_qp] * _wu[_qp] / (_Va * _Va * _ksu[_qp]) - _wu[_qp] * _csu_eq[_qp] / _Va +
          -0.5 * _wi[_qp] * _wi[_qp] / (_Va * _Va * _ksi[_qp]) - _wi[_qp] * _csi_eq[_qp] / _Va;
      _domegasdwu[_qp] = -_rhosu[_qp];
      _d2omegasdwu2[_qp] = -_drhosudw[_qp];
      _domegasdwi[_qp] = -_rhosi[_qp];
      _d2omegasdwi2[_qp] = -_drhosidw[_qp];

      // // bodyforce and matreact coefficients for strict mass conservation case
      // _hv_c_min[_qp] = _hv[_qp] * 1.0;
      // _hs_c_min[_qp] = _hs[_qp] * _cs_eq[_qp];
      // _hv_over_kVa[_qp] = _hv[_qp] / (_Va * _kv[_qp]);
      // _hs_over_kVa[_qp] = _hs[_qp] / (_Va * _ks[_qp]);

      for (unsigned int i = 0; i < _neta; ++i)
      {
        (*_drhosu[i])[_qp] = (*_dcsu_eq[i])[_qp] / _Va;
        (*_d2rhosudwdeta[i])[_qp] = 0.0;
        (*_drhosi[i])[_qp] = (*_dcsi_eq[i])[_qp] / _Va;
        (*_d2rhosidwdeta[i])[_qp] = 0.0;
        (*_domegasdeta[i])[_qp] =
            -_wu[_qp] * (*_dcsu_eq[i])[_qp] / _Va - _wi[_qp] * (*_dcsi_eq[i])[_qp] / _Va;
        (*_d2omegasdwudeta[i])[_qp] = -(*_dcsu_eq[i])[_qp] / _Va;
        (*_d2omegasdwideta[i])[_qp] = -(*_dcsi_eq[i])[_qp] / _Va;
        for (unsigned int j = i; j < _neta; ++j)
        {
          (*_d2rhosu[i][j])[_qp] = (*_d2csu_eq[i][j])[_qp] / _Va;
          (*_d2rhosi[i][j])[_qp] = (*_d2csi_eq[i][j])[_qp] / _Va;
          (*_d2omegasdetadeta[i][j])[_qp] =
              -_wu[_qp] * (*_d2csu_eq[i][j])[_qp] / _Va - _wi[_qp] * (*_d2csi_eq[i][j])[_qp] / _Va;
        }
      }
      break;
    }       // case 0; // PARABOLIC
    case 1: // DILUTE
    {
      mooseError("Cannot run Dilute");
      break;
    }       // case 1: // DILUTE
    case 2: // IDEAL
    {
      mooseError("Cannot run Ideal");
      break;
    } // case 2: // IDEAL
  }   // switch (_solid_energy)

  // Calculate the susceptibility
  // Vacancy susceptibility
  _chiu[_qp] = _hs[_qp] * _drhosudw[_qp] + _hv[_qp] * _drhovudw[_qp];
  _dchiudphi[_qp] = _dhs[_qp] * _drhosudw[_qp] + _dhv[_qp] * _drhovudw[_qp];
  _dchiudw[_qp] = _hs[_qp] * _d2rhosudw2[_qp];
  _d2chiudphi2[_qp] = _d2hs[_qp] * _drhosudw[_qp] + _d2hv[_qp] * _drhovudw[_qp];
  _d2chiudw2[_qp] = _hs[_qp] * d3rhosudw3;
  _d2chiudphidw[_qp] = _dhs[_qp] * _d2rhosudw2[_qp];
  // Interstitial susceptibility
  _chii[_qp] = _hs[_qp] * _drhosidw[_qp] + _hv[_qp] * _drhovidw[_qp];
  _dchiidphi[_qp] = _dhs[_qp] * _drhosidw[_qp] + _dhv[_qp] * _drhovidw[_qp];
  _dchiidw[_qp] = _hs[_qp] * _d2rhosidw2[_qp];
  _d2chiidphi2[_qp] = _d2hs[_qp] * _drhosidw[_qp] + _d2hv[_qp] * _drhovidw[_qp];
  _d2chiidw2[_qp] = _hs[_qp] * d3rhosidw3;
  _d2chiidphidw[_qp] = _dhs[_qp] * _d2rhosidw2[_qp];

  // Calculate the gb and s components here now that sigma is a material
  Real mu_s = 6.0 * _sigma_s[_qp] / _int_width;
  Real mu_gb = 6.0 * _sigma_gb[_qp] / _int_width;
  Real kappa_s = 0.75 * _sigma_s[_qp] * _int_width;
  Real kappa_gb = 0.75 * _sigma_gb[_qp] * _int_width;

  // thermodynamic parameters
  _mu[_qp] = mu_gb + (mu_s - mu_gb) * f;
  _kappa[_qp] = kappa_gb + (kappa_s - kappa_gb) * f;
  _dmu[_qp] = (mu_s - mu_gb) * df;
  _dkappa[_qp] = (kappa_s - kappa_gb) * df;
  _d2mu[_qp] = (mu_s - mu_gb) * d2f;
  _d2kappa[_qp] = (kappa_s - kappa_gb) * d2f;
  _gamma[_qp] = 1.5;
} // void GrandPotentialSinteringMaterial::computeQpProperties()
