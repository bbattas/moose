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

#include "GrandPotentialRadiationDamage.h"
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

registerMooseObject("MarmotApp", GrandPotentialRadiationDamage);

template <>
InputParameters validParams<GrandPotentialRadiationDamage>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Includes switching functions and thermodynamic properties for a grand potential 2-phase model with interstitials");
  // Coupled variables
  params.addRequiredCoupledVarWithAutoBuild("etas", "var_name_base", "op_num", "Array of order parameters that define the solid phase");
  params.addRequiredCoupledVar("chemical_potentials", "Vector of chemical potential variable names");
  params.addRequiredCoupledVar("void_op", "Order paramter that defines the void phase");
  params.addRequiredCoupledVar("Temperature", "Name of the absolute temperature variable");
  // Naming conventions
  params.addRequiredParam<std::vector<std::string>>("potential_subscripts", "vector of strings to differentiate material properties for the different chemical potential variables");
  params.addParam<std::string>("solid_subscript", "s", "Subscript to distinguish material properties for the solid phase");
  params.addParam<std::string>("void_subscript", "v", "Subscript to distinguish material properties for the void phase");
  params.addParam<std::string>("switching_base", "h", "Base name for switching functions");
  params.addParam<std::string>("susceptibility_base", "chi", "Base name for susceptibilities");
  params.addParam<std::string>("density_base", "rho", "Base name for densities");
  params.addParam<std::string>("potential_base", "omega", "Base name for potential density functions");
  params.addParam<std::string>("kappa_name", "kappa", "Name of the kappa function");
  params.addParam<std::string>("gamma_name", "gamma", "Name of the gamma function");
  // Material Properties
  params.addParam<std::vector<MaterialPropertyName>>("solid_energy_coefficient", "Vector of names of the parabolic energy coefficients for the solid phase free energies. Only needed when 'solid_energy_model = PARABOLIC'.");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("void_energy_coefficient", "Vector of names of the parabolic energy coefficients for the void phase free energies");
  params.addRequiredParam<std::vector<MaterialPropertyName>>("solid_equilibrium_concentrations", "Vector of names of materials that defines the equilibrium concentrations for the solid phase. The entries in the vector correspond to the chemical potential vector.");
  MooseEnum solid_energy_model("PARABOLIC DILUTE IDEAL", "PARABOLIC");
  params.addParam<MooseEnum>("solid_energy_model", solid_energy_model, "Type of energy function to use for the solid phase.");
  // Constants
  params.addRequiredParam<std::vector<Real>>("void_equilibrium_concentrations", "Vector of constant values for the void-phase equilibrium concentration values associated with the chemical potential variables.");
  params.addRequiredParam<Real>("surface_energy", "Surface energy (energy/area)");
  params.addRequiredParam<Real>("GB_energy", "Grain boundary energy (energy/area)");
  params.addParam<Real>("int_width", 1, "interface width (length)");
  params.addParam<Real>("interface_switch_value", 0.3, "Value between 0 and 1 that controls at what value of the void_op that mu and kappa begin transition between surface and GB values");
  params.addRequiredParam<Real>("atomic_volume", "Atomic volume of the species (volume)");

  return params;
}

GrandPotentialRadiationDamage::GrandPotentialRadiationDamage(const InputParameters & parameters) : DerivativeMaterialInterface<Material>(parameters),
  // get properties from inputs
  _nop(coupledComponents("etas")),
  _ops(_nop),
  _op_names(_nop),
  _phi(coupledValue("void_op")),
  _phi_name(getVar("void_op", 0)->name()),
  _nw(coupledComponents("chemical_potentials")),
  _w(_nw),
  _w_names(_nw),
  _T(coupledValue("Temperature")),
  _sub_w(getParam<std::vector<std::string>>("potential_subscripts")),
  _sub_s(getParam<std::string>("solid_subscript")),
  _sub_v(getParam<std::string>("void_subscript")),
  _h_base(getParam<std::string>("switching_base")),
  _chi_base(getParam<std::string>("susceptibility_base")),
  _rho_base(getParam<std::string>("density_base")),
  _omega_base(getParam<std::string>("potential_base")),
  _kappa_name(getParam<std::string>("kappa_name")),
  _gamma_name(getParam<std::string>("gamma_name")),
  _ksolid_names(getParam<std::vector<MaterialPropertyName>>("solid_energy_coefficient")),
  _ksolid(_nw),
  _kvoid_names(getParam<std::vector<MaterialPropertyName>>("void_energy_coefficient")),
  _kvoid(_nw),
  _cs_eq_names(getParam<std::vector<MaterialPropertyName>>("solid_equilibrium_concentrations")),
  _cs_eq(_nw),
  _dcs_eqdeta(_nw),
  _d2cs_eqdeta2(_nw),
  _solid_energy(getParam<MooseEnum>("solid_energy_model")),
  _cv_eq(getParam<std::vector<Real>>("void_equilibrium_concentrations")),
  _sigma_s(getParam<Real>("surface_energy")),
  _sigma_gb(getParam<Real>("GB_energy")),
  _int_width(getParam<Real>("int_width")),
  _phi0(getParam<Real>("interface_switch_value")),
  _Va(getParam<Real>("atomic_volume")),
  // declare material properties & derivatives
  // switching functions
  _hs(declareProperty<Real>(_h_base + _sub_s)),
  _dhsdphi(declarePropertyDerivative<Real>(_h_base + _sub_s, _phi_name)),
  _d2hsdphi2(declarePropertyDerivative<Real>(_h_base + _sub_s, _phi_name, _phi_name)),
  _hv(declareProperty<Real>(_h_base + _sub_v)),
  _dhvdphi(declarePropertyDerivative<Real>(_h_base + _sub_v, _phi_name)),
  _d2hvdphi2(declarePropertyDerivative<Real>(_h_base + _sub_v, _phi_name, _phi_name)),
  // susceptibilities
  _chi(_nw),
  _dchidw(_nw),
  _dchidphi(_nw),
  _d2chidw2(_nw),
  _d2chidphidw(_nw),
  _d2chidphi2(_nw),
  _dchideta(_nw),
  _d2chidetadphi(_nw),
  _d2chidetadw(_nw),
  _d2chideta2(_nw),
  // densities
  _rhos(_nw),
  _drhosdw(_nw),
  _d2rhosdw2(_nw),
  _drhosdeta(_nw),
  _d2rhosdetadw(_nw),
  _d2rhosdeta2(_nw),
  _rhov(_nw),
  _drhovdw(_nw),
  // potential densities
  _omegas(declareProperty<Real>(_omega_base + _sub_s)),
  _domegasdw(_nw),
  _d2omegasdw2(_nw),
  _domegasdeta(_nop),
  _d2omegasdetadw(_nw),
  _d2omegasdeta2(_nop),
  _omegav(declareProperty<Real>(_omega_base + _sub_v)),
  _domegavdw(_nw),
  _d2omegavdw2(_nw),
  // energy coefficients DONE
  _mu(declareProperty<Real>("mu")),
  _dmudphi(declarePropertyDerivative<Real>("mu", _phi_name)),
  _d2mudphi2(declarePropertyDerivative<Real>("mu", _phi_name, _phi_name)),
  _kappa(declareProperty<Real>(_kappa_name)),
  _dkappadphi(declarePropertyDerivative<Real>(_kappa_name, _phi_name)),
  _d2kappadphi2(declarePropertyDerivative<Real>(_kappa_name, _phi_name, _phi_name)),
  _gamma(declareProperty<Real>(_gamma_name))
{
  // Check vector lengths
  if (_sub_w.size() != _nw)
    mooseError("GrandPotentialRadiationDamage: 'potential_subscripts' should be a vector of the same length as 'chemical_potentials'.");
  if (_cs_eq_names.size() != _nw)
    mooseError("GrandPotentialRadiationDamage: 'solid_equilibrium_concentrations' should be a vector of the same length as 'chemical_potentials'.");
  if (_cv_eq.size() != _nw)
    mooseError("GrandPotentialRadiationDamage: 'void_equilibrium_concentrations' should be a vector of the same length as 'chemical_potentials'.");
  if (_solid_energy == "PARABOLIC")
    if (_ksolid_names.size() != _nw)
      mooseError("GrandPotentialRadiationDamage: 'solid_energy_coefficient' should be a vector of the same length as 'chemical_potentials'.");
  if (_kvoid_names.size() != _nw)
    mooseError("GrandPotentialRadiationDamage: 'void_energy_coefficient' should be a vector of the same length as 'chemical_potentials'.");
  // additional checks
  if (_phi0 < 0.0 || _phi0 > 1.0)
    mooseError("GrandPotentialRadiationDamage: The value 'interface_switch_value' should be between 0 and 1.");
  // declare order parameter vectors
  for (unsigned int iop = 0; iop < _nop; ++iop)
  {
    _ops[iop] = &coupledValue("etas", iop);
    _op_names[iop] = getVar("etas", iop)->name();
    _domegasdeta[iop] = &declarePropertyDerivative<Real>(_omega_base + _sub_s, _op_names[iop]);
    _d2omegasdeta2[iop].resize(iop+1);
    for (unsigned int jop = 0; jop <= iop; ++jop)
      _d2omegasdeta2[iop][jop] = &declarePropertyDerivative<Real>(_omega_base + _sub_s, _op_names[iop], _op_names[jop]);
  }

  // declare chemical potential vectors
  for (unsigned int i = 0; i < _nw; ++i)
  {
    // chemical potentials
    _w[i] = &coupledValue("chemical_potentials", i);
    _w_names[i] = getVar("chemical_potentials", i)->name();
    if (_solid_energy == "PARABOLIC")
      _ksolid[i] = &getMaterialProperty<Real>(_ksolid_names[i]);
    _kvoid[i] = &getMaterialProperty<Real>(_kvoid_names[i]);
    _cs_eq[i] = &getMaterialProperty<Real>(_cs_eq_names[i]);
    _dcs_eqdeta[i].resize(_nop);
    _d2cs_eqdeta2[i].resize(_nop);

    // susceptibilities
    std::string chiname = _chi_base + _sub_w[i];
    _chi[i] = &declareProperty<Real>(chiname);
    _dchidw[i] = &declarePropertyDerivative<Real>(chiname, _w_names[i]);
    _dchidphi[i] = &declarePropertyDerivative<Real>(chiname, _phi_name);
    _d2chidw2[i] = &declarePropertyDerivative<Real>(chiname, _w_names[i], _w_names[i]);
    _d2chidphidw[i] = &declarePropertyDerivative<Real>(chiname, _phi_name, _w_names[i]);
    _d2chidphi2[i] = &declarePropertyDerivative<Real>(chiname, _phi_name, _phi_name);
    _dchideta[i].resize(_nop);
    _d2chidetadphi[i].resize(_nop);
    _d2chidetadw[i].resize(_nop);
    _d2chideta2[i].resize(_nop);

    // densities
    std::string rhos_name = _rho_base + _sub_s + _sub_w[i];
    std::string rhov_name = _rho_base + _sub_v + _sub_w[i];
    _rhos[i] = &declareProperty<Real>(rhos_name);
    _drhosdw[i] = &declarePropertyDerivative<Real>(rhos_name, _w_names[i]);
    _d2rhosdw2[i] = &declarePropertyDerivative<Real>(rhos_name, _w_names[i], _w_names[i]);
    _drhosdeta[i].resize(_nop);
    _d2rhosdetadw[i].resize(_nop);
    _d2rhosdeta2[i].resize(_nop);
    _rhov[i] = &declareProperty<Real>(rhov_name);
    _drhovdw[i] = &declarePropertyDerivative<Real>(rhov_name, _w_names[i]);

    // potential densities
    _domegasdw[i] = &declarePropertyDerivative<Real>(_omega_base + _sub_s, _w_names[i]);
    _d2omegasdw2[i] = &declarePropertyDerivative<Real>(_omega_base + _sub_v, _w_names[i], _w_names[i]);
    _d2omegasdetadw[i].resize(_nop);
    _domegavdw[i] = &declarePropertyDerivative<Real>(_omega_base + _sub_v, _w_names[i]);
    _d2omegavdw2[i] = &declarePropertyDerivative<Real>(_omega_base + _sub_v, _w_names[i], _w_names[i]);

    //subloops
    for (unsigned int iop = 0; iop < _nop; ++iop)
    {
      _dcs_eqdeta[i][iop] = &getMaterialPropertyDerivativeByName<Real>(_cs_eq_names[i], _op_names[iop]);
      _d2cs_eqdeta2[i][iop].resize(iop+1);
      _dchideta[i][iop] = &declarePropertyDerivative<Real>(chiname, _op_names[iop]);
      _d2chidetadphi[i][iop] = &declarePropertyDerivative<Real>(chiname, _phi_name, _op_names[iop]);
      _d2chidetadw[i][iop] = &declarePropertyDerivative<Real>(chiname, _w_names[i], _op_names[iop]);
      _d2chideta2[i][iop].resize(iop+1);
      _drhosdeta[i][iop] = &declarePropertyDerivative<Real>(rhos_name, _op_names[iop]);
      _d2rhosdetadw[i][iop] = &declarePropertyDerivative<Real>(rhos_name, _w_names[i], _op_names[iop]);
      _d2rhosdeta2[i][iop].resize(iop+1);
      _d2omegasdetadw[i][iop] = &declarePropertyDerivative<Real>(_omega_base + _sub_s, _w_names[i], _op_names[iop]);
      for (unsigned int jop = 0; jop <= iop; ++jop)
      {
        _d2cs_eqdeta2[i][iop][jop] = &getMaterialPropertyDerivativeByName<Real>(_cs_eq_names[i], _op_names[iop], _op_names[jop]);
        _d2chideta2[i][iop][jop] = &declarePropertyDerivative<Real>(chiname, _op_names[iop], _op_names[jop]);
        _d2rhosdeta2[i][iop][jop] = &declarePropertyDerivative<Real>(rhos_name, _op_names[iop], _op_names[jop]);
      }
    }
  } // for (unsigned int i = 0; i < _nw; ++i)
}

void
GrandPotentialRadiationDamage::computeQpProperties()
{
  // Boltzmann's constant
  Real kB = 8.617343e-5;

  // Preallocate omega terms
  _omegas[_qp] = 0.0;
  _omegav[_qp] = 0.0;
  for (unsigned int iop = 0; iop < _nop; ++iop)
  {
    (*_domegasdeta[iop])[_qp] = 0.0;
    for (unsigned int jop = 0; jop <= iop; ++jop)
      (*_d2omegasdeta2[iop][jop])[_qp] = 0.0;
  }

  // Calculate sum_i of eta[i]^2
  Real sumeta = 0.0;
  for (unsigned int iop = 0; iop < _nop; ++iop)
    sumeta += (*_ops[iop])[_qp] * (*_ops[iop])[_qp];

  // Calculate phase switching functions
  std::vector<Real> hvec = switchingBase(_phi[_qp], 1.0);
  _hv[_qp] = hvec[0];
  _dhvdphi[_qp] = hvec[1];
  _d2hvdphi2[_qp] = hvec[2];
  _hs[_qp] = 1.0 - _hv[_qp];
  _dhsdphi[_qp] = - _dhvdphi[_qp];
  _d2hsdphi2[_qp] = -_d2hvdphi2[_qp];

  // Calculate interface switching functions
  hvec = switchingBase(_phi[_qp], _phi0);
  Real hS = hvec[0];
  Real dhS = hvec[1];
  Real d2hS = hvec[2];

  // Calculate the energy terms
  _mu[_qp] = 6.0 / _int_width * (hS * _sigma_s + (1.0 - hS) * _sigma_gb);
  _dmudphi[_qp] = 6.0 / _int_width * dhS * (_sigma_s - _sigma_gb);
  _d2mudphi2[_qp] = 6.0 / _int_width * d2hS * (_sigma_s - _sigma_gb);
  _kappa[_qp] = 0.75 * _int_width * (hS * _sigma_s + (1.0 - hS) * _sigma_gb);
  _dkappadphi[_qp] = 0.75 * _int_width * dhS * (_sigma_s - _sigma_gb);
  _d2kappadphi2[_qp] = 0.75 * _int_width * d2hS * (_sigma_s - _sigma_gb);
  _gamma[_qp] = 1.5;

  // chemical potential-based terms
  for (unsigned int i = 0; i < _nw; ++i)
  {
    // void phase density
    (*_rhov[i])[_qp] = (*_w[i])[_qp] / (_Va * _Va * (*_kvoid[i])[_qp]) + _cv_eq[i] / _Va;
    (*_drhovdw[i])[_qp] = 1.0 / (_Va * _Va * (*_kvoid[i])[_qp]);
    // void potential density
    _omegav[_qp] += -0.5 * (*_w[i])[_qp] * (*_w[i])[_qp] / (_Va * _Va * (*_kvoid[i])[_qp]) - (*_w[i])[_qp] * _cv_eq[i] / _Va;
    (*_domegavdw[i])[_qp] = -(*_rhov[i])[_qp];
    (*_d2omegavdw2[i])[_qp] = -(*_drhovdw[i])[_qp];
    // solid potential density


    Real d3rhosdw3;
    std::vector<Real> d3rhosdetadw2;
    d3rhosdetadw2.resize(_nop);
    std::vector<std::vector<Real>> d3rhosdeta2dw;
    d3rhosdeta2dw.resize(_nop);
    for (unsigned int iop = 0; iop < _nop; ++iop)
      d3rhosdeta2dw[iop].resize(iop+1);

    switch (_solid_energy)
    {
      case 0: //PARABOLIC
      {
        (*_rhos[i])[_qp] = (*_w[i])[_qp] / (_Va * _Va * (*_ksolid[i])[_qp]) + (*_cs_eq[i])[_qp] / _Va;
        (*_drhosdw[i])[_qp] = 1.0 / (_Va * _Va * (*_ksolid[i])[_qp]);
        (*_d2rhosdw2[i])[_qp] = 0.0;
        d3rhosdw3 = 0.0;

        _omegas[_qp] += -0.5 * (*_w[i])[_qp] * (*_w[i])[_qp] / (_Va * _Va * (*_ksolid[i])[_qp]) - (*_w[i])[_qp] * (*_cs_eq[i])[_qp] / _Va;
        (*_domegasdw[i])[_qp] = -(*_rhos[i])[_qp];
        (*_d2omegasdw2[i])[_qp] = -(*_drhosdw[i])[_qp];

        for (unsigned int iop = 0; iop < _nop; ++iop)
        {
          (*_drhosdeta[i][iop])[_qp] = (*_dcs_eqdeta[i][iop])[_qp] / _Va;
          (*_d2rhosdetadw[i][iop])[_qp] = 0.0;
          d3rhosdetadw2[iop] = 0.0;
          (*_domegasdeta[iop])[_qp] += -(*_w[i])[_qp] * (*_dcs_eqdeta[i][iop])[_qp] / _Va;
          (*_d2omegasdetadw[i][iop])[_qp] = -(*_drhosdeta[i][iop])[_qp];

          for (unsigned int jop = 0; jop <= iop; ++jop)
          {
            d3rhosdeta2dw[iop][jop] = 0.0;
            (*_d2rhosdeta2[i][iop][jop])[_qp] = (*_d2cs_eqdeta2[i][iop][jop])[_qp] / _Va;
            (*_d2omegasdeta2[iop][jop])[_qp] += -(*_w[i])[_qp] * (*_d2cs_eqdeta2[i][iop][jop])[_qp];
          }
        }
        break;
      } // case 0:
      case 1: // DILUTE
      {
        Real rho_exp = std::exp((*_w[i])[_qp] / kB / _T[_qp]);
        (*_rhos[i])[_qp] = (*_cs_eq[i])[_qp] / _Va * rho_exp;
        (*_drhosdw[i])[_qp] = (*_rhos[i])[_qp] / kB / _T[_qp];
        (*_d2rhosdw2[i])[_qp] = (*_drhosdw[i])[_qp] / kB / _T[_qp];
        d3rhosdw3 = (*_d2rhosdw2[i])[_qp] / kB / _T[_qp];

        _omegas[_qp] += kB * _T[_qp] * ((*_cs_eq[i])[_qp] / _Va - (*_rhos[i])[_qp]);
        (*_domegasdw[i])[_qp] = -(*_rhos[i])[_qp];
        (*_d2omegasdw2[i])[_qp] = -(*_drhosdw[i])[_qp];

        for (unsigned int iop = 0; iop < _nop; ++iop)
        {
          (*_drhosdeta[i][iop])[_qp] = (*_dcs_eqdeta[i][iop])[_qp] / _Va * rho_exp;
          (*_d2rhosdetadw[i][iop])[_qp] = (*_drhosdeta[i][iop])[_qp] / kB / _T[_qp];
          d3rhosdetadw2[iop] = (*_d2rhosdetadw[i][iop])[_qp] / kB /_T[_qp];
          (*_domegasdeta[iop])[_qp] += kB * _T[_qp] * (*_dcs_eqdeta[i][iop])[_qp] / _Va * (1.0 - rho_exp);
          (*_d2omegasdetadw[i][iop])[_qp] = -1.0 / _Va * (*_dcs_eqdeta[i][iop])[_qp] * rho_exp;
          for (unsigned int jop = 0; jop <= iop; ++jop)
          {
            (*_d2rhosdeta2[i][iop][jop])[_qp] = (*_d2cs_eqdeta2[i][iop][jop])[_qp] * rho_exp / _Va;
            d3rhosdeta2dw[iop][jop] = (*_d2rhosdeta2[i][iop][jop])[_qp] / kB / _T[_qp];
            (*_d2omegasdeta2[iop][jop])[_qp] += kB * _T[_qp] * (*_d2cs_eqdeta2[i][iop][jop])[_qp] / _Va * (1.0 - rho_exp);
          }
        }
        break;
      } // case 1:
      case 2: //IDEAL
      {
        Real Ef = -kB * _T[_qp] * std::log((*_cs_eq[i])[_qp] / (1.0 - (*_cs_eq[i])[_qp]));
        std::vector<Real> dEfdeta;
        std::vector<std::vector<Real>> d2Efdeta2;
        dEfdeta.resize(_nop);
        d2Efdeta2.resize(_nop);

        Real x = std::exp(((*_w[i])[_qp] - Ef) / (kB * _T[_qp]));
        Real x0 = std::exp(-Ef / (kB * _T[_qp]));
        (*_rhos[i])[_qp] = x / ((1.0 + x) * _Va);
        Real rhos0 = x0 / ((1.0 + x0) * _Va);
        (*_drhosdw[i])[_qp] = x / (Utility::pow<2>(1.0 + x) * _Va * kB * _T[_qp]);
        (*_d2rhosdw2[i])[_qp] = x * (1.0 - x) / (_Va * Utility::pow<2>(kB * _T[_qp]) * Utility::pow<3>(1.0 + x));
        d3rhosdw3 = x * (1.0 - 4.0 * x + x * x) / (_Va * Utility::pow<3>(kB * _T[_qp]) * Utility::pow<4>(1.0 + x));
        Real logx0 = x0;
        Real logx = x;
        if (x0 > 1e-16) //In cases of very small values it is better to use a linear approximation
          logx0 = std::log(1.0 + x0);
        if (x > 1e-16)
          logx = std::log(1.0 + x);

        _omegas[_qp] += kB * _T[_qp] / _Va * (logx0 - logx);
        (*_domegasdw[i])[_qp] = -(*_rhos[i])[_qp];
        (*_d2omegasdw2[i])[_qp] = -(*_drhosdw[i])[_qp];

        for (unsigned int iop = 0; iop < _nop; ++iop)
        {
          dEfdeta[iop] = -kB * _T[_qp] * (*_dcs_eqdeta[i][iop])[_qp] * (1.0 / (*_cs_eq[i])[_qp] + 1.0 / (1.0 - (*_cs_eq[i])[_qp]));
          d2Efdeta2[iop].resize(iop+1);

          (*_drhosdeta[i][iop])[_qp] = -dEfdeta[iop] * (*_drhosdw[i])[_qp];
          (*_d2rhosdetadw[i][iop])[_qp] = -dEfdeta[iop] * (*_d2rhosdw2[i])[_qp];
          d3rhosdetadw2[iop] = -dEfdeta[iop] * d3rhosdw3;
          (*_domegasdeta[iop])[_qp] += dEfdeta[iop] * ((*_rhos[i])[_qp] - rhos0);
          (*_d2omegasdetadw[i][iop])[_qp] = dEfdeta[iop] * (*_drhosdw[i])[_qp];

          for (unsigned int jop = 0; jop <= iop; ++jop)
          {
            d2Efdeta2[iop][jop] = -kB * _T[_qp] *
                         ((*_d2cs_eqdeta2[i][iop][jop])[_qp] * (1.0 / (*_cs_eq[i])[_qp] + 1.0 / (1.0 - (*_cs_eq[i])[_qp])) +
                          (*_dcs_eqdeta[i][iop])[_qp] * (*_dcs_eqdeta[i][jop])[_qp] *
                          (1.0 / ((1.0 - (*_cs_eq[i])[_qp]) * (1.0 - (*_cs_eq[i])[_qp])) -
                           1.0 / ((*_cs_eq[i])[_qp] * (*_cs_eq[i])[_qp])));
            (*_d2rhosdeta2[i][iop][jop])[_qp] = -d2Efdeta2[iop][jop] * (*_drhosdw[i])[_qp] +
                         dEfdeta[iop] * dEfdeta[jop] * (*_d2rhosdw2[i])[_qp];
            d3rhosdeta2dw[iop][jop] = - d2Efdeta2[iop][jop] * (*_d2rhosdw2[i])[_qp] +
                         dEfdeta[iop] * dEfdeta[jop] * d3rhosdw3;
            (*_d2omegasdeta2[iop][jop])[_qp] += d2Efdeta2[iop][jop] *
                         ((*_rhos[i])[_qp] - rhos0) + dEfdeta[iop] * dEfdeta[jop] /
                          (_Va * kB * _T[_qp]) * (x / Utility::pow<2>(1.0 + x) -
                           x0 / Utility::pow<2>(1.0 + x0));
          }
        }
        break;
      } // case 2:
    } // switch (_solid_energy)

    // susceptibility
    (*_chi[i])[_qp] = _hs[_qp] * (*_drhosdw[i])[_qp] + _hv[_qp] * (*_drhovdw[i])[_qp];
    (*_dchidphi[i])[_qp] = _dhsdphi[_qp] * (*_drhosdw[i])[_qp] + _dhvdphi[_qp] * (*_drhovdw[i])[_qp];
    (*_dchidw[i])[_qp] = _hs[_qp] * (*_d2rhosdw2[i])[_qp];
    (*_d2chidphi2[i])[_qp] = _d2hsdphi2[_qp] * (*_drhosdw[i])[_qp] + _d2hvdphi2[_qp] * (*_drhovdw[i])[_qp];
    (*_d2chidw2[i])[_qp] = _hs[_qp] * d3rhosdw3;
    (*_d2chidphidw[i])[_qp] = _dhsdphi[_qp] * (*_d2rhosdw2[i])[_qp];

    for (unsigned int iop = 0; iop < _nop; ++iop)
    {
      (*_dchideta[i][iop])[_qp] = _hs[_qp] * (*_d2rhosdetadw[i][iop])[_qp];
      (*_d2chidetadphi[i][iop])[_qp] = _dhsdphi[_qp] * (*_d2rhosdetadw[i][iop])[_qp];
      (*_d2chidetadw[i][iop])[_qp] = _hs[_qp] * d3rhosdetadw2[iop];
      for (unsigned int jop = 0; jop <= iop; ++jop)
        (*_d2chideta2[i][iop][jop])[_qp] = _hs[_qp] * d3rhosdeta2dw[iop][jop];
    }

  } // for (unsigned int i = 0; i < _nw; ++i)
} // GrandPotentialRadiationDamage::computeQpProperties()

std::vector<Real>
GrandPotentialRadiationDamage::switchingBase(Real x, Real x0)
{
  // Base function returns switching function value and first and second derivatives
  std::vector<Real> vec;
  vec.resize(3);
  vec[0] = 0.0;
  vec[1] = 0.0;
  vec[2] = 0.0;

  if (x >= x0)
  {
    vec[0] = 1.0;
  }
  else if (x > 0.0)
  {
    Real t = x / x0;
    vec[0] = t * t * t * (6.0 * t * t - 15.0 * t + 10.0);
    vec[1] = t * t * (30.0 * t * t - 60.0 * t + 30.0) / x0;
    vec[2] = t * (120.0 * t * t - 180.0 * t + 60.0) / (x0 * x0);
  }
  return vec;
}
