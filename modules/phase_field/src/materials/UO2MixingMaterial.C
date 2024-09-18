#include "UO2MixingMaterial.h"

// libMesh includes might not need with computeQpProp?
#include "libmesh/quadrature.h"

registerMooseObject("PhaseFieldApp", UO2MixingMaterial);

InputParameters
UO2MixingMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("UO2 ballistic mixing materials for vacancies "
                             "and interstitials.");
  params.addRequiredCoupledVarWithAutoBuild(
      "etas", "var_name_base", "op_num", "Array of coupled variables"); //
  // Variable inputs
  params.addRequiredCoupledVar("chemical_potential_vac",
                               "The name of the vacancy chemical potential variable");
  params.addRequiredCoupledVar("chemical_potential_int",
                               "The name of the interstitial chemical potential variable");
  params.addRequiredCoupledVar("concentration_vac",
                               "The name of the vacancy concentration variable");
  params.addRequiredCoupledVar("concentration_int",
                               "The name of the interstitial concentration variable");
  params.addRequiredCoupledVar("void_op", "The name of the void phase order parameter");
  // params.addCoupledVar("T", "Temperature variable in Kelvin");
  // Material Inputs
  params.addRequiredParam<MaterialPropertyName>("vac_chi",
                                                "Vacancy susceptability (* Va from GPMSM)");
  params.addRequiredParam<MaterialPropertyName>("int_chi",
                                                "Interstitial susceptability (* Va from GPMSM)");
  // Mixing Output Names
  params.addParam<std::string>(
      "vac_mix_out", "cv_mixing", "Name for the ballistic mixing output material for vacancies");
  params.addParam<std::string>("int_mix_out",
                               "ci_mixing",
                               "Name for the ballistic mixing output material for interstitials");
  // General Constants
  params.addParam<Real>("f_dot", 1e-8, "Fission rate in units of problem (fissions/length^3*s)");
  params.addParam<Real>("atomic_volume", 0.04092, "Atomic volume of material");
  // Cascade constants
  params.addParam<Real>("Nc", 2, "Number of cascades per fission");
  params.addParam<Real>("Vc", 268, "Cascade volume (length^3)");
  params.addParam<Real>("tc", 1e-11, "Cascade time (s)");
  params.addParam<Real>("Dc", 1e12, "Diffusivity within a cascade (length^2/s)");

  return params;
}

UO2MixingMaterial::UO2MixingMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    // Grain OPs
    _neta(coupledComponents("etas")),
    _eta(_neta),
    _eta_name(_neta),
    // chemical potentials
    _wu(coupledValue("chemical_potential_vac")),
    _wu_name(coupledName("chemical_potential_vac", 0)),
    _wi(coupledValue("chemical_potential_int")),
    _wi_name(coupledName("chemical_potential_int", 0)),
    // Concentrations
    _cu(coupledValue("concentration_vac")),
    _cu_name(coupledName("concentration_vac", 0)),
    _ci(coupledValue("concentration_int")),
    _ci_name(coupledName("concentration_int", 0)),
    // Void Phase OP
    _phi(coupledValue("void_op")),
    _phi_name(coupledName("void_op", 0)),
    // Chi*Va from GPMultiSinteringMaterial
    _chiu(getMaterialProperty<Real>("vac_chi")),
    _dchiudphi(getMaterialPropertyDerivativeByName<Real>("vac_chi", _phi_name)),
    _dchiudwu(getMaterialPropertyDerivativeByName<Real>("vac_chi", _wu_name)),
    _chii(getMaterialProperty<Real>("int_chi")),
    _dchiidphi(getMaterialPropertyDerivativeByName<Real>("int_chi", _phi_name)),
    _dchiidwi(getMaterialPropertyDerivativeByName<Real>("int_chi", _wi_name)),
    // Output Mixing Names
    _bmv_name(getParam<std::string>("vac_mix_out")),
    _bmv(declareProperty<Real>(_bmv_name)),
    _dbmvdphi(declarePropertyDerivative<Real>(_bmv_name, _phi_name)),
    _dbmvdwu(declarePropertyDerivative<Real>(_bmv_name, _wu_name)),
    _dbmvdcu(declarePropertyDerivative<Real>(_bmv_name, _cu_name)),
    _bmi_name(getParam<std::string>("int_mix_out")),
    _bmi(declareProperty<Real>(_bmi_name)),
    _dbmidphi(declarePropertyDerivative<Real>(_bmi_name, _phi_name)),
    _dbmidwi(declarePropertyDerivative<Real>(_bmi_name, _wi_name)),
    _dbmidci(declarePropertyDerivative<Real>(_bmi_name, _ci_name)),
    // CONSTANTS
    _fdot(getParam<Real>("f_dot")),
    _Va(getParam<Real>("atomic_volume")),
    // Cascade Constants
    _Nc(getParam<Real>("Nc")),
    _Vc(getParam<Real>("Vc")),
    _tc(getParam<Real>("tc")),
    _Dc(getParam<Real>("Dc"))
{
  if (_neta == 0)
    mooseError("Model requires op_num > 0");

  for (unsigned int i = 0; i < _neta; ++i) //_op_num
  {
    _eta[i] = &coupledValue("etas", i);    // Grain OP Values
    _eta_name[i] = coupledName("etas", i); // Grain OP Names
  }
}

void
UO2MixingMaterial::computeQpProperties()
{
  // f_dot * noise * Nc * tc * Vc * Dc * chi * Va
  // Mass cons chiu = chi*Va

  // Ballistic Mixing Values
  _bmv[_qp] = _fdot * _Nc * _tc * _Vc * _Dc * _chiu[_qp];
  _bmi[_qp] = _fdot * _Nc * _tc * _Vc * _Dc * _chii[_qp];

  // Derivatives- actually all 0 as long as ks = kv because then chiu is constant
  // Vacancy
  _dbmvdphi[_qp] = _fdot * _Nc * _tc * _Vc * _Dc * _dchiudphi[_qp];
  _dbmvdwu[_qp] = _fdot * _Nc * _tc * _Vc * _Dc * _dchiudwu[_qp];
  _dbmvdcu[_qp] = _fdot * _Nc * _tc * _Vc * _Dc * _dchiudwu[_qp] / (_chiu[_qp] / _Va);
  // Interstitial
  _dbmidphi[_qp] = _fdot * _Nc * _tc * _Vc * _Dc * _dchiidphi[_qp];
  _dbmidwi[_qp] = _fdot * _Nc * _tc * _Vc * _Dc * _dchiidwi[_qp];
  _dbmidci[_qp] = _fdot * _Nc * _tc * _Vc * _Dc * _dchiidwi[_qp] / (_chii[_qp] / _Va);
}
