#include "UO2RecombinationMaterial.h"

// libMesh includes might not need with computeQpProp?
#include "libmesh/quadrature.h"

registerMooseObject("PhaseFieldApp", UO2RecombinationMaterial);

InputParameters
UO2RecombinationMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("UO2 material for vacancies and interstitial recombination (1600K).");
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
      "rec_out", "cvi_recomb", "Name for the recombination output material");
  // General Constants
  params.addParam<Real>(
      "hv_thresh", 1e-6, "Upper threshold for hv to allow recombination (in addition to *hs)");
  params.addParam<Real>("atomic_volume", 0.04092, "Atomic volume of material");
  // Cascade constants
  params.addParam<Real>("Z", 250, "Recombination number");
  params.addParam<Real>("a_0", 0.25, "Atomic jump distance (length)");
  params.addParam<Real>("Di_0", 4.0767e11, "Bulk interstitial diffusivity prefactor (length^2/s)");
  params.addParam<Real>("Ei_B", 4.08453089, "Bulk interstitial diffusivity migration energy (eV)");

  return params;
}

UO2RecombinationMaterial::UO2RecombinationMaterial(const InputParameters & parameters)
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
    _chii(getMaterialProperty<Real>("int_chi")),
    // Output Recombination Names
    _rec_name(getParam<std::string>("rec_out")),
    _rec(declareProperty<Real>(_rec_name)),
    _drecdphi(declarePropertyDerivative<Real>(_rec_name, _phi_name)),
    _drecdwu(declarePropertyDerivative<Real>(_rec_name, _wu_name)),
    _drecdwi(declarePropertyDerivative<Real>(_rec_name, _wi_name)),
    _drecdcu(declarePropertyDerivative<Real>(_rec_name, _cu_name)),
    _drecdci(declarePropertyDerivative<Real>(_rec_name, _ci_name)),
    // CONSTANTS
    _hv_thresh(getParam<Real>("hv_thresh")),
    _Va(getParam<Real>("atomic_volume")),
    _kB(8.617343e-5), // eV/K
    // a_r Constants
    _znum(getParam<Real>("Z")),
    _a0(getParam<Real>("a_0")),
    _Di0(getParam<Real>("Di_0")),
    _EiB(getParam<Real>("Ei_B")),
    // Switching Functions
    _hv(getMaterialProperty<Real>("hv")),
    _hs(getMaterialProperty<Real>("hs")),
    _dhs(getMaterialPropertyDerivativeByName<Real>("hs", _phi_name))

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
UO2RecombinationMaterial::computeQpProperties()
{
  // Interstitial Bulk D
  Real Dib = _Di0 * std::exp(-_EiB / _kB / 1600);
  // Recombination rate / prefactor term
  Real ar = _hs[_qp] * _Va * _znum * Dib / (_a0 * _a0);
  Real dar_dphi = _dhs[_qp] * _Va * _znum * Dib / (_a0 * _a0);

  // Recombination
  Real rec = ar * _cu[_qp] * _ci[_qp] / _Va;

  // Set recombination if > 0 and hv is below the threshold value specified
  if ((_hv[_qp] < _hv_thresh) && (rec > 0.0))
  {
    // Negative recombination rate
    _rec[_qp] = -rec;
    // remember from GPMultiSinteringMat chi in masscons is chi*Va
    _drecdphi[_qp] = -dar_dphi * _cu[_qp] * _ci[_qp] / _Va;
    _drecdwu[_qp] = -(_chiu[_qp] / _Va) * ar * _ci[_qp] / _Va;
    _drecdwi[_qp] = -(_chii[_qp] / _Va) * ar * _cu[_qp] / _Va;
    _drecdcu[_qp] = -ar * _ci[_qp] / _Va;
    _drecdci[_qp] = -ar * _cu[_qp] / _Va;
  }
  else
  {
    // Set recombination rate and derivatives to zero
    _rec[_qp] = 0.0;
    _drecdphi[_qp] = 0.0;
    _drecdwu[_qp] = 0.0;
    _drecdwi[_qp] = 0.0;
    _drecdcu[_qp] = 0.0;
    _drecdci[_qp] = 0.0;
  }
}
