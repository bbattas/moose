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
  // ENUM to set which version of the recombination i want
  MooseEnum thresh_case("FULL_IF SHORTSWITCH_NONEGIF SHORTSWITCH NONE", "FULL_IF");
  params.addRequiredParam<MooseEnum>(
      "if_case",
      thresh_case,
      "Which combination of if statements to use (if any), between the hs and"
      " keeping the recombination from being negative (growth).");
  params.addParam<Real>("phi_0",
                        0.3,
                        "Shortened switching function divisor (0,1]. Smaller sharpens interface. "
                        "Only used if thresh_case = SHORTSWITCH_NONEGIF or SHORTSWITCH.");
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
    _dhs(getMaterialPropertyDerivativeByName<Real>("hs", _phi_name)),
    // Toggle
    _if_case(getParam<MooseEnum>("if_case")),
    _switch(getParam<Real>("phi_0"))

{
  if (_neta == 0)
    mooseError("Model requires op_num > 0");

  if ((_switch > 1.0) || (_switch <= 0.0))
    mooseError("UO2RecombinationMaterial: phi_0 should be between 0 and 1, but not = 0");

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
  Real ar = _Va * _znum * Dib / (_a0 * _a0);
  // Real dar_dphi = _dhs[_qp] * _Va * _znum * Dib / (_a0 * _a0);

  // Empty placeholders for the values in the case statement based on enum switch
  Real rec, drecdphi, drecdwu, drecdwi, drecdcu, drecdci;

  // Which if any if statements to use
  switch (_if_case)
  {
    case 0: // FULL_IF: Previous approach of if on both hv and the - recombination value
    {
      // Recombination
      rec = _hs[_qp] * ar * _cu[_qp] * _ci[_qp] / _Va;

      // Set recombination if > 0 and hv is below the threshold value specified
      if ((_hv[_qp] < _hv_thresh) && (rec > 0.0))
      {
        // Negative recombination rate
        rec = -rec;
        // remember from GPMultiSinteringMat chi in masscons is chi*Va
        drecdphi = -ar * _dhs[_qp] * _cu[_qp] * _ci[_qp] / _Va;
        drecdwu = -(_chiu[_qp] / _Va) * ar * _hs[_qp] * _ci[_qp] / _Va;
        drecdwi = -(_chii[_qp] / _Va) * ar * _hs[_qp] * _cu[_qp] / _Va;
        drecdcu = -ar * _hs[_qp] * _ci[_qp] / _Va;
        drecdci = -ar * _hs[_qp] * _cu[_qp] / _Va;
      }
      else
      {
        // Set recombination rate and derivatives to zero
        rec = 0.0;
        drecdphi = 0.0;
        drecdwu = 0.0;
        drecdwi = 0.0;
        drecdcu = 0.0;
        drecdci = 0.0;
      }
      break;
    }

    case 1: // SHORTSWITCH_NONEGIF: Using cutdown switching function with only negative rec if
    {
      // Narrow switching function
      Real phi = _phi[_qp] / _switch;
      Real f = 1.0;
      Real df = 0.0;
      if (phi >= 1.0) // f is hs as a function of phi
        f = 0.0;
      else if (phi > 0.0)
      {
        f = 1 - (phi * phi * phi * (10.0 + phi * (-15.0 + phi * 6.0)));
        df = -30.0 / _switch * phi * phi * (phi - 1.0) * (phi - 1.0);
      }

      // Recombination
      rec = f * ar * _cu[_qp] * _ci[_qp] / _Va;

      // Set recombination if > 0
      if (rec > 0.0)
      {
        // Negative recombination rate
        rec = -rec;
        // remember from GPMultiSinteringMat chi in masscons is chi*Va
        drecdphi = -ar * df * _cu[_qp] * _ci[_qp] / _Va;
        drecdwu = -(_chiu[_qp] / _Va) * ar * f * _ci[_qp] / _Va;
        drecdwi = -(_chii[_qp] / _Va) * ar * f * _cu[_qp] / _Va;
        drecdcu = -ar * f * _ci[_qp] / _Va;
        drecdci = -ar * f * _cu[_qp] / _Va;
      }
      else
      {
        // Set recombination rate and derivatives to zero
        rec = 0.0;
        drecdphi = 0.0;
        drecdwu = 0.0;
        drecdwi = 0.0;
        drecdcu = 0.0;
        drecdci = 0.0;
      }
      break;
    }

    case 2: // SHORTSWITCH: no ifs just the narrow hs
    {
      // Narrow switching function
      Real phi = _phi[_qp] / _switch;
      Real f = 1.0;
      Real df = 0.0;
      if (phi >= 1.0) // f is hs as a function of phi
        f = 0.0;
      else if (phi > 0.0)
      {
        f = 1 - (phi * phi * phi * (10.0 + phi * (-15.0 + phi * 6.0)));
        df = -30.0 / _switch * phi * phi * (phi - 1.0) * (phi - 1.0);
      }

      // Recombination
      rec = -f * ar * _cu[_qp] * _ci[_qp] / _Va;
      // remember from GPMultiSinteringMat chi in masscons is chi*Va
      drecdphi = -ar * df * _cu[_qp] * _ci[_qp] / _Va;
      drecdwu = -(_chiu[_qp] / _Va) * ar * f * _ci[_qp] / _Va;
      drecdwi = -(_chii[_qp] / _Va) * ar * f * _cu[_qp] / _Va;
      drecdcu = -ar * f * _ci[_qp] / _Va;
      drecdci = -ar * f * _cu[_qp] / _Va;
      break;
    }

    case 3: // NONE: no if or narrowed hs
    {
      // Recombination as is
      rec = -_hs[_qp] * ar * _cu[_qp] * _ci[_qp] / _Va;
      drecdphi = -ar * _dhs[_qp] * _cu[_qp] * _ci[_qp] / _Va;
      drecdwu = -(_chiu[_qp] / _Va) * ar * _hs[_qp] * _ci[_qp] / _Va;
      drecdwi = -(_chii[_qp] / _Va) * ar * _hs[_qp] * _cu[_qp] / _Va;
      drecdcu = -ar * _hs[_qp] * _ci[_qp] / _Va;
      drecdci = -ar * _hs[_qp] * _cu[_qp] / _Va;
      break;
    }

    default:
      paramError("if_case", "Incorrect if case enum");
  }

  // Set recombination rate and derivatives the real values from cases
  _rec[_qp] = rec;
  _drecdphi[_qp] = drecdphi;
  _drecdwu[_qp] = drecdwu;
  _drecdwi[_qp] = drecdwi;
  _drecdcu[_qp] = drecdcu;
  _drecdci[_qp] = drecdci;
}
