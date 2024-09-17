#include "UO2CeqMaterial.h"

// libMesh includes might not need with computeQpProp?
#include "libmesh/quadrature.h"

registerMooseObject("PhaseFieldApp", UO2CeqMaterial);

InputParameters
UO2CeqMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("UO2 equilibrium concentration with or without irradiation "
                             "for vacancies or interstitials (at 1600K)");
  params.addRequiredCoupledVarWithAutoBuild(
      "etas", "var_name_base", "op_num", "Array of coupled variables"); //
  // params.addCoupledVar("T", "Temperature variable in Kelvin");
  // params.addRequiredCoupledVar("void_phase", "Vacancy/void phase variable");
  params.addParam<std::string>(
      "ceq_name", "cv_eq", "Name for the solid phase eq conc to be saved as");
  params.addParam<MaterialPropertyName>("hgb", "hgb", "Name of GB switching function material");
  MooseEnum vac_int("VAC_TH INT_TH VAC_IRR INT_IRR", "VAC_TH");
  params.addRequiredParam<MooseEnum>(
      "vi", vac_int, "Vacancies or interstitials with or without irradiaton");
  return params;
}

UO2CeqMaterial::UO2CeqMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    // Grain OPs
    _neta(coupledComponents("etas")),
    _eta(_neta),
    _eta_name(_neta),
    // _op_num(coupledComponents("v")), // Number of grains
    // _vals_name(_op_num),
    // _T(coupledValue("T")),
    // _c(coupledValue("void_phase")),
    // _c_name(getVar("void_phase", 0)->name()),
    // Output Ceq Name
    _ceq_name(getParam<std::string>("ceq_name")),
    _ceqs(declareProperty<Real>(_ceq_name)),
    _dceqsdeta(_neta),
    _d2ceqsdeta2(_neta),
    // hgb
    _hgb_name(getParam<MaterialPropertyName>("hgb")),
    _hgb(getMaterialProperty<Real>(_hgb_name)),
    _dhgbdeta(_neta),
    _d2hgbdeta2(_neta),
    // Toggle
    _vac_int(getParam<MooseEnum>("vi"))

{
  if (_neta == 0)
    mooseError("Model requires op_num > 0 for ceq material");

  for (unsigned int i = 0; i < _neta; ++i) //_op_num
  {
    _eta[i] = &coupledValue("etas", i);    // Grain OP Values
    _eta_name[i] = coupledName("etas", i); // Grain OP Names
    // First Derivatives
    _dhgbdeta[i] = &getMaterialPropertyDerivativeByName<Real>(_hgb_name, _eta_name[i]);
    _dceqsdeta[i] = &declarePropertyDerivative<Real>(_ceq_name, _eta_name[i]);
    // No Cross term second derivatives
    _d2hgbdeta2[i].resize(_neta);
    _d2hgbdeta2[i][i] =
        &getMaterialPropertyDerivativeByName<Real>(_hgb_name, _eta_name[i], _eta_name[i]);
    _d2ceqsdeta2[i].resize(_neta);
    _d2ceqsdeta2[i][i] = &declarePropertyDerivative<Real>(_ceq_name, _eta_name[i], _eta_name[i]);

    // Commented out getting and declaring second deriv cross terms
    // _d2ceqsdeta2[i].resize(_op_num)
    // for (unsigned int j = 0; j <= i; ++j)
    // {
    //   _d2hgbdeta2[j][i] =
    //       &getMaterialPropertyDerivativeByName<Real>(_csu_eq_name, _vals_name[j], _vals_name[i]);
    //   _d2ceqsdeta2[j][i] = &declarePropertyDerivative<Real>(_ceq_name, _vals_name[j],
    //   _vals_name[i]);
    // }
  }
}

void
UO2CeqMaterial::computeQpProperties()
{
  Real cb = 0.0;
  Real cgb = 0.0;
  // Switch based on vac/int thermal/irradiation
  switch (_vac_int)
  {
    case 0: // Vacancies Thermal
      cb = 2.424e-06;
      cgb = 5.130e-03;
      break;
    case 1: // Interstitials Thermal
      cb = 1.667e-32;
      cgb = 6.170e-08;
      break;
    case 2: // Vacancies Irradiation
      cb = 3.877e-04;
      cgb = 4.347e-03;
      break;
    case 3: // Interstitials Irradiation
      cb = 7.258e-09;
      cgb = 5.900e-06;
      break;
    default:
      paramError("vi", "Incorrect Vacancy or Interstitial");
  }
  // Solid phase EQ Conc
  // _ceqs[_qp] = cgb * _hgb[_qp] + (1 - _hgb[_qp]) * cb;
  _ceqs[_qp] = _hgb[_qp] * (cgb - cb) + cb; // Rearranged for easy derivs
  // Grain OP Derivatives
  for (unsigned int i = 0; i < _neta; ++i)
  {
    (*_dceqsdeta[i])[_qp] = (*_dhgbdeta[i])[_qp] * (cgb - cb) + cb;
    (*_d2ceqsdeta2[i][i])[_qp] = (*_d2hgbdeta2[i][i])[_qp] * (cgb - cb) + cb;
  }
}
