#include "GrandPotentialIsoMaterial.h"

// libMesh includes
#include "libmesh/quadrature.h"

registerMooseObject("PhaseFieldApp", GrandPotentialIsoMaterial);

InputParameters
GrandPotentialIsoMaterial::validParams()
{
  // InputParameters params = PolycrystalDiffusivityTensorBase::validParams();
  InputParameters params = Material::validParams();
  params.addClassDescription("Diffusion and mobility parameters for grand potential model "
                             "governing equations. Uses a scalar diffusivity");
  // From PCDTB
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables"); //
  params.addCoupledVar("T", "Temperature variable in Kelvin");
  params.addRequiredParam<Real>(
      "D0", "Diffusion prefactor for vacancies in m^2/s"); // paper value is nm2/s??????
  params.addRequiredParam<Real>("Em", "Vacancy migration energy in eV");
  params.addRequiredCoupledVar("c", "Vacancy phase variable");
  params.addParam<Real>("surfindex", 1.0, "Surface diffusion index weight");
  params.addParam<Real>("gbindex",
                        -1.0,
                        "Grain boundary diffusion index weight"
                        "(-1 uses LANL GB D)");
  params.addParam<Real>("bulkindex", 1.0, "Bulk diffusion index weight");
  params.addParam<Real>("vaporindex",
                        1e-2,
                        "Vapor transport (void phase) diffusion weight"
                        "to multiply the bulk D by.");
  // from GPTensorMaterial
  params.addRequiredParam<Real>("int_width",
                                "The interfacial width in the lengthscale of the problem");
  params.addParam<MaterialPropertyName>("chi", "Coefficient to multiply by D");
  params.addParam<Real>("GBmob0", 0.0, "Grain boundary mobility prefactor");
  params.addRequiredParam<Real>("Q", "Grain boundary migration activation energy in eV");
  params.addParam<Real>(
      "GBMobility", -1, "GB mobility input that overrides the temperature dependent calculation");
  params.addParam<std::string>("f_name", "chiD", "Name for the mobility material property");
  // params.addRequiredParam<MaterialPropertyName>("surface_energy", "Surface energy of material");
  params.addParam<std::string>("solid_mobility", "L", "Name of grain mobility for solid phase");
  params.addParam<std::string>("void_mobility", "Lv", "Name of void phase mobility");
  params.addParam<Real>("GBwidth",
                        1.0,
                        "Real grain boundary width in units"
                        "of problem (length -> nm)");
  params.addParam<Real>("surf_thickness",
                        0.5,
                        "Surface diffusion layer thickness"
                        "in units of problem (length -> nm)");
  params.addParam<bool>("iw_scaling", true, "Enable the iw based scaling for GB and Surface D.");
  // ON OFF for D Scaling
  // MooseEnum iw_scaling("TRUE FALSE", "TRUE"); //ADDED
  // params.addParam<MooseEnum>("iw_scaling", iw_scaling,
  //                             "Whether or not to scale D based on IW."); //ADDED
  return params;
}

GrandPotentialIsoMaterial::GrandPotentialIsoMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    //  : PolycrystalDiffusivityTensorBase(parameters),
    // From PCDTB  (all but grad_c)
    _T(coupledValue("T")),
    _c(coupledValue("c")),
    _c_name(getVar("c", 0)->name()),
    _D(declareProperty<Real>("D")),
    _dDdc(declarePropertyDerivative<Real>("D", _c_name)),
    _D0(getParam<Real>("D0")),
    _Em(getParam<Real>("Em")),
    _s_index(getParam<Real>("surfindex")),
    _gb_index(getParam<Real>("gbindex")),
    _b_index(getParam<Real>("bulkindex")),
    _vap_index(getParam<Real>("vaporindex")),
    _kb(8.617343e-5), // Boltzmann constant in eV/K
    _op_num(coupledComponents("v")),
    // From GPTensorMaterial
    _D_name(getParam<std::string>("f_name")),
    _chiD(declareProperty<Real>(_D_name)),
    _dchiDdc(declarePropertyDerivative<Real>(_D_name, _c_name)),
    _Ls_name(getParam<std::string>("solid_mobility")),
    _Ls(declareProperty<Real>(_Ls_name)),
    _Lv_name(getParam<std::string>("void_mobility")),
    _Lv(declareProperty<Real>(_Lv_name)),
    _GBwidth(getParam<Real>("GBwidth")),
    _surf_thickness(getParam<Real>("surf_thickness")),
    // _sigma_s(getMaterialProperty<Real>("surface_energy")),
    _int_width(getParam<Real>("int_width")),
    _chi_name(getParam<MaterialPropertyName>("chi")),
    _chi(getMaterialProperty<Real>(_chi_name)),
    _dchidc(getMaterialPropertyDerivative<Real>("chi", _c_name)),
    _dchideta(_op_num),
    _dchiDdeta(_op_num),
    _GBMobility(getParam<Real>("GBMobility")),
    _GBmob0(getParam<Real>("GBmob0")),
    _Q(getParam<Real>("Q")),
    _vals_name(_op_num),
    _D_out(declareProperty<Real>("diffusivity")),
    _iw_scaling_bool(getParam<bool>("iw_scaling"))
// _iw_scaling(getParam<MooseEnum>("iw_scaling")) //ADDED
{
  if (_op_num == 0)                          //
    mooseError("Model requires op_num > 0"); //

  _vals.resize(_op_num); //
  for (unsigned int i = 0; i < _op_num; ++i)
  {
    _vals_name[i] = getVar("v", i)->name();
    _vals[i] = &coupledValue("v", i); //?????? for DGB
    _dchideta[i] = &getMaterialPropertyDerivative<Real>(_chi_name, _vals_name[i]);
    _dchiDdeta[i] = &declarePropertyDerivative<Real>(_D_name, _vals_name[i]);
  }
}

void
GrandPotentialIsoMaterial::computeProperties()
{
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    Real c = _c[_qp];  // phi -- void phase
    Real mc = 1.0 - c; // 1 - phi

    // Calculate Switching Functions
    // Make hgb -- gb switching function
    // Maintaining his style of calculation so it remains consistent to the
    // tensor form
    Real hgb = 0.0;
    for (unsigned int i = 0; i < _op_num; ++i)
      for (unsigned int j = i + 1; j < _op_num; ++j)
      {
        hgb += (*_vals[i])[_qp] * (*_vals[j])[_qp];
        hgb += (*_vals[j])[_qp] * (*_vals[i])[_qp];
      }
    hgb = 8 * hgb; // end result is 16 * gr1^2 * gr2^2

    // Make hsurf -- free surface switching function
    Real hsurf = 16 * c * c * mc * mc;
    Real dhsurfdc = 32.0 * c * mc * mc - 32.0 * c * c * mc;

    // void phase switching function
    Real hv = c * c * c * (10.0 + c * (-15.0 + c * 6.0));
    Real dhv = 30.0 * c * c * (c - 1.0) * (c - 1.0);

    // Caclulate Diffusivities
    // _Dbulk -- bulk diffusivity
    _Dbulk = _D0 * std::exp(-_Em / _kb / _T[_qp]); // * _b_index
    // Real dDbulkdc = 0;  //bulk diffusivity is independent of the void phase
    // vapor transport diffusivity -- this is a poor approximation -- temporary
    // Possible it might cause convergence issues if its too low
    Real Dv = _Dbulk * _vap_index; // 1E-2;

    // Dgb -- grain boundary diffusivity
    // Real Dgb = _Dbulk * _gb_index;
    Real Dgb = 0;
    if (_gb_index == -1)
    {
      Dgb = 4.74E14 * std::exp(-2.72 / _kb / _T[_qp]);
    }
    else
    {
      Dgb = _Dbulk * _gb_index;
    }

    // Real Dgb = 4.74E14 * std::exp(-2.72 / _kb / _T[_qp]);  //Single line option without the
    // gb_index toggle
    // Real dDgbdc = 0;  //gb diffusivity is independent of the void phase so this
    // was commented out

    // Dsurf -- free surface diffusivity
    Real Dsurf = _Dbulk * _s_index;

    // Define the variables for interface diffusivities
    Real newDgb = 0.0;
    Real newDsurf = 0.0;
    if (_iw_scaling_bool) // Use the IW scaling factor
    {
      // Define Scaling Factor
      Real gbScale = 1.5 * _GBwidth / _int_width;
      Real surfScale = 1.5 * _surf_thickness / _int_width;
      // Apply scaling factor to Dgb and Ds
      newDgb = Dgb * gbScale + _Dbulk * (1 - gbScale);
      newDsurf = Dsurf * surfScale + _Dbulk * (1 - surfScale);
    }
    else // Dont use iw scaling- just normal iso D
    {
      // Set the variable values = the unscaled version
      newDgb = Dgb;
      newDsurf = Dsurf;
    }

    // Calculate the scalar D from the determined components
    // not final- check perpendicular boundary in tonks jupyter notebook for surface version?
    _D[_qp] = (1 - hgb - hsurf) * (hv * Dv + (1 - hv) * _Dbulk) + hgb * newDgb + hsurf * newDsurf;
    _dDdc[_qp] = (1 - hgb - hsurf) * (dhv * Dv - dhv * _Dbulk) -
                 dhsurfdc * (hv * Dv + (1 - hv) * _Dbulk) + dhsurfdc * newDsurf;

    // Chemical susceptibility * Diffusivity
    _chiD[_qp] = _D[_qp] * _chi[_qp];
    _dchiDdc[_qp] = _dDdc[_qp] * _chi[_qp] + _D[_qp] * _dchidc[_qp];
    for (unsigned int i = 0; i < _op_num; ++i)
      (*_dchiDdeta[i])[_qp] = _D[_qp] * (*_dchideta[i])[_qp];

    // _Dmag[_qp] = _chiD[_qp].norm();

    Real GBmob;
    if (_GBMobility < 0)
      GBmob = _GBmob0 * std::exp(-_Q / (_kb * _T[_qp]));
    else
      GBmob = _GBMobility;

    _Ls[_qp] = 4.0 / 3.0 * GBmob / _int_width;
    _Lv[_qp] = 40 * _Ls[_qp];
    // Test output for Diffusivity
    _D_out[_qp] = _D[_qp];
  }
}
