#include "GrandPotentialInterstitialMaterial.h"

registerMooseObject("PhaseFieldApp", GrandPotentialInterstitialMaterial);

InputParameters
GrandPotentialInterstitialMaterial::validParams()
{
  InputParameters params = Material::validParams();
  // params.addParam<MaterialPropertyName>("gen_mat", "Name of parsed material with the vacancy
  // generation term");
  params.addRequiredParam<MaterialPropertyName>("mat_prop",
                                                "Name of the property this material defines");
  params.addRequiredParam<MaterialPropertyName>(
      "combined_vac_mat", "Name of the vacancy material property to couple into this material");
  params.addParam<bool>(
      "use_old_prop",
      true,
      "Boolean indicating whether to use the old coupled property instead of the current property");
  return params;
}

GrandPotentialInterstitialMaterial::GrandPotentialInterstitialMaterial(
    const InputParameters & parameters)
  : Material(parameters),
    _mat_prop_name(getParam<MaterialPropertyName>("mat_prop")),
    _mat_prop(declareProperty<Real>(_mat_prop_name)),
    _coupled_mat_prop(getParam<bool>("use_old_prop")
                          ? getMaterialPropertyOld<Real>("combined_vac_mat")
                          : getMaterialProperty<Real>("combined_vac_mat")),
    _rhov(getMaterialProperty<Real>("rhov")),
    _rhos(getMaterialProperty<Real>("rhos")),
    _hv(getMaterialProperty<Real>("hv")),
    _hs(getMaterialProperty<Real>("hs")),
    _self_old(getMaterialPropertyOld<Real>(_mat_prop_name))

{
}

void
GrandPotentialInterstitialMaterial::computeQpProperties()
{
  // Current vacancy density
  Real vacancies = _hs[_qp] * _rhos[_qp] + _hv[_qp] * _rhov[_qp];
  // Vacancy density from previous timestep
  Real old_vac = _coupled_mat_prop[_qp];
  // if the previous timestep is nan use the current instead
  if (std::isnan(_coupled_mat_prop[_qp]))
    old_vac = vacancies;
  // Change in vacancy concentration from previous to current timestep
  Real vac_change = vacancies - old_vac;

  // Pull the previous timestep value for this material (0.0 if nan)
  Real self_old = _self_old[_qp];
  if (std::isnan(_self_old[_qp]))
    self_old = 0.0; // vacancies;
  // Limit to positive total (no negative interstitial density)
  Real output_val = self_old + vac_change;
  if ((output_val < 0) || std::isnan(output_val))
    output_val = 0.0;
  _mat_prop[_qp] = output_val;
  // 4.0 / _coupled_mat_prop[_qp]; // This will produce a NaN if evaluated out of order
}
