#pragma once

#include "Material.h"

/**
 * A material that couples a material property
 */
class GrandPotentialInterstitialMaterial : public Material
{
public:
  static InputParameters validParams();

  GrandPotentialInterstitialMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override { _mat_prop[_qp] = 0.0; } // was 1.0
  virtual void computeQpProperties() override;

  std::string _mat_prop_name;
  MaterialProperty<Real> & _mat_prop;

  const MaterialProperty<Real> & _coupled_mat_prop;
  const MaterialProperty<Real> & _rhov;
  const MaterialProperty<Real> & _rhos;
  const MaterialProperty<Real> & _hv;
  const MaterialProperty<Real> & _hs;
  const MaterialProperty<Real> & _self_old;
};
