#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"
#include "RankThreeTensor.h"

/**
 * Generates a diffusion tensor to distinguish between the bulk, grain boundary,
 * and surface diffusion rates.
 */
class ElementalGammaWellMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ElementalGammaWellMaterial(const InputParameters & parameters);

protected:
  virtual void computeProperties();

  const MaterialProperty<Real> & _gamma_in;
  MaterialProperty<Real> & _gamma_out;

  const MaterialProperty<Real> & _int_width_in;
  MaterialProperty<Real> & _int_width_out;

  const MaterialProperty<Real> & _L_in;
  MaterialProperty<Real> & _L_out;

  const bool _skip;
};
