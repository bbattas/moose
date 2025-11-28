#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"
#include "RankThreeTensor.h"

/**
 * Generates a diffusion tensor to distinguish between the bulk, grain boundary,
 * and surface diffusion rates.
 */
class ElementalGammaMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  ElementalGammaMaterial(const InputParameters & parameters);

protected:
  virtual void computeProperties();

  const MaterialProperty<Real> & _gamma_in;
  MaterialProperty<Real> & _gamma_out;

  const MaterialProperty<bool> & _no_ij_pairs;
  MaterialProperty<bool> & _elem_no_ij;
  MaterialProperty<Real> & _int_noij;

  const MaterialProperty<Real> & _gbe_iso;
  const MaterialProperty<Real> & _kappa;
  Real _const_m;

  const bool _aniso_L;
  const MaterialProperty<Real> & _L0;
  const MaterialProperty<Real> * _L_qp;
  MaterialProperty<Real> * _L_elem;

  const bool _well;
  const MaterialProperty<Real> & _int_width_in;
  MaterialProperty<Real> & _int_width_out;
};
