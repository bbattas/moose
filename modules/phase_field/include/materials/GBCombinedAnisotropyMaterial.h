#pragma once

#include "GBInclinationBase.h"
#include "GBMisorientationHelper.h"

class GBCombinedAnisotropyMaterial : public GBInclinationBase, protected GBMisorientationHelper
{
public:
  static InputParameters validParams();

  GBCombinedAnisotropyMaterial(const InputParameters & parameters);

protected:
  struct AngleFunctionResult
  {
    Real f = 0.0;
    Real df_dtheta = 0.0;
    Real d2f_dtheta2 = 0.0;
    Real df_dpolar = 0.0;
    Real d2f_dpolar2 = 0.0;
    Real d2f_dthetadpolar = 0.0;
  };

  virtual void computeQpProperties() override;

  enum GBMode
  {
    COS = 0,
    INC = 1,
    MISO = 2,
    FULL = 3,
    ISO = 4
  };

  enum AvgType
  {
    AVG = 0,
    WEIGHTED = 1,
  };

  AngleFunctionResult computeIsoGBE(Real bulk_multiplier) const;

  AngleFunctionResult
  computeCosineOnlyGBE(Real theta_inc, Real polar_inc, Real a, Real b, Real c) const;

  AngleFunctionResult computeInclinationGBE(Real theta_inc, Real polar_inc) const;

  AngleFunctionResult computeMisorientationGBE(const MisorientationData & miso) const;

  AngleFunctionResult computeFullGBE(Real theta_inc,
                                     Real polar_inc,
                                     const MisorientationData & miso,
                                     Real w_inc,
                                     Real w_miso) const;

protected:
  const int _gb_mode;
  const int _tj_mode;
  const bool _aniso_L;
  const bool _stiffness;
  const bool _aniso_mob;

  // Iso Inputs
  const MaterialProperty<Real> & _kappa;
  const MaterialProperty<Real> & _gbe_iso;
  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _L0;

  // Aniso Outputs
  MaterialProperty<std::vector<Real>> & _fgbe;
  MaterialProperty<Real> & _gbe_norm;
  MaterialProperty<Real> & _gbe_gb;
  MaterialProperty<Real> & _gb_notj_mask;
  MaterialProperty<Real> & _gb_tj_mask;
  // MaterialProperty<std::vector<Real>> & _L_ij; // Array version for kernel? Not being used?
  MaterialProperty<Real> & _L;

  // gamma_ij and derivatives flatpacked to a vector with GBPairPacking
  MaterialProperty<std::vector<Real>> & _gamma_ij;
  MaterialProperty<std::vector<RealGradient>> & _dgamma_dgradeta;
  MaterialProperty<std::vector<RealTensorValue>> & _d2gamma_dgradeta2;

  // Condensed Forms
  MaterialProperty<Real> & _gamma_asymm;
  MaterialProperty<Real> & _int_width;

  // MaterialProperty<std::vector<Real>> & _gamma_ij;

  const Real _iso_gbe;   // normalized iso gbe
  const Real _bulk_mult; // Bulk normalized gbe value to use on non-gbs
  const Real _w_inc;
  const Real _w_miso;

  // COS Inclination function constants
  const Real _if_a;
  const Real _if_rot;
  // const Real _if_b;
  // const Real _if_c;

  // TEMPORARY
  // MaterialProperty<bool> & _elem_noij;
  MaterialProperty<Real> & _testout1;
  MaterialProperty<Real> & _thetaout;
  MaterialProperty<Real> & _noij_out;

  /// Moelans L Derivatives
  std::vector<MaterialProperty<Real> *> _dL_deta;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2L_deta2;
  // d grad eta
  std::vector<MaterialPropertyName> _dL_dgradeta_name;
  std::vector<MaterialProperty<RealGradient> *> _dL_dgradeta;
  std::vector<MaterialPropertyName> _d2L_dgradeta2_name;
  std::vector<MaterialProperty<std::vector<RealTensorValue>> *> _d2L_dgradeta2;
  std::vector<std::vector<MaterialProperty<RealGradient> *>> _d2L_dgradetadeta;
};
