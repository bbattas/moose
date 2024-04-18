#pragma once

// #include "PolycrystalDiffusivityTensorBase.h"
// #include "DerivativeMaterialPropertyNameInterface.h"  //i think depreciated to
// DerivativeMaterialInterface???
#include "Material.h"
#include "DerivativeMaterialInterface.h"

/**
 * Calculates mobilities for grand potential model. The potential mobility (\chi*D)
 * is a scalar, while the Allen Cahn mobilities for the solid and void phases are
 * also scalars.
 */
// class GrandPotentialIsoMaterial : public PolycrystalDiffusivityTensorBase
class GrandPotentialIsoMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  GrandPotentialIsoMaterial(const InputParameters & parameters);

  //  virtual void computeProperties() override;

protected:
  // From PCDTB
  virtual void computeProperties();
  const VariableValue & _T;
  std::vector<const VariableValue *> _vals;
  const VariableValue & _c;
  VariableName _c_name;

  MaterialProperty<Real> & _D;
  MaterialProperty<Real> & _dDdc;

  Real _D0;
  Real _Em;
  Real _s_index;
  Real _gb_index;
  Real _b_index;
  Real _vap_index;
  Real _Dbulk;

  const Real _kb;
  const unsigned int _op_num;

  // From GPTensorMaterial
  /// mobility tensor
  std::string _D_name;
  MaterialProperty<Real> & _chiD;
  MaterialProperty<Real> & _dchiDdc;

  /// grain boundary mobility
  std::string _Ls_name;
  MaterialProperty<Real> & _Ls;

  /// void mobility
  std::string _Lv_name;
  MaterialProperty<Real> & _Lv;

  // Actual grain boundary width
  Real _GBwidth;

  // Surface diffusion layer thickness
  Real _surf_thickness;

  /// magnitude of mobility tensor
  // MaterialProperty<Real> & _Dmag;

  /// surface energy
  // const MaterialProperty<Real> & _sigma_s;

  /// interface width
  Real _int_width;

  /// susceptibility
  const MaterialPropertyName _chi_name;
  const MaterialProperty<Real> & _chi;
  const MaterialProperty<Real> & _dchidc;
  std::vector<const MaterialProperty<Real> *> _dchideta;
  std::vector<MaterialProperty<Real> *> _dchiDdeta;

  Real _GBMobility;
  Real _GBmob0;
  const Real _Q;

  /// solid phase order parameters
  std::vector<NonlinearVariableName> _vals_name;

  // Test output for diffusivity
  std::string _D_out_name;
  MaterialProperty<Real> & _D_out;

  // Toggle to for the iw based scaling
  bool _iw_scaling_bool;

  // Toggle to use interstitial values for GB and surface D
  bool _interstitials;

  // const MooseEnum _iw_scaling;
};
