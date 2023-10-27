/****************************************************************/
/*                  DO NOT MODIFY THIS HEADER                   */
/*                           Marmot                             */
/*                                                              */
/*            (c) 2017 Battelle Energy Alliance, LLC            */
/*                     ALL RIGHTS RESERVED                      */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*             Under Contract No. DE-AC07-05ID14517             */
/*             With the U. S. Department of Energy              */
/*                                                              */
/*             See COPYRIGHT for full restrictions              */
/****************************************************************/

#ifndef GRANDPOTENTIALRADIATIONDAMAGE_H
#define GRANDPOTENTIALRADIATIONDAMAGE_H

#include "Material.h"
#include "DerivativeMaterialInterface.h"

class GrandPotentialRadiationDamage;

template<>
InputParameters validParams<GrandPotentialRadiationDamage>();

/**
 * This material calculates the following terms and 1st and 2nd derivatives for
 * the grand potential sintering model with an arbitrary number of order
 * parameters and chemical potentials:
 * - switching functions (h)
 * - susceptibilities (chi)
 * - number densities (rho)
 * - potential density functions (omega)
 * - energy barrier coefficient (mu)
 * - gradient energy coefficient (kappa)
 * - interface symmetry coefficient (gamma = 1.5)
 */

class GrandPotentialRadiationDamage : public DerivativeMaterialInterface<Material>
{
public:
  GrandPotentialRadiationDamage(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  /// solid-phase order parameters
  const unsigned int _nop;
  std::vector<const VariableValue *> _ops;
  std::vector<VariableName> _op_names;
  /// void-phase order parameter
  const VariableValue & _phi;
  const NonlinearVariableName _phi_name;
  /// chemical potentials
  const unsigned int _nw;
  std::vector<const VariableValue *> _w;
  std::vector<VariableName> _w_names;
  /// Temperature
  const VariableValue & _T;

  /// Naming conventions
  const std::vector<std::string> _sub_w;
  const std::string _sub_s;
  const std::string _sub_v;
  const std::string _h_base;
  const std::string _chi_base;
  const std::string _rho_base;
  const std::string _omega_base;
  const std::string _kappa_name;
  const std::string _gamma_name;
  /// Parabolic energy coefficients
  const std::vector<MaterialPropertyName> _ksolid_names;
  std::vector<const MaterialProperty<Real> *> _ksolid;
  const std::vector<MaterialPropertyName> _kvoid_names;
  std::vector<const MaterialProperty<Real> *> _kvoid;
  /// Equilibrium concentrations
  const std::vector<MaterialPropertyName> _cs_eq_names;
  std::vector<const MaterialProperty<Real> *> _cs_eq;
  std::vector<std::vector<const MaterialProperty<Real> *>> _dcs_eqdeta;
  std::vector<std::vector<std::vector<const MaterialProperty<Real> *>>> _d2cs_eqdeta2;
  /// Solid-phase energy function. Three options available according to Plapp (2011).
  const MooseEnum _solid_energy;
  /// Interface energies
  const std::vector<Real> _cv_eq;
  const Real _sigma_s;
  const Real _sigma_gb;
  /// Other constants
  const Real _int_width;
  const Real _phi0;
  const Real _Va;

  /// Switching Functions
  MaterialProperty<Real> & _hs;
  MaterialProperty<Real> & _dhsdphi;
  MaterialProperty<Real> & _d2hsdphi2;
  MaterialProperty<Real> & _hv;
  MaterialProperty<Real> & _dhvdphi;
  MaterialProperty<Real> & _d2hvdphi2;
  /// Susceptibilities
  std::vector<MaterialProperty<Real> *> _chi;
  std::vector<MaterialProperty<Real> *> _dchidw;
  std::vector<MaterialProperty<Real> *> _dchidphi;
  std::vector<MaterialProperty<Real> *> _d2chidw2;
  std::vector<MaterialProperty<Real> *> _d2chidphidw;
  std::vector<MaterialProperty<Real> *> _d2chidphi2;
  std::vector<std::vector<MaterialProperty<Real> *>> _dchideta;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2chidetadphi;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2chidetadw;
  std::vector<std::vector<std::vector<MaterialProperty<Real> *>>> _d2chideta2;
  /// Densities
  std::vector<MaterialProperty<Real> *> _rhos;
  std::vector<MaterialProperty<Real> *> _drhosdw;
  std::vector<MaterialProperty<Real> *> _d2rhosdw2;
  std::vector<std::vector<MaterialProperty<Real> *>> _drhosdeta;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2rhosdetadw;
  std::vector<std::vector<std::vector<MaterialProperty<Real> *>>> _d2rhosdeta2;
  std::vector<MaterialProperty<Real> *> _rhov;
  std::vector<MaterialProperty<Real> *> _drhovdw;
  // Potential Densities
  MaterialProperty<Real> & _omegas;
  std::vector<MaterialProperty<Real> *> _domegasdw;
  std::vector<MaterialProperty<Real> *> _d2omegasdw2;
  std::vector<MaterialProperty<Real> *> _domegasdeta;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2omegasdetadw;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2omegasdeta2;
  MaterialProperty<Real> & _omegav;
  std::vector<MaterialProperty<Real> *> _domegavdw;
  std::vector<MaterialProperty<Real> *> _d2omegavdw2;
  // Energy coefficients
  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> & _dmudphi;
  MaterialProperty<Real> & _d2mudphi2;
  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _dkappadphi;
  MaterialProperty<Real> & _d2kappadphi2;
  MaterialProperty<Real> & _gamma;

private:
  std::vector<Real> switchingBase(Real x, Real x0);
};

#endif //GRANDPOTENTIALRADIATIONDAMAGE_H
