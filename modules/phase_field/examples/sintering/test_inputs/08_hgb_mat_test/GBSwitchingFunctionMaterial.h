//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"

// Forward Declarations

/**
 * Material class to provide the switching function \f$ h(\eta) \f$ for
 * the KKS system.
 *
 * \see KKSPhaseChemicalPotential
 * \see KKSCHBulk
 */
class GBSwitchingFunctionMaterial : public DerivativeMaterialInterface<Material>
{
public:
  static InputParameters validParams();

  GBSwitchingFunctionMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  const unsigned int _op_num;

  std::vector<const VariableValue *> _vals;
  /// solid phase order parameters
  std::vector<NonlinearVariableName> _vals_name;

  std::string _hgb_name;
  MaterialProperty<Real> & _hgb;
  std::vector<MaterialProperty<Real> *> _dhgbdeta;
};
