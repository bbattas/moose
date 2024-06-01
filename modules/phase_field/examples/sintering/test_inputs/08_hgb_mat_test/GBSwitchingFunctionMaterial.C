//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "GBSwitchingFunctionMaterial.h"

registerMooseObject("PhaseFieldApp", GBSwitchingFunctionMaterial);

InputParameters
GBSwitchingFunctionMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Helper material to provide h(eta) and its derivative in one of two "
                             "polynomial forms. MIX234 and MIX246");
  params.addParam<std::string>("property_name", "hgb", "actual name for f(eta), i.e. 'h' or 'g'");
  return params;
}

GBSwitchingFunctionMaterial::GBSwitchingFunctionMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    // Break line
    _op_num(coupledComponents("v")),
    _vals_name(_op_num),
    _hgb_name(getParam<std::string>("property_name")),
    _hgb(declareProperty<Real>(_hgb_name)),
    _dhgbdeta(_op_num)
{
  if (_op_num == 0)                          //
    mooseError("Model requires op_num > 0"); //

  _vals.resize(_op_num); //
  for (unsigned int i = 0; i < _op_num; ++i)
  {
    _vals_name[i] = getVar("v", i)->name();
    _vals[i] = &coupledValue("v", i);
    _dhgbdeta[i] = &declarePropertyDerivative<Real>(_hgb_name, _vals_name[i]);
  }
}

void
GBSwitchingFunctionMaterial::computeQpProperties()
{
  Real hgb = 0;
}
