#pragma once

#include "DerivativeParsedMaterialHelper.h"
#include "ExpressionBuilder.h"

/**
 * Calculates the equilibrium solid phase vacancy concentration in UO2+/-x as a
 * function of temperature and the stoichiometric deviation (O/U ratio).  For
 * use with GrandPotentialSinteringMaterial- output (f_name) here is the
 * "equilibrium_vacancy_concentration" in GPSM
 */
class UO2CvMaterial : public DerivativeParsedMaterialHelper, public ExpressionBuilder
{
public:
  static InputParameters validParams();

  UO2CvMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  const VariableValue & _T;
  const VariableValue & _OU;
  const VariableValue & _phi;

  const Real _kb;

  const unsigned int _op_num;
  std::vector<const VariableValue *> _vals;
  std::vector<NonlinearVariableName> _vals_name;

  std::string _cv_name;
  MaterialProperty<Real> & _cv_out;

  std::vector<MaterialProperty<Real> *> _dcv;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2cv;

  const std::vector<FileName> _gb_csv;
  std::vector<std::vector<Real>> _gb_data;

  const Real _gb_max;
  const Real _b_max;

  // std::vector<Real> _sigma11_data;
  // std::vector<Real> _sigma9_data;

};
