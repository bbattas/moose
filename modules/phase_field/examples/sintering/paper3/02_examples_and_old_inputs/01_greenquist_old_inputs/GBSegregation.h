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

#ifndef GBSEGREGATION_H
#define GBSEGREGATION_H

#include "Material.h"
#include "DerivativeMaterialInterface.h"

class GBSegregation;

template<>
InputParameters validParams<GBSegregation>();

class GBSegregation : public DerivativeMaterialInterface<Material>
{
public:
  GBSegregation(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  const Real _A;
  const Real _Ef;
  const Real _Egb;
  const Real _max_bulk;
  const Real _max_gb;
  const VariableValue & _T;
  const unsigned int _nop;
  std::vector<const VariableValue *> _ops;
  std::vector<VariableName> _op_names;

  const std::string _ceq_name;
  MaterialProperty<Real> & _ceq;
  std::vector<MaterialProperty<Real> *> _dceqdopa;
  std::vector<MaterialProperty<Real> *> _d2ceqdopa2;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2ceqdopadopb;
};

#endif // GBSEGREGATION_H
