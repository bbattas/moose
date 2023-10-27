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

#ifndef GRANDPOTENTIALSCALARMATERIAL_H
#define GRANDPOTENTIALSCALARMATERIAL_H

#include "Material.h"
#include "DerivativeMaterialInterface.h"

class GrandPotentialScalarMaterial;

template<>
InputParameters validParams<GrandPotentialScalarMaterial>();

class GrandPotentialScalarMaterial : public DerivativeMaterialInterface<Material>
{
public:
  GrandPotentialScalarMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  const unsigned int _nop;
  std::vector<const VariableValue *> _ops;
  std::vector<VariableName> _op_names;

  const VariableValue & _phi;
  const NonlinearVariableName _phi_name;

  const NonlinearVariableName _w_name;

  const VariableValue & _T;

  const Real _D0;
  const Real _Em;
  const Real _wB;
  const Real _wGB;
  const Real _wS;

  const MaterialPropertyName _chi_name;
  const MaterialProperty<Real> & _chi;
  const MaterialProperty<Real> & _dchidw;
  const MaterialProperty<Real> & _dchidphi;
  const MaterialProperty<Real> & _d2chidw2;
  const MaterialProperty<Real> & _d2chidphidw;
  const MaterialProperty<Real> & _d2chidphi2;

  const std::string _Dname;
  MaterialProperty<Real> & _D;
  MaterialProperty<Real> & _dDdphi;
  std::vector<MaterialProperty<Real> *> _dDdeta;
  MaterialProperty<Real> & _d2Ddphi2;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2Ddeta2;

  const std::string _Mname;
  MaterialProperty<Real> & _M;
  MaterialProperty<Real> & _dMdw;
  MaterialProperty<Real> & _dMdphi;
  std::vector<MaterialProperty<Real> *> _dMdeta;
  MaterialProperty<Real> & _d2Mdw2;
  MaterialProperty<Real> & _d2Mdphidw;
  std::vector<MaterialProperty<Real> *> _d2Mdetadw;
  MaterialProperty<Real> & _d2Mdphi2;
  std::vector<MaterialProperty<Real> *> _d2Mdetadphi;
  std::vector<std::vector<MaterialProperty<Real> *>> _d2Mdeta2;
};

#endif //GRANDPOTENTIALSCALARMATERIAL_H
