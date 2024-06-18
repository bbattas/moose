#pragma once

#include "InitialCondition.h"
#include "EBSDReader.h"
#include "PolycrystalICTools.h"

/**
 * IC based on ReconPhaseVarIC to read additional column(s) in the EBSD
 * data as a nodal variable IC.
 */
class EBSDConsVarIC : public InitialCondition
{
public:
  static InputParameters validParams();

  EBSDConsVarIC(const InputParameters & parameters);

  virtual Real value(const Point & /*p*/);

private:
  MooseMesh & _mesh;
  const EBSDReader & _ebsd_reader;
  unsigned int _cons;
  const std::map<dof_id_type, std::vector<Real>> & _node_to_custom_weight_map;
};
