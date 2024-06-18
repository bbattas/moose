#include "EBSDConsVarIC.h"

registerMooseObject("PhaseFieldApp", EBSDConsVarIC);

InputParameters
EBSDConsVarIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addClassDescription("IC based on ReconPhaseVarIC to read additional column(s) in the EBSD "
                             "data as a nodal variable IC.");
  params.addRequiredParam<UserObjectName>("ebsd_reader",
                                          "The EBSDReader object holding the EBSD data");
  params.addRequiredParam<unsigned int>("column", "Extra column idx (0,1,2...)");
  return params;
}

EBSDConsVarIC::EBSDConsVarIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _mesh(_fe_problem.mesh()),
    _ebsd_reader(getUserObject<EBSDReader>("ebsd_reader")),
    _cons(getParam<unsigned int>("column")),
    _node_to_custom_weight_map(_ebsd_reader.getNodeToCustomWeightMap())
{
}

Real
EBSDConsVarIC::value(const Point & /*p*/)
{
  // Return error if current node is NULL
  if (_current_node == nullptr)
    mooseError("_current_node is reporting NULL");

  // Make sure the _current_node is in the _node_to_phase_weight_map (return error if not in map)
  std::map<dof_id_type, std::vector<Real>>::const_iterator it =
      _node_to_custom_weight_map.find(_current_node->id());
  if (it == _node_to_custom_weight_map.end())
    mooseError("The following node id is not in the node map: ", _current_node->id());

  // make sure we have enough ophase weights
  if (_cons >= it->second.size())
    mooseError("Requested an out-of-range column number");

  return it->second[_cons];
}
