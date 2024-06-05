##############################################################################
# File: 00_small2D_d3d_viz.i
# File Location: /examples/sintering/paper3/18_newCeq/00_d3d_txt
# Created Date: Wednesday June 5th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Wednesday June 5th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Testing the new smaller but still largeish 2D ic for testing
#  310,947 DoFs
#
#
##############################################################################

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = 2D_50x50um_8umavg_5pore.txt
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = ebsd_mesh
    combinatorial_geometry = 'x > 48000'
    block_id = 1
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 6
  var_name_base = gr
  int_width = 2000 #1000 #min radius is like 2250, element size of 250
  # profile = TANH # not used at the moment? only in circleic?
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
    []
  []
  [VoidIC] #External
    type = BoundingBoxIC
    variable = phi
    block = 1
    inside = 1
    outside = 0.01
    x1 = 50000
    x2 = 80000
    y1 = -5000
    y2 = 55000
  []
  [voidIC2] #Internal
    type = SpecifiedSmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0.01
    radii = '2843.2  2952.34 3306.29 3235.79 3089.41'
    x_positions = '17637.6  17065.8  38663.37 30171.23 29699.3'
    y_positions = '13307.48 34958.57 19841.07 36142.58  9987.43'
    z_positions = '    0.       0.       0.       0.       0.  '
    block = 0
  []
[]

[Variables]
  [phi]
  []
  [PolycrystalVariables]
  []
[]

[AuxVariables]
  [bnds]
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
  []
[]

[UserObjects]
  [ebsd_reader]
    type = EBSDReader
  []
  [ebsd]
    type = PolycrystalEBSD
    # coloring_algorithm = bt
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    output_adjacency_matrix = true
    phase = 1
  []
  [grain_tracker]
    type = GrainTracker
    threshold = 0.1 #0.2
    connecting_threshold = 0.09 #0.08
    compute_halo_maps = false #true#false
    verbosity_level = 1
  []
[]

[Problem]
  type = FEProblem
  solve = false
  kernel_coverage_check = false
[]

[Executioner]
  type = Transient
  num_steps = 0
[]

[Outputs]
  checkpoint = false
  exodus = false
  nemesis = true
[]
