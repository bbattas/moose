##############################################################################
# File: 00_2D_d3d_viz.i
# File Location: /examples/sintering/paper3/18_newCeq/00_d3d_txt
# Created Date: Friday May 31st 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Friday May 31st 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Bare minimum input for visualization of the d3d input
#  Just this simple input has 1,350,167 DoFs
#
#
##############################################################################

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = 2D_100x100um_8umavg_20pore.txt
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = ebsd_mesh
    combinatorial_geometry = 'x > 99000'
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
    x1 = 100000
    x2 = 150000
    y1 = -5000
    y2 = 125000
  []
  [voidIC2] #Internal
    type = SpecifiedSmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0.01
    radii = '3387.41 3204.63 3097.33 3157.62 2608.26 2573.66 3162.96 2580.68 2922.15 3499.
    3438.79 3442.07 3025.42 2697.36 3107.91 3175.05 3355.8  2983.18 2935.   3191.81'
    x_positions = '65392.28 87521.94 74200.75 33056.45 38578.91 73717.72 79407.04 61057.73 65502.19 34105.4
    15495.84 39069.52 47270.14 81182.91 19481.53 53804.82 16891.68 80337.59 52473.89 51584.05'
    y_positions = '57040.12 66912.04 21248.38 73556.73 87021.67 91228.2  35742.18 89239.61 9333.9
    11173.35 81164.89 38512.26 52276.59 55357.99 23911.9  24603.55 55058.71 80676.15 65521.33 78256.87'
    z_positions = '0.       0.       0.       0.       0.       0.       0.       0. 0.       0.
    0.       0.       0.       0.       0.       0. 0.       0.       0.       0.'
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
