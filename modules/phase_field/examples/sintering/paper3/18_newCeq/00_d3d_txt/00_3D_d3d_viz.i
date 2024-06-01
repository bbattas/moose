##############################################################################
# File: 00_3D_d3d_viz.i
# File Location: /examples/sintering/paper3/18_newCeq/00_d3d_txt
# Created Date: Friday May 31st 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Friday May 31st 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Bare minimum input for visualization of the d3d input
#  Just this simple input has 1,586,610 DoFs
#  Too large for autoformat on save, have to skip it
#
##############################################################################

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = 3D_20x20x20um_8umavg_20pore_shortExternal.txt
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = ebsd_mesh
    combinatorial_geometry = 'x > 19000'
    block_id = 1
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 9
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
    x1 = 25000
    x2 = 50000
    y1 = -5000
    y2 = 30000
  []
  [voidIC2] #Internal
    type = SpecifiedSmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0.01
    radii = '1804.12 1567.47 1902.78 1929.8  1576.57 1692.6  1963.42 1794.93 1925.8
    1875.61 1924.63 1571.62 1716.74 1737.71 1609.51 1807.06 1975.33 1807.03
    1749.88 1660.38'
    x_positions = '15347.09 14654.89 10212.18  9023.26  8108.55  9288.76  5540.02 13133.5
    5396.92 15704.63  4400.98  8903.65  3877.76 12881.87  4046.32 16040.23
   16040.4   4075.51 15879.02 10002.46'
    y_positions = '9239.32  7064.04 15460.85  8120.34  4337.69 11517.3   5497.48  5337.64
    15423.17 14274.87 13493.45 12680.29 10663.43 16308.02 11185.2  12067.31
    15457.17  3742.25  3638.    7889.2'
    z_positions = '4257.46 14931.82  5279.29  3869.06 15257.39 16156.82  9271.6   8473.5
    12734.09  9881.96  3762.12  9832.09 15919.8  14262.3   9451.76 14931.1
     3773.7   3622.63  3552.39 12105.96'
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
