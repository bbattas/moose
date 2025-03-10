##############################################################################
# File: testbed1.i
# File Location: /examples/agg/01_initial_testing/01_misc
# Created Date: Wednesday March 5th 2025
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Wednesday March 5th 2025
# Modified By: Brandon Battas
# -----
# Description:
#  First misc input for testing all kinds of things and getting a feel for
#   the GG problem input
#
#
##############################################################################

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  nz = 0
  xmax = 1000
  ymax = 1000
  zmax = 0
  elem_type = QUAD4
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 4
  var_name_base = 'gr'
[]

[Variables]
  [PolycrystalVariables]
  []
[]

[UserObjects]
  [voronoi]
    type = PolycrystalVoronoi
    rand_seed = 102
    grain_num = 4
    coloring_algorithm = bt
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
  []
[]

[AuxVariables]
  [bnds]
    order = FIRST
    family = LAGRANGE
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
  # Halos
  [halo0]
    order = CONSTANT
    family = MONOMIAL
  []
  [halo1]
    order = CONSTANT
    family = MONOMIAL
  []
  [halo2]
    order = CONSTANT
    family = MONOMIAL
  []
  [halo3]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [PolycrystalKernel]
  []
[]

[AuxKernels]
  [BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'timestep_end'
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = timestep_begin
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = timestep_begin
  []
  # HALOS
  [halo0]
    type = FeatureFloodCountAux
    variable = halo0
    map_index = 0
    field_display = HALOS
    flood_counter = grain_tracker
  []
  [halo1]
    type = FeatureFloodCountAux
    variable = halo1
    map_index = 1
    field_display = HALOS
    flood_counter = grain_tracker
  []
  [halo2]
    type = FeatureFloodCountAux
    variable = halo2
    map_index = 2
    field_display = HALOS
    flood_counter = grain_tracker
  []
  [halo3]
    type = FeatureFloodCountAux
    variable = halo3
    map_index = 3
    field_display = HALOS
    flood_counter = grain_tracker
  []
[]

[BCs]
  [Periodic]
    [All]
      auto_direction = 'x y'
    []
  []
[]

[Materials]
  [Moly_GB]
    type = GBEvolution
    time_scale = 1.0
    GBmob0 = 3.986e-6
    T = 500 # K
    wGB = 60 # nm
    Q = 1.0307
    GBenergy = 2.4
  []
  [incl_test01]
    type = GGInclinationMaterial
    inclination_name = inclination_mat01
    i_value = 0
    j_value = 1
  []
  [incl_test02]
    type = GGInclinationMaterial
    inclination_name = inclination_mat02
    i_value = 0
    j_value = 2
  []
  [incl_test12]
    type = GGInclinationMaterial
    inclination_name = inclination_mat12
    i_value = 1
    j_value = 2
    outputs = 'exodus'
  []
  [incl_01]
    type = ParsedMaterial
    property_name = incl_01
    material_property_names = 'inclination_mat01'
    expression = 'inclination_mat01'
    outputs = 'exodus'
  []
  [incl_02]
    type = ParsedMaterial
    property_name = incl_02
    material_property_names = 'inclination_mat02'
    expression = 'inclination_mat02'
    outputs = 'exodus'
  []
  [incl_12]
    type = ParsedMaterial
    property_name = incl_12
    material_property_names = 'inclination_mat12'
    expression = 'inclination_mat12'
    outputs = 'exodus'
  []
[]

[UserObjects]
  [grain_tracker]
    type = GrainTracker
    # variable = eta
    threshold = 0.2
    connecting_threshold = 0.08
    compute_halo_maps = true
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
  []
[]

[Postprocessors]
  [gr1area]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
    execute_on = 'initial timestep_end'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  nl_max_its = 30
  l_max_its = 60
  l_tol = 1e-06 #4
  nl_rel_tol = 1e-8 #default is 1e-8
  # nl_abs_tol = 1e-14 #only needed when near equilibrium or veeeery small dt

  start_time = 0.0
  num_steps = 2
  dt = 1
[]

[Outputs]
  exodus = true
[]
