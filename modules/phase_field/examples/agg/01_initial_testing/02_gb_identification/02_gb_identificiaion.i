##############################################################################
# File: 02_gb_identificiaion.i
# File Location: /examples/agg/01_initial_testing/02_gb_identification
# Created Date: Friday March 7th 2025
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Monday March 10th 2025
# Modified By: Brandon Battas
# -----
# Description:
#  Testing with ebsd mesh for more grains than ops and trying to recreate
#   the GB identification for ij from ComputeGBMisorientationType
#
#
##############################################################################

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = '2D_500x500_ggTest3.txt'
  []
  # uniform_refine = 1
  parallel_type = DISTRIBUTED
[]

[GlobalParams]
  op_num = 8
  var_name_base = 'gr'
[]

[Variables]
  [PolycrystalVariables]
  []
[]

[UserObjects]
  [ebsd_reader]
    type = EBSDReader
  []
  [ebsd]
    type = PolycrystalEBSD
    coloring_algorithm = bt
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    # output_adjacency_matrix = true
  []
  [grain_tracker]
    type = GrainTracker
    # variable = eta
    threshold = 0.001
    connecting_threshold = 0.008
    compute_halo_maps = true
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
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
  [ghost_regions]
    order = CONSTANT
    family = MONOMIAL
  []
  [ebsd_ic]
    order = CONSTANT
    family = MONOMIAL
  []
  [ebsd_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  # Halos
  [halos]
    order = CONSTANT
    family = MONOMIAL
  []
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
  [ghosted_entities]
    type = FeatureFloodCountAux
    variable = ghost_regions
    flood_counter = grain_tracker
    field_display = GHOSTED_ENTITIES
    execute_on = 'initial timestep_end'
  []
  [ebsd_ic_aux]
    type = EBSDReaderPointDataAux
    variable = ebsd_ic
    ebsd_reader = ebsd_reader
    data_name = 'feature_id'
    execute_on = 'initial timestep_end'
  []
  [ebsd_grains]
    type = EBSDReaderAvgDataAux
    data_name = 'feature_id'
    ebsd_reader = ebsd_reader
    grain_tracker = grain_tracker
    variable = ebsd_grains
    execute_on = 'initial timestep_end'
  []
  [halos]
    type = FeatureFloodCountAux
    variable = halos
    flood_counter = grain_tracker
    field_display = HALOS
    execute_on = 'initial timestep_end'
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

# [BCs]
#   [Periodic]
#     [All]
#       auto_direction = 'x y'
#     []
#   []
# []

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
    ebsd_reader = ebsd_reader
    grain_tracker = grain_tracker
    output_properties = 'gb_id grains_on_gb'
    outputs = 'exodus'
  []
  # [incl_test02]
  #   type = GGInclinationMaterial
  #   inclination_name = inclination_mat02
  #   i_value = 0
  #   j_value = 2
  # []
  # [incl_test12]
  #   type = GGInclinationMaterial
  #   inclination_name = inclination_mat12
  #   i_value = 1
  #   j_value = 2
  #   outputs = 'exodus'
  # []
  [incl_01]
    type = ParsedMaterial
    property_name = incl_01
    material_property_names = 'inclination_mat01'
    expression = 'inclination_mat01'
    outputs = 'exodus'
  []
  # [incl_02]
  #   type = ParsedMaterial
  #   property_name = incl_02
  #   material_property_names = 'inclination_mat02'
  #   expression = 'inclination_mat02'
  #   outputs = 'exodus'
  # []
  # [incl_12]
  #   type = ParsedMaterial
  #   property_name = incl_12
  #   material_property_names = 'inclination_mat12'
  #   expression = 'inclination_mat12'
  #   outputs = 'exodus'
  # []
  [sumgr]
    type = ParsedMaterial
    property_name = sumgr
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    expression = 'gr0 + gr1 + gr2 + gr3 + gr4 + gr5 + gr6 + gr7'
    outputs = 'exodus'
  []
  [Lij_test]
    type = ParsedMaterial
    property_name = Lij_test
    constant_names = 'L01 L02 L03 L04 L05 L06 L07'
    constant_expressions = '1 2 3 4 5 6 7'
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    expression = '(L01*(gr0^2*gr1^2) + L02*(gr0^2*gr2^2) + L03*(gr0^2*gr3^2)
     + L04*(gr0^2*gr4^2) + L05*(gr0^2*gr5^2) + L06*(gr0^2*gr6^2) + L07*(gr0^2*gr7^2)) / (
       gr0^2*gr1^2 + gr0^2*gr2^2 + gr0^2*gr3^2 + gr0^2*gr4^2 + gr0^2*gr5^2 + gr0^2*gr6^2 + gr0^2*gr7^2)'
    outputs = 'exodus'
  []
[]

# [Postprocessors]
#   [gr1area]
#     type = ElementIntegralVariablePostprocessor
#     variable = gr1
#     execute_on = 'initial timestep_end'
#   []
# []

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
  num_steps = 5
  dt = 0.1
[]

[Outputs]
  exodus = true
[]
