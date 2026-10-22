##############################################################################
# File: 02_halo12.i
# File Location: /examples/agg/09_miso_testing/03_ic_duplicate_grain_testing/02_halo12
# Created Date: Thursday October 22nd 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Thursday October 22nd 2026
# Modified By: Brandon Battas
# -----
# Description:
#  testing a nosolve/barely solve 1 step just to check loading from IC
#  9 (gr20), 12 (gr22), 15 (gr25)
#
#
##############################################################################

# randseed = 20
[Problem]
  solve = false
[]

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 768
    ny = 768
    xmin = 0
    xmax = 512
    ymin = 0
    ymax = 512
  []
  parallel_type = DISTRIBUTED
  second_order = false
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 26
  var_name_base = gr
[]

[Variables]
  [PolycrystalVariables]
  []
[]


[ICs]
  [ic_gr0]
    type = SolutionIC
    from_variable = gr0
    variable = gr0
    solution_uo = solution_userobject
  []
  [ic_gr1]
    type = SolutionIC
    from_variable = gr1
    variable = gr1
    solution_uo = solution_userobject
  []
  [ic_gr2]
    type = SolutionIC
    from_variable = gr2
    variable = gr2
    solution_uo = solution_userobject
  []
  [ic_gr3]
    type = SolutionIC
    from_variable = gr3
    variable = gr3
    solution_uo = solution_userobject
  []
  [ic_gr4]
    type = SolutionIC
    from_variable = gr4
    variable = gr4
    solution_uo = solution_userobject
  []
  [ic_gr5]
    type = SolutionIC
    from_variable = gr5
    variable = gr5
    solution_uo = solution_userobject
  []
  [ic_gr6]
    type = SolutionIC
    from_variable = gr6
    variable = gr6
    solution_uo = solution_userobject
  []
  [ic_gr7]
    type = SolutionIC
    from_variable = gr7
    variable = gr7
    solution_uo = solution_userobject
  []
  [ic_gr8]
    type = SolutionIC
    from_variable = gr8
    variable = gr8
    solution_uo = solution_userobject
  []
  [ic_gr9]
    type = SolutionIC
    from_variable = gr9
    variable = gr9
    solution_uo = solution_userobject
  []
  [ic_gr10]
    type = SolutionIC
    from_variable = gr10
    variable = gr10
    solution_uo = solution_userobject
  []
  [ic_gr11]
    type = SolutionIC
    from_variable = gr11
    variable = gr11
    solution_uo = solution_userobject
  []
  [ic_gr12]
    type = SolutionIC
    from_variable = gr12
    variable = gr12
    solution_uo = solution_userobject
  []
  [ic_gr13]
    type = SolutionIC
    from_variable = gr13
    variable = gr13
    solution_uo = solution_userobject
  []
  [ic_gr14]
    type = SolutionIC
    from_variable = gr14
    variable = gr14
    solution_uo = solution_userobject
  []
  [ic_gr15]
    type = SolutionIC
    from_variable = gr15
    variable = gr15
    solution_uo = solution_userobject
  []
  [ic_gr16]
    type = SolutionIC
    from_variable = gr16
    variable = gr16
    solution_uo = solution_userobject
  []
  [ic_gr17]
    type = SolutionIC
    from_variable = gr17
    variable = gr17
    solution_uo = solution_userobject
  []
  [ic_gr18]
    type = SolutionIC
    from_variable = gr18
    variable = gr18
    solution_uo = solution_userobject
  []
  [ic_gr19]
    type = SolutionIC
    from_variable = gr19
    variable = gr19
    solution_uo = solution_userobject
  []
  [ic_gr20]
    type = SolutionIC
    from_variable = gr20
    variable = gr20
    solution_uo = solution_userobject
  []
  [ic_gr21]
    type = SolutionIC
    from_variable = gr21
    variable = gr21
    solution_uo = solution_userobject
  []
  [ic_gr22]
    type = SolutionIC
    from_variable = gr22
    variable = gr22
    solution_uo = solution_userobject
  []
  [ic_gr23]
    type = ConstantIC
    variable = gr23
    value = 0.0
  []
  [ic_gr24]
    type = ConstantIC
    variable = gr24
    value = 0.0
  []
  [ic_gr25]
    type = ConstantIC
    variable = gr25
    value = 0.0
  []
[]

[UserObjects]
  # [euler_file]
  #   type = EulerAngleTxtFileReader
  #   file_name = '../../00_sub/euler_1k_mackenzie.txt'
  # []
  # [voronoi]
  #   type = PolycrystalVoronoi
  #   coloring_algorithm = bt #jp #bt
  #   grain_num = 256
  #   rand_seed = ${randseed}
  #   int_width = 3 #6
  # []
  [solution_userobject]
    type = SolutionUserObject
    mesh = '../00_sub/16_iso_relax_2400_1000_12halo.e'
    timestep = 'LATEST'
    system_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13 gr14 gr15 gr16 gr17 gr18 gr19 gr20 gr21 gr22'
  []
  [grain_tracker]
    type = GrainTracker
    # variable = 'gr0 gr1 gr2 gr3 gr4'
    threshold = 0.01
    connecting_threshold = 0.01
    compute_var_to_feature_map = true
    compute_halo_maps = false
    execute_on = 'initial timestep_end'
    halo_level = 6
    verbosity_level = 2
    # use_less_than_threshold_comparison = true
  []
  [term]
    type = Terminator
    expression = 'grain_tracker < 100'
  []
[]

[AuxVariables]
  [bnds]
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  # [halos]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
[]

[Kernels]
  [PolycrystalKernel]
    variable_mobility = false
    # order = SECOND
  []
  # [PolycrystalInclinationKernel]
  #   variable_mobility = false
  #   hgb_mask = hgb
  # []
[]

[AuxKernels]
  [bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  []
  # [halos]
  #   type = FeatureFloodCountAux
  #   variable = halos
  #   flood_counter = grain_tracker
  #   field_display = HALOS
  #   execute_on = 'initial timestep_end'
  # []
[]

[BCs]
  [Periodic]
    [All]
      auto_direction = 'x y'
    []
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L0             kappa_op   mu      gbe_iso   '
    prop_values = '2.7815e-6    2.100e07    1.305e7    1.125e7   '
  []
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    hgb_threshold = 0
    func_type = COMBINED
    # output_properties = 'hgb'
    # outputs = 'exodus'
  []
  [new_aniso_mat]
    type = GBCombinedAnisotropyMaterial
    gb_id_method = graintracker
    grain_tracker = grain_tracker
    gb_mode = ISO
    iso_gbe = 0.95
    bulk_scalar = 0.75
    alpha_tol = 10
    hgbalpha_tol = 5
    # ifunc_a = 0.1
    kappa_name = kappa_op
    gbe_iso_name = gbe_iso
    enable_misorientation = false
    # euler_angle_provider = euler_file
    w_inc = 1
    w_miso = 1
    output_properties = 'gt_num L gamma_asymm int_width gbe_norm gbe_gb gb_notj_mask gb_tj_mask'
    outputs = 'exodus'
  []
  [gbe_notj_masked]
    type = ParsedMaterial
    property_name = gbe_notj_masked
    material_property_names = 'gbe_gb gb_notj_mask'
    expression = 'gbe_gb * gb_notj_mask'
    outputs = 'exodus'
  []
  [gbe_tj_masked]
    type = ParsedMaterial
    property_name = gbe_tj_masked
    material_property_names = 'gbe_gb gb_tj_mask'
    expression = 'gbe_gb * gb_tj_mask'
    outputs = 'exodus'
  []
[]

[Postprocessors]
  [runtime]
    type = PerfGraphData
    section_name = "Root"
    data_type = TOTAL
  []
  [nl_its]
    type = NumNonlinearIterations
  []
  [l_its]
    type = NumLinearIterations
  []
  [timestep]
    type = TimestepSize
  []
  [n_elements]
    type = NumElements
    execute_on = 'initial timestep_end'
  []
  [n_nodes]
    type = NumNodes
    execute_on = 'initial timestep_end'
  []
  [DOFs]
    type = NumDOFs
    execute_on = 'initial timestep_end'
  []
  # IW and Gamma_asymm
  [avg_iw]
    type = ElementAverageMaterialProperty
    mat_prop = int_width
  []
  [avg_gamma]
    type = ElementAverageMaterialProperty
    mat_prop = gamma_asymm
  []
  [gamma_min]
    type = ElementExtremeMaterialProperty
    mat_prop = gamma_asymm
    value_type = MIN
  []
  [gamma_max]
    type = ElementExtremeMaterialProperty
    mat_prop = gamma_asymm
    value_type = MAX
  []
  [iw_min]
    type = ElementExtremeMaterialProperty
    mat_prop = int_width
    value_type = MIN
  []
  [iw_max]
    type = ElementExtremeMaterialProperty
    mat_prop = int_width
    value_type = MAX
  []
    # GBE
  [gbe_notj_total]
    type = ElementIntegralMaterialProperty
    mat_prop = gbe_notj_masked
  []
  [gbe_tj_total]
    type = ElementIntegralMaterialProperty
    mat_prop = gbe_tj_masked
  []
  [gb_notj_area]
    type = ElementIntegralMaterialProperty
    mat_prop = gb_notj_mask
  []
  [gb_tj_area]
    type = ElementIntegralMaterialProperty
    mat_prop = gb_tj_mask
  []
  [gbe_notj_avg]
    type = ParsedPostprocessor
    pp_names = 'gbe_notj_total gb_notj_area'
    expression = 'gbe_notj_total / gb_notj_area'
  []
  [gbe_tj_avg]
    type = ParsedPostprocessor
    pp_names = 'gbe_tj_total gb_tj_area'
    expression = 'gbe_tj_total / gb_tj_area'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre boomeramg 101 ds'
  # petsc_options_iname = '-pc_type -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  # petsc_options_value = 'asm        preonly       lu           2'

  l_max_its = 60 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 12 # Max number of nonlinear iterations
  nl_rel_tol = 1e-8 #1e-10 # Relative tolerance for nonlienar solves
  # nl_abs_tol = 1e-10

  start_time = 0.0
  num_steps = 1
  # end_time = 3
  # dt = 0.1
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01 #0.001
    # cutback_factor = 0.9
    # growth_factor = 1.1
    optimal_iterations = 6
    linear_iteration_ratio = 1e5
  []

  automatic_scaling = true
  compute_scaling_once = false
  # line_search = 'none'
[]

[Outputs]
  exodus = true # Exodus file will be outputted
  #nemesis = true
  console = true
  csv = true
  #
  print_linear_residuals = true
  print_nonlinear_residuals = true
  #
  # checkpoint = false
  [checkpoint]
    type = Checkpoint
    num_files = 3
  []
  # [pgraph]
  #   type = PerfGraphOutput
  #   execute_on = 'initial final' # Default is "final"
  #   level = 2 # Default is 1
  #   heaviest_branch = true # Default is false
  #   heaviest_sections = 7 # Default is 0
  # []
  # file_base = 01_base_a${a_tol}_i${i_tol}_out
[]
