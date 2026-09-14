##############################################################################
# File: 02_1k_iso_tex.i
# File Location: /examples/agg/09_miso_testing/02_1k_isomatch_relaxed/02_1k_iso_tex
# Created Date: Monday August 10th 2026
# Author: Battas,Brandon Scott (bbattas@ufl.edu)
# -----
# Last Modified: Monday September 14th 2026
# Modified By: Battas,Brandon Scott
# -----
# Description:
#  iso with matching gbe avg (Step 1 no change from relaxed is 0.83(89))
#  Might want to wait for miso to run then check a few timesteps in and readjust
#
#
##############################################################################

# randseed = 20

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
  op_num = 24
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
    type = ConstantIC
    variable = gr21
    value = 0.0
  []
  [ic_gr22]
    type = ConstantIC
    variable = gr22
    value = 0.0
  []
  [ic_gr23]
    type = ConstantIC
    variable = gr23
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
    mesh = '../../00_sub/relaxed_2400_1000.e'
    timestep = 'LATEST'
    system_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13 gr14 gr15 gr16 gr17 gr18 gr19 gr20'
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
    iso_gbe = 0.84
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
  # num_steps = 3
  end_time = 3
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
  # exodus = true # Exodus file will be outputted
  #nemesis = true
  console = true
  # csv = true
  [exodus]
    type = Exodus
    sync_only = true
    sync_times = '0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15
                  0.16 0.17 0.18 0.19 0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3 0.31
                  0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47
                  0.48 0.49 0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6 0.61 0.62 0.63
                  0.64 0.65 0.66 0.67 0.68 0.69 0.7 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
                  0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95
                  0.96 0.97 0.98 0.99 1 1.01 1.02 1.03 1.04 1.05 1.06 1.07 1.08 1.09 1.1 1.11
                  1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.2 1.21 1.22 1.23 1.24 1.25 1.26 1.27
                  1.28 1.29 1.3 1.31 1.32 1.33 1.34 1.35 1.36 1.37 1.38 1.39 1.4 1.41 1.42 1.43
                  1.44 1.45 1.46 1.47 1.48 1.49 1.5 1.51 1.52 1.53 1.54 1.55 1.56 1.57 1.58 1.59
                  1.6 1.61 1.62 1.63 1.64 1.65 1.66 1.67 1.68 1.69 1.7 1.71 1.72 1.73 1.74 1.75
                  1.76 1.77 1.78 1.79 1.8 1.81 1.82 1.83 1.84 1.85 1.86 1.87 1.88 1.89 1.9 1.91
                  1.92 1.93 1.94 1.95 1.96 1.97 1.98 1.99 2 2.01 2.02 2.03 2.04 2.05 2.06 2.07
                  2.08 2.09 2.1 2.11 2.12 2.13 2.14 2.15 2.16 2.17 2.18 2.19 2.2 2.21 2.22 2.23
                  2.24 2.25 2.26 2.27 2.28 2.29 2.3 2.31 2.32 2.33 2.34 2.35 2.36 2.37 2.38 2.39
                  2.4 2.41 2.42 2.43 2.44 2.45 2.46 2.47 2.48 2.49 2.5 2.51 2.52 2.53 2.54 2.55
                  2.56 2.57 2.58 2.59 2.6 2.61 2.62 2.63 2.64 2.65 2.66 2.67 2.68 2.69 2.7 2.71
                  2.72 2.73 2.74 2.75 2.76 2.77 2.78 2.79 2.8 2.81 2.82 2.83 2.84 2.85 2.86 2.87
                  2.88 2.89 2.9 2.91 2.92 2.93 2.94 2.95 2.96 2.97 2.98 2.99 3'
  []
  [csv]
    type = CSV
    sync_only = true
    sync_times = '0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15
                  0.16 0.17 0.18 0.19 0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29 0.3 0.31
                  0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47
                  0.48 0.49 0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.6 0.61 0.62 0.63
                  0.64 0.65 0.66 0.67 0.68 0.69 0.7 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79
                  0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89 0.9 0.91 0.92 0.93 0.94 0.95
                  0.96 0.97 0.98 0.99 1 1.01 1.02 1.03 1.04 1.05 1.06 1.07 1.08 1.09 1.1 1.11
                  1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 1.2 1.21 1.22 1.23 1.24 1.25 1.26 1.27
                  1.28 1.29 1.3 1.31 1.32 1.33 1.34 1.35 1.36 1.37 1.38 1.39 1.4 1.41 1.42 1.43
                  1.44 1.45 1.46 1.47 1.48 1.49 1.5 1.51 1.52 1.53 1.54 1.55 1.56 1.57 1.58 1.59
                  1.6 1.61 1.62 1.63 1.64 1.65 1.66 1.67 1.68 1.69 1.7 1.71 1.72 1.73 1.74 1.75
                  1.76 1.77 1.78 1.79 1.8 1.81 1.82 1.83 1.84 1.85 1.86 1.87 1.88 1.89 1.9 1.91
                  1.92 1.93 1.94 1.95 1.96 1.97 1.98 1.99 2 2.01 2.02 2.03 2.04 2.05 2.06 2.07
                  2.08 2.09 2.1 2.11 2.12 2.13 2.14 2.15 2.16 2.17 2.18 2.19 2.2 2.21 2.22 2.23
                  2.24 2.25 2.26 2.27 2.28 2.29 2.3 2.31 2.32 2.33 2.34 2.35 2.36 2.37 2.38 2.39
                  2.4 2.41 2.42 2.43 2.44 2.45 2.46 2.47 2.48 2.49 2.5 2.51 2.52 2.53 2.54 2.55
                  2.56 2.57 2.58 2.59 2.6 2.61 2.62 2.63 2.64 2.65 2.66 2.67 2.68 2.69 2.7 2.71
                  2.72 2.73 2.74 2.75 2.76 2.77 2.78 2.79 2.8 2.81 2.82 2.83 2.84 2.85 2.86 2.87
                  2.88 2.89 2.9 2.91 2.92 2.93 2.94 2.95 2.96 2.97 2.98 2.99 3'
  []
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
