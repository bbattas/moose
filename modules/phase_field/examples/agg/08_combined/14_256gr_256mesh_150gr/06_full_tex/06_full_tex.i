##############################################################################
# File: 06_full_tex.i
# File Location: /examples/agg/08_combined/14_256gr_256mesh_150gr/06_full_tex
# Created Date: Tuesday July 28th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Tuesday July 28th 2026
# Modified By: Brandon Battas
# -----
# Description:
#
#
#
#
##############################################################################

randseed = 20

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 384
    ny = 384
    xmin = 0
    xmax = 256
    ymin = 0
    ymax = 256
  []
  parallel_type = DISTRIBUTED
  second_order = false
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 25
  var_name_base = gr
[]

[Variables]
  [PolycrystalVariables]
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
  []
[]

[UserObjects]
  [euler_file]
    type = EulerAngleTxtFileReader
    file_name = '../../00_sub/euler_256_textured.txt'
  []
  [voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = bt #jp #bt
    grain_num = 256
    rand_seed = ${randseed}
    int_width = 3 #6
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
    expression = 'grain_tracker < 50'
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
  [PolycrystalInclinationKernel]
    variable_mobility = false
    hgb_mask = hgb
  []
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
    gb_mode = FULL
    bulk_scalar = 0.75
    alpha_tol = 8 #10
    hgbalpha_tol = 5
    # ifunc_a = 0.1
    kappa_name = kappa_op
    gbe_iso_name = gbe_iso
    enable_misorientation = true
    euler_angle_provider = euler_file
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
  nl_rel_tol = 1e-6#1e-8 #1e-10 # Relative tolerance for nonlienar solves
  # nl_abs_tol = 1e-10

  start_time = 0.0
  # num_steps = 3
  end_time = 6.24
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
  line_search = 'none'
[]

[Outputs]
  # exodus = true # Exodus file will be outputted
  #nemesis = true
  console = true
  # csv = true
  [exodus]
    type = Exodus
    sync_times = '0 0.052 0.104 0.156 0.208 0.26 0.312 0.364 0.416 0.468 0.52 0.572 0.624 0.676 0.728 0.78
    0.832 0.884 0.936 0.988 1.04 1.092 1.144 1.196 1.248 1.3 1.352 1.404 1.456 1.508 1.56 1.612
    1.664 1.716 1.768 1.82 1.872 1.924 1.976 2.028 2.08 2.132 2.184 2.236 2.288 2.34 2.392 2.444
    2.496 2.548 2.6 2.652 2.704 2.756 2.808 2.86 2.912 2.964 3.016 3.068 3.12 3.172 3.224 3.276
    3.328 3.38 3.432 3.484 3.536 3.588 3.64 3.692 3.744 3.796 3.848 3.9 3.952 4.004 4.056 4.108
    4.16 4.212 4.264 4.316 4.368 4.42 4.472 4.524 4.576 4.628 4.68 4.732 4.784 4.836 4.888 4.94
    4.992 5.044 5.096 5.148 5.2 5.252 5.304 5.356 5.408 5.46 5.512 5.564 5.616 5.668 5.72 5.772
    5.824 5.876 5.928 5.98 6.032 6.084 6.136 6.188 6.24'
  []
  [csv]
    type = CSV
    sync_times = '0 0.052 0.104 0.156 0.208 0.26 0.312 0.364 0.416 0.468 0.52 0.572 0.624 0.676 0.728 0.78
    0.832 0.884 0.936 0.988 1.04 1.092 1.144 1.196 1.248 1.3 1.352 1.404 1.456 1.508 1.56 1.612
    1.664 1.716 1.768 1.82 1.872 1.924 1.976 2.028 2.08 2.132 2.184 2.236 2.288 2.34 2.392 2.444
    2.496 2.548 2.6 2.652 2.704 2.756 2.808 2.86 2.912 2.964 3.016 3.068 3.12 3.172 3.224 3.276
    3.328 3.38 3.432 3.484 3.536 3.588 3.64 3.692 3.744 3.796 3.848 3.9 3.952 4.004 4.056 4.108
    4.16 4.212 4.264 4.316 4.368 4.42 4.472 4.524 4.576 4.628 4.68 4.732 4.784 4.836 4.888 4.94
    4.992 5.044 5.096 5.148 5.2 5.252 5.304 5.356 5.408 5.46 5.512 5.564 5.616 5.668 5.72 5.772
    5.824 5.876 5.928 5.98 6.032 6.084 6.136 6.188 6.24'
  []
  #
  print_linear_residuals = true
  print_nonlinear_residuals = true
  #
  checkpoint = false
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 3
  # []
  # [pgraph]
  #   type = PerfGraphOutput
  #   execute_on = 'initial final' # Default is "final"
  #   level = 2 # Default is 1
  #   heaviest_branch = true # Default is false
  #   heaviest_sections = 7 # Default is 0
  # []
  # file_base = 01_base_a${a_tol}_i${i_tol}_out
[]
