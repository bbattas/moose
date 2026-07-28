##############################################################################
# File: 04_miso_tex.i
# File Location: /examples/agg/08_combined/14_256gr_256mesh_150gr/04_miso_tex
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
    gb_mode = MISO
    bulk_scalar = 0.75
    alpha_tol = 10
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
  nl_rel_tol = 1e-8 #1e-10 # Relative tolerance for nonlienar solves
  # nl_abs_tol = 1e-10

  start_time = 0.0
  # num_steps = 3
  end_time = 3.6
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
    sync_times = '0 0.03 0.06 0.09 0.12 0.15 0.18 0.21 0.24 0.27 0.3 0.33 0.36 0.39 0.42 0.45
    0.48 0.51 0.54 0.57 0.6 0.63 0.66 0.69 0.72 0.75 0.78 0.81 0.84 0.87 0.9 0.93
    0.96 0.99 1.02 1.05 1.08 1.11 1.14 1.17 1.2 1.23 1.26 1.29 1.32 1.35 1.38 1.41
    1.44 1.47 1.5 1.53 1.56 1.59 1.62 1.65 1.68 1.71 1.74 1.77 1.8 1.83 1.86 1.89
    1.92 1.95 1.98 2.01 2.04 2.07 2.1 2.13 2.16 2.19 2.22 2.25 2.28 2.31 2.34 2.37
    2.4 2.43 2.46 2.49 2.52 2.55 2.58 2.61 2.64 2.67 2.7 2.73 2.76 2.79 2.82 2.85
    2.88 2.91 2.94 2.97 3 3.03 3.06 3.09 3.12 3.15 3.18 3.21 3.24 3.27 3.3 3.33
    3.36 3.39 3.42 3.45 3.48 3.51 3.54 3.57 3.6'
  []
  [csv]
    type = CSV
    sync_times = '0 0.03 0.06 0.09 0.12 0.15 0.18 0.21 0.24 0.27 0.3 0.33 0.36 0.39 0.42 0.45
    0.48 0.51 0.54 0.57 0.6 0.63 0.66 0.69 0.72 0.75 0.78 0.81 0.84 0.87 0.9 0.93
    0.96 0.99 1.02 1.05 1.08 1.11 1.14 1.17 1.2 1.23 1.26 1.29 1.32 1.35 1.38 1.41
    1.44 1.47 1.5 1.53 1.56 1.59 1.62 1.65 1.68 1.71 1.74 1.77 1.8 1.83 1.86 1.89
    1.92 1.95 1.98 2.01 2.04 2.07 2.1 2.13 2.16 2.19 2.22 2.25 2.28 2.31 2.34 2.37
    2.4 2.43 2.46 2.49 2.52 2.55 2.58 2.61 2.64 2.67 2.7 2.73 2.76 2.79 2.82 2.85
    2.88 2.91 2.94 2.97 3 3.03 3.06 3.09 3.12 3.15 3.18 3.21 3.24 3.27 3.3 3.33
    3.36 3.39 3.42 3.45 3.48 3.51 3.54 3.57 3.6'
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
