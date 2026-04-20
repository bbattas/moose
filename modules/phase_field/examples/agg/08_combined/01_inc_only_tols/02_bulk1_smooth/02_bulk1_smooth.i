##############################################################################
# File: 02_bulk1_smooth.i
# File Location: /examples/agg/08_combined/01_inc_only_tols/02_bulk1_smooth
# Created Date: Monday April 20th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Monday April 20th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  array test for bulk/gb scale of 1 and smoothed version of lins combo func
#  running 0/1/5/10/20/50/100/1000 for ai tols
#
#  5 cpus
##############################################################################

randseed = 20
a_tol = 0
i_tol = 0

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 80
    ny = 80
    xmin = 0
    xmax = 40
    ymin = 0
    ymax = 40
  []
  parallel_type = DISTRIBUTED
  second_order = false
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 8 #20
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
  [voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = bt #jp #bt
    grain_num = 8
    rand_seed = ${randseed}
    int_width = 4
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
    expression = 'grain_tracker < 4'
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
  [sum_inc]
    order = CONSTANT
    family = MONOMIAL
  []
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
  [sum_inc]
    type = MatVectorijSum
    property = inclination
    variable = sum_inc
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
  # [incl_vec]
  #   type = GGInclinationVector
  #   # grain_tracker = grain_tracker
  #   output_properties = 'inclination_vector ang_dist'
  #   gb_id_method = SWITCH
  #   # ffc = gr_flood_uo
  #   grain_tracker = grain_tracker
  #   hgb = hgb
  #   hgb_threshold = 0.5
  #   outputs = exodus
  # []
  # [constants]
  #   type = GenericConstantMaterial
  #   prop_names = 'L0         gamma_iso kappa_op  iw_iso gbe_iso'
  #   prop_values = '1.15382e-6   1.5    2.07337e7   6   4.60748e6'
  # []
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L0             kappa_op  int_width_iso   const_m      gbe_iso   '
    prop_values = '2.7823e-6    2.590909e7        6      1.282828e7  1.388384e7   '
  []
  # [constants_2]
  #   type = GenericConstantMaterial
  #   prop_names = 'L        '#gamma_asymm'# int_width'#   mu'
  #   prop_values = '1.15382e-6 '#' 1.5 '#'    6        '#5.521269e6'
  # []
  [GBInc_cos]
    type = GBInclination
    gb_id_method = GRAINTRACKER
    # ffc = gr_flood_uo
    grain_tracker = grain_tracker
    angular_func = ATAN_2D
    intol = ${i_tol}
    altol = ${a_tol}
    limit_umag = true
    bulk_scale = 1.0
    nongb_scale = 1.0
    # Inclination function
    inc_func = SMOOTH_LIN
    ifunc_a = 0 #${amag}
    ifunc_b = 2
    ifunc_c = 0
    ifunc_d = 0
    # L
    combine_gb_form = avg
    aniso_L = false
    noDeriv_L = true
    aniso_gbmob = false
    # Other Properties
    gb_energy_iso_name = gbe_iso
    kappa = kappa_op
    free_energy_m = 1.282828e7 #5.521269e6
    output_properties = 'gtnum int_width gamma_asymm L'
    outputs = 'exodus'
  []
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    hgb_threshold = 0
    func_type = COMBINED
    # output_properties = 'hgb'
    # outputs = 'exodus'
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
  [avg_grain_volumes]
    type = AverageGrainVolume
    feature_counter = grain_tracker
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
[]

# [VectorPostprocessors]
#   [grain_volumes]
#     type = FeatureVolumeVectorPostprocessor
#     flood_counter = grain_tracker
#     execute_on = 'initial timestep_end'
#     # output_centroids = true
#   []
#   # [./mem]
#   #   type = VectorMemoryUsage
#   #   execute_on = 'INITIAL TIMESTEP_END NONLINEAR LINEAR'
#   #   report_peak_value = true
#   #   #mem_units = kilobytes # or bytes, megabytes, gigabytes
#   #   mem_units = gigabytes
#   # [../]
# []

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
  # end_time = 30
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
  file_base = 02_bulk1_smooth_a${a_tol}_i${i_tol}_out
[]
