##############################################################################
# File: 04_miso_only_newparams.i
# File Location: /examples/agg/07_DOE_examples/04_miso_only_newparams
# Created Date: Tuesday March 31st 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Tuesday March 31st 2026
# Modified By: Brandon Battas
# -----
# Description:
#  refined mesh version of miso only using the new params for iw of like
#   2.7-8.6 ish maximum range
#
#
##############################################################################

# file = 'PFTrainingDataset_ISOTEST'
randseed = 10
# amag = 0.05

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 256
    ny = 256
    xmin = 0
    xmax = 256
    ymin = 0
    ymax = 256
  []
  parallel_type = DISTRIBUTED
  second_order = false
  uniform_refine = 1
[]

[GlobalParams]
  op_num = 30 #20
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
    file_name = '../00_sub/euler_487.txt'
  []
  [voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = bt #jp #bt
    grain_num = 256
    rand_seed = ${randseed}
    int_width = 6
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
  # [sum_inc]
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
  # [sum_inc]
  #   type = MatVectorijSum
  #   property = inclination
  #   variable = sum_inc
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
  # [constants]
  #   type = GenericConstantMaterial
  #   prop_names = 'L0             kappa_op  int_width_iso   mu      sigma' #'gamma_asymm'
  #   prop_values = '1.15382e-6  2.07337e7        6      5.521269e6  4.60748e6' #'1.5    '
  # []
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L0             kappa_op  int_width_iso   mu      sigma   '
    prop_values = '2.7823e-6    2.590909e7        6      1.282828e7  1.388384e7   '
  []
  # [constants_2]
  #   type = GenericConstantMaterial
  #   prop_names = 'L        '#gamma_asymm'# int_width'#   mu'
  #   prop_values = '1.15382e-6 '#' 1.5 '#'    6        '#5.521269e6'
  # []
  # [GBInc_cos]
  #   type = GBInclination
  #   gb_id_method = GRAINTRACKER
  #   # ffc = gr_flood_uo
  #   grain_tracker = grain_tracker
  #   angular_func = ATAN_2D
  #   intol = 0
  #   altol = 5
  #   limit_umag = true
  #   # Inclination function
  #   inc_func = COS
  #   ifunc_a = ${amag}
  #   ifunc_b = 2
  #   ifunc_c = 0
  #   ifunc_d = 0
  #   # L
  #   combine_gb_form = avg
  #   aniso_L = false
  #   noDeriv_L = true
  #   aniso_gbmob = false
  #   # Other Properties
  #   gb_energy_iso_name = gbe_iso
  #   kappa = kappa_op
  #   free_energy_m = 5.521269e6
  #   # output_properties = 'gtnum int_width gamma_asymm L'
  #   # outputs = 'exodus'
  # []
  [GBM]
    type = GBMisorientationTxt
    euler_angle_provider = euler_file
    grain_tracker = grain_tracker
    kappa = kappa_op
    output_properties = 'misorientation twist_energy tilt_energy f_miso gamma_asymm int_width'
    outputs = exodus
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
  [dt]
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
  print_linear_residuals = false
  print_nonlinear_residuals = false
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
  # file_base = 02_inc_only_a${amag}_cos/02_inc_only_a${amag}_cos_out
  # [csv]
  #   type = CSV
  #   execute_on = 'INITIAL TIMESTEP_END'
  #   # execute_on = 'CUSTOM'
  #   sync_only = true
  #   sync_times = '0.0 0.3 0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5
  #   4.8 5.1 5.4 5.7 6.0 6.3 6.6 6.9 7.2 7.5 7.8 8.1 8.4 8.7 9.0 9.3 9.6
  #   9.9 10.2 10.5 10.8 11.1 11.4 11.7 12.0 12.3 12.6 12.9 13.2 13.5 13.8 14.1
  #   14.4 14.7 15.0 15.3 15.6 15.9 16.2 16.5 16.8 17.1 17.4 17.7 18.0 18.3 18.6 18.9
  #   19.2 19.5 19.8 20.1 20.4 20.7 21.0 21.3 21.6 21.9 22.2 22.5 22.8 23.1 23.4 23.7
  #   24.0 24.3 24.6 24.9 25.2 25.5 25.8 26.1 26.4 26.7 27.0 27.3 27.6 27.9 28.2 28.5
  #   28.8 29.1 29.4 29.7 30.0'
  # []
  # [exodus]
  #   type = Exodus
  #   execute_on = 'INITIAL TIMESTEP_END'
  #   # execute_on = 'CUSTOM'
  #   sync_only = true
  #   sync_times = '0.0 0.3 0.6 0.9 1.2 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5
  #   4.8 5.1 5.4 5.7 6.0 6.3 6.6 6.9 7.2 7.5 7.8 8.1 8.4 8.7 9.0 9.3 9.6
  #   9.9 10.2 10.5 10.8 11.1 11.4 11.7 12.0 12.3 12.6 12.9 13.2 13.5 13.8 14.1
  #   14.4 14.7 15.0 15.3 15.6 15.9 16.2 16.5 16.8 17.1 17.4 17.7 18.0 18.3 18.6 18.9
  #   19.2 19.5 19.8 20.1 20.4 20.7 21.0 21.3 21.6 21.9 22.2 22.5 22.8 23.1 23.4 23.7
  #   24.0 24.3 24.6 24.9 25.2 25.5 25.8 26.1 26.4 26.7 27.0 27.3 27.6 27.9 28.2 28.5
  #   28.8 29.1 29.4 29.7 30.0'
  # []
[]
