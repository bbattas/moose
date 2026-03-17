##############################################################################
# File: 09_custom_noStiff.i
# File Location: /examples/agg/05_large_bicr_L/16_inc_more_aspect/09_custom_noStiff
# Created Date: Monday March 16th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Tuesday March 17th 2026
# Modified By: bbattas
# -----
# Description:
#  custom ellipse IC for no stiffness bicrystal
#
#
#
##############################################################################

i_tol = 0
a_tol = 5
amag = 0.0
yic = 80
xic = 80

[Mesh]
  # [ebsd_mesh]
  #   type = EBSDMeshGenerator
  #   filename = '../../00_sub/bicr_200_r80.txt'
  # []
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 200
    ny = 200
    xmin = 0
    xmax = 200
    ymin = 0
    ymax = 200
  []
  parallel_type = DISTRIBUTED
  # second_order = true
  # uniform_refine = 1
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [PolycrystalVariables]
    # order = SECOND
  []
[]

[ICs]
  # [PolycrystalICs]
  #   [PolycrystalColoringIC]
  #     polycrystal_ic_uo = ebsd
  #   []
  # []
  [gr0_IC]
    type = SmoothSuperellipsoidIC
    invalue = 1
    outvalue = 0
    variable = gr0
    a = ${xic}
    b = ${yic}
    n = 2
    x1 = 100
    y1 = 100
  []
  [gr1_IC]
    type = SmoothSuperellipsoidIC
    invalue = 0
    outvalue = 1
    variable = gr1
    a = ${xic}
    b = ${yic}
    n = 2
    x1 = 100
    y1 = 100
  []
[]

[UserObjects]
  # [ebsd_reader]
  #   type = EBSDReader
  # []
  # [ebsd]
  #   type = PolycrystalEBSD
  #   coloring_algorithm = bt #jp
  #   ebsd_reader = ebsd_reader
  #   enable_var_coloring = true
  # []
  [grain_tracker]
    type = GrainTracker
    # variable = 'gr0 gr1 gr2 gr3 gr4'
    threshold = 0.01
    connecting_threshold = 0.01
    compute_var_to_feature_map = true
    compute_halo_maps = true # Only necessary for displaying HALOS
    execute_on = 'initial timestep_end'
    halo_level = 6
    # use_less_than_threshold_comparison = true
    # polycrystal_ic_uo = ebsd
  []
  [term]
    type = Terminator
    expression = 'grain_tracker < 2'
  []
[]

[AuxVariables]
  [bnds]
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [ebsd_numbers]
    order = CONSTANT
    family = MONOMIAL
  []
  [contour]
  []
[]

[Kernels]
  [PolycrystalKernel]
    variable_mobility = false
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
  # Import the unique grain ID from ebsd data, and the data structure
  # will change with the guide from grain_tracker
  # [ebsd_numbers]
  #   type = EBSDReaderAvgDataAux
  #   data_name = feature_id
  #   ebsd_reader = ebsd_reader
  #   grain_tracker = grain_tracker
  #   variable = ebsd_numbers
  #   execute_on = 'initial timestep_end'
  # []
  [contour]
    type = ParsedAux
    variable = contour
    coupled_variables = 'gr0 gr1'
    expression = 'gr0 * gr0 / ((gr0 * gr0) + (gr1 * gr1))'
    execute_on = 'initial timestep_end'
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
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L0             kappa_op  int_width_iso   const_m      sigma0' #'gamma_asymm'
    prop_values = '1.15382e-6  2.07337e7        6      5.521269e6  4.60748e6' #'1.5    '
  []
  # [constants2]
  #   type = GenericConstantMaterial
  #   prop_names = 'int_width  gamma_asymm  sigma      L'
  #   prop_values = '6.0208     1.2110    4.6075e6 1.153820e-6' #1.15382e-6'
  # []
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    func_type = COMBINED
    hgb_threshold = 0
    output_properties = 'hgb'
    outputs = 'exodus'
  []
  [inc_mat]
    type = GBInclination
    gb_id_method = GRAINTRACKER
    grain_tracker = grain_tracker
    angular_func = ATAN_2D
    intol = ${i_tol}
    altol = ${a_tol}
    limit_umag = true
    # Inclination function
    inc_func = COS
    ifunc_a = ${amag}
    ifunc_b = 2
    ifunc_c = 0
    ifunc_d = 0
    # L
    combine_gb_form = avg
    aniso_L = false
    noDeriv_L = true
    aniso_gbmob = false
    stiffness = false
    # Other Properties
    gb_energy_iso_name = sigma0
    kappa = kappa_op
    free_energy_m = 5.521269e6
    output_properties = 'int_width gamma_asymm L'
    outputs = 'exodus'
  []
  #
  # [GBM]
  #   type = GBMisorientation
  #   ebsd_reader = ebsd_reader
  #   grain_tracker = grain_tracker
  #   kappa = kappa_op
  #   output_properties = 'eul_a eul_b eul_c quat_a quat_b quat_c quat_d quat_mag
  #   misorientation miso_axis_polar miso_axis_azimuth other_out miso_ang_energy miso_ax_energy
  #   twist_energy tilt_energy f_miso gamma_asymm int_width'
  #   outputs = exodus
  # []
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
  [max_mpi_memory]
    type = MemoryUsage
    value_type = max_process
    report_peak_value = True
    mem_units = MEGABYTES
    execute_on = 'NONLINEAR TIMESTEP_END'
  []
  [gr0_total]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  []
  [gr1_total]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  []
  # [mis_angle_max]
  #   type = ElementExtremeMaterialProperty
  #   mat_prop = misorientation
  #   value_type = MAX
  # []
  # [mis_ax_polar_max]
  #   type = ElementExtremeMaterialProperty
  #   mat_prop = miso_axis_polar
  #   value_type = MAX
  # []
  # [mis_ax_azim_max]
  #   type = ElementExtremeMaterialProperty
  #   mat_prop = miso_axis_azimuth
  #   value_type = MAX
  # []
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

[VectorPostprocessors]
  [horizontal]
    type = LineValueSampler
    end_point = '200 100 0'
    num_points = 200
    sort_by = X
    start_point = '0 100 0'
    variable = 'gr0 gr1 contour'
    outputs = 'vpp'
  []
  [vertical]
    type = LineValueSampler
    end_point = '100 200 0'
    num_points = 200
    sort_by = Y
    start_point = '100 0 0'
    variable = 'gr0 gr1 contour'
    outputs = 'vpp'
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
  end_time = 120 #5
  # dtmin = 0.1
  # end_time = 1000000.0
  # num_steps = 1
  # dt = 1
  automatic_scaling = true
  compute_scaling_once = false
  # dt = 0.05
  # dtmax = 0.5
  # dt = 2e-5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5 #0.001
    # cutback_factor = 0.9
    # growth_factor = 1.1
    optimal_iterations = 6
    linear_iteration_ratio = 1e5
  []

  # start_time = 0.0
  # dt = 0.1
  # end_time = 100000
  #
  # [./TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 0.1 # Initial time step.  In this simulation it changes.
  #   optimal_iterations = 6 # Time step will adapt to maintain this number of nonlinear iterations
  #   cutback_factor = 0.9
  #   growth_factor = 1.1
  # [../]

  # [./Adaptivity]
  #   # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
  #   initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
  #   refine_fraction = 0.7 # Fraction of high error that will be refined
  #   coarsen_fraction = 0.1 # Fraction of low error that will coarsened
  #   max_h_level = 4 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  # [../]
[]

[Outputs]
  # file_base = MLPaperActualSimulations/Case3
  exodus = true # Exodus file will be outputted
  #nemesis = true
  console = true
  # csv = true
  checkpoint = false
  [vpp]
    type = CSV
  []
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 5
  # []
  # [console]
  #   type = Console
  #   max_rows = 20 # Will print the 20 most recent postprocessor values to the screen
  # []
  # [pgraph]
  #   type = PerfGraphOutput
  #   execute_on = 'initial final' # Default is "final"
  #   level = 2 # Default is 1
  #   heaviest_branch = true # Default is false
  #   heaviest_sections = 7 # Default is 0
  # []
  file_base = 09_a${amag}_custom_noStiff/09_a${amag}_ciustom_noStiff_out
[]
