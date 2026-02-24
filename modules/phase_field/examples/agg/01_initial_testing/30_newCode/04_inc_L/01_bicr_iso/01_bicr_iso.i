##############################################################################
# File: 01_bicr_iso.i
# File Location: /examples/agg/01_initial_testing/30_newCode/04_inc_L/01_bicr_iso
# Created Date: Monday February 23rd 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Monday February 23rd 2026
# Modified By: Brandon Battas
# -----
# Description:
#
#
#
#
##############################################################################

amag = 0.0
a_tol = 5
i_tol = 0

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 120
    ny = 120
    xmin = 0
    xmax = 120
    ymin = 0
    ymax = 120
  []
  parallel_type = DISTRIBUTED
  # uniform_refine = 1
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 2 #15 # Number of order parameters used
  var_name_base = gr # Base name of grains
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [PolycrystalVariables]
  []
[]

[UserObjects]
  [grain_tracker]
    type = GrainTracker
    # variable = 'gr0 gr1 gr2 gr3 gr4'
    threshold = 0.01
    connecting_threshold = 0.01
    compute_var_to_feature_map = true
    compute_halo_maps = true # Only necessary for displaying HALOS
    execute_on = 'initial timestep_end'
    # use_less_than_threshold_comparison = true
  []
  # [term]
  #   type = Terminator
  #   expression = 'grain_tracker < 2'
  # []
[]

[ICs]
  [gr0_IC]
    type = SmoothCircleIC
    invalue = 1
    outvalue = 0
    radius = 45
    variable = gr0
    x1 = 60
    y1 = 60
    int_width = 9 #6
  []
  [gr1_IC]
    type = SmoothCircleIC
    invalue = 0
    outvalue = 1
    radius = 45
    variable = gr1
    x1 = 60
    y1 = 60
    int_width = 9 #6
  []
[]

[AuxVariables]
  [bnds]
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [halos]
    order = CONSTANT
    family = MONOMIAL
  []
  [sum_inc]
    order = CONSTANT
    family = MONOMIAL
  []
  [contour]
  []
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
    variable_mobility = false
  []
  # [PolycrystalInclinationKernel]
  #   variable_mobility = false
  # []
[]

[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  [bnds_aux]
    # AuxKernel that calculates the GB term
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
  [halos]
    type = FeatureFloodCountAux
    variable = halos
    flood_counter = grain_tracker
    field_display = HALOS
    execute_on = 'initial timestep_end'
  []
  [sum_inc]
    type = MatVectorijSum
    property = inclination
    variable = sum_inc
  []
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
    prop_names = 'L0             kappa_op  int_width_iso  gbe_iso' #'gamma_asymm'
    prop_values = '1.15382e-6  2.07337e7        6       4.60748e6' #'1.5    '
  []
  # [constants2]
  #   type = GenericConstantMaterial
  #   prop_names = 'int_width  gamma_asymm  sigma      L     mu  '
  #   prop_values = '6.0208     1.2110    4.6075e6 1.153820e-6  5.521269e6
  # []
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    func_type = COMBINED
    hgb_threshold = 0
    output_properties = 'hgb'
    outputs = 'exodus'
  []
  [GBInc_cos]
    type = GBInclination
    gb_id_method = GRAINTRACKER
    # ffc = gr_flood_uo
    grain_tracker = grain_tracker
    angular_func = ATAN_2D
    intol = ${i_tol} #100 #10 #100 #200 # cut if alpha > intol
    altol = ${a_tol} #100 #10 #1.5 # cut if h*alpha > altol
    limit_umag = true
    # Inclination function
    inc_func = COS
    ifunc_a = ${amag} #0.3
    ifunc_b = 2
    ifunc_c = 0
    ifunc_d = 0
    # L
    combine_gb_form = weighted
    aniso_L = false
    noDeriv_L = false
    # Other Properties
    gb_energy_iso_name = gbe_iso
    kappa = kappa_op
    free_energy_m = 5.521269e6
    output_properties = 'int_width gamma_asymm L'
    outputs = 'exodus'
    # # testout1 testout2 testoutgrad testouttens
    # outputs = 'out'
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
  [max_mpi_memory]
    type = MemoryUsage
    value_type = max_process
    report_peak_value = True
    mem_units = MEGABYTES
    execute_on = 'NONLINEAR TIMESTEP_END'
  []
  [avg_iw]
    type = ElementAverageMaterialProperty
    mat_prop = int_width
  []
  [avg_gamma]
    type = ElementAverageMaterialProperty
    mat_prop = gamma_asymm
  []
  [gr0]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  []
  [gr1]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  []
[]

[VectorPostprocessors]
  [horizontal]
    type = LineValueSampler
    end_point = '120 60 0'
    num_points = 120
    sort_by = X
    start_point = '60 60 0'
    variable = 'gr0 gr1 contour'
    outputs = 'vpp'
  []
  [vertical]
    type = LineValueSampler
    end_point = '60 120 0'
    num_points = 120
    sort_by = Y
    start_point = '60 60 0'
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
  end_time = 40
  # dt = 0.1
  # dtmin = 0.1
  # end_time = 1000000.0
  # num_steps = 10
  automatic_scaling = true
  compute_scaling_once = false
  # line_search = NONE
  # dt = 0.05
  # dtmax = 0.5
  dt = 0.4
  # [TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 0.1
  #   # cutback_factor = 0.9
  #   # growth_factor = 1.1
  #   optimal_iterations = 6
  #   linear_iteration_ratio = 1e5
  # []

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
  # exodus = true # Exodus file will be outputted
  [exodus]
    type = Exodus
    # time_step_interval = 4
    sync_times = '0.0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6 6.0
    6.4 6.8 7.2 7.6 8.0 8.4 8.8 9.2 9.6 10.0 10.4 10.8 11.2 11.6 12.0 12.4 12.8
    13.2 13.6 14.0 14.4 14.8 15.2 15.6 16.0 16.4 16.8 17.2 17.6 18.0 18.4 18.8 19.2
    19.6 20.0 20.4 20.8 21.2 21.6 22.0 22.4 22.8 23.2 23.6 24.0 24.4 24.8 25.2 25.6
    26.0 26.4 26.8 27.2 27.6 28.0 28.4 28.8 29.2 29.6 30.0 30.4 30.8 31.2 31.6 32.0
    32.4 32.8 33.2 33.6 34.0 34.4 34.8 35.2 35.6 36.0 36.4 36.8 37.2 37.6 38.0 38.4
    38.8 39.2 39.6 40.0'
  []
  [vpp]
    type = CSV
    # time_step_interval = 4
    sync_times = '0.0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2 3.6 4.0 4.4 4.8 5.2 5.6 6.0
    6.4 6.8 7.2 7.6 8.0 8.4 8.8 9.2 9.6 10.0 10.4 10.8 11.2 11.6 12.0 12.4 12.8
    13.2 13.6 14.0 14.4 14.8 15.2 15.6 16.0 16.4 16.8 17.2 17.6 18.0 18.4 18.8 19.2
    19.6 20.0 20.4 20.8 21.2 21.6 22.0 22.4 22.8 23.2 23.6 24.0 24.4 24.8 25.2 25.6
    26.0 26.4 26.8 27.2 27.6 28.0 28.4 28.8 29.2 29.6 30.0 30.4 30.8 31.2 31.6 32.0
    32.4 32.8 33.2 33.6 34.0 34.4 34.8 35.2 35.6 36.0 36.4 36.8 37.2 37.6 38.0 38.4
    38.8 39.2 39.6 40.0'
  []
  console = true
  csv = true
  checkpoint = false
  # [checkpoint]
  #   type = Checkpoint
  #   num_files = 5
  # []
  # [console]
  #   type = Console
  #   max_rows = 20 # Will print the 20 most recent postprocessor values to the screen
  # []
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final' # Default is "final"
    level = 2 # Default is 1
    heaviest_branch = true # Default is false
    heaviest_sections = 7 # Default is 0
  []
  # file_base = 01_3op_tols_amag${amag}_i${i_tol}_a${a_tol}
  # file_base = 17_bicr_large_withnewKernel
  # file_base = test
[]
