##############################################################################
# File: 01_iso_iso.i
# File Location: /examples/agg/05_large_bicr_L/10_miso_iso/01_iso_iso
# Created Date: Monday March 9th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Monday March 9th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  Isotropic reference for the iso material version of the aniso low/high miso
#
#
#
##############################################################################

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = '../../00_sub/bicr_200_r80.txt'
  []
  # [gmg]
  #   type = DistributedRectilinearMeshGenerator
  #   dim = 2
  #   nx = 256
  #   ny = 256
  #   xmin = 0
  #   xmax = 256
  #   ymin = 0
  #   ymax = 256
  # []
  parallel_type = DISTRIBUTED # Periodic BCs
  # second_order = false
  # uniform_refine = 1
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 2 #15 # Number of order parameters used
  var_name_base = gr # Base name of grains
  # T = 1400 # Constant temperature of the simulation (for mobility calculation)
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [PolycrystalVariables]
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
    []
  []
[]

[UserObjects]
  [ebsd_reader]
    type = EBSDReader
  []
  [ebsd]
    type = PolycrystalEBSD
    coloring_algorithm = bt #jp
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
  []
  [grain_tracker]
    type = GrainTracker
    # variable = 'gr0 gr1 gr2 gr3 gr4'
    threshold = 0.01
    connecting_threshold = 0.01
    compute_var_to_feature_map = true
    compute_halo_maps = false # Only necessary for displaying HALOS
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
  [contour]
  []
[]

[Kernels]
  [PolycrystalKernel]
  []
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
    prop_names = 'L0             kappa_input  int_width_iso   const_m    sigma0 sigma_low  sigma_high'
    prop_values = '1.15382e-6  2.07337e7        6          5.521269e6  4.60748e6  2.6293e6   4.5612e6'
  []
  [gbe_backconvert]
    type = GenericConstantMaterial
    prop_names = 'gbej_iso      gbej_low     gbej_high'
    prop_values = '7.38200e-13  4.21260e-13  7.30785e-13'
  []
  # [constants2]
  #   type = GenericConstantMaterial
  #   prop_names = 'int_width          gamma_asymm     sigma      L'
  #   prop_values = '10.740284831071 0.64666592158764 2.6293e6 6.584442e-07 ' #1.15382e-6'
  # []
  [GBMat]
    type = GBEvolution
    T = 1300 # shouldnt actually be used
    wGB = 6
    GBMobility = 3.2407e+13
    length_scale = 1
    time_scale = 1
    GBenergy = 7.38200e-13
    output_properties = 'L l_GB gamma_asymm sigma M_GB kappa_op'
    outputs = 'exodus'
  []
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    func_type = COMBINED
    hgb_threshold = 0
    output_properties = 'hgb'
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
    mat_prop = l_GB #int_width
    value_type = MIN
  []
  [iw_max]
    type = ElementExtremeMaterialProperty
    mat_prop = l_GB #int_width
    value_type = MAX
  []
[]

[VectorPostprocessors]
  [radial]
    type = LineValueSampler
    end_point = '200 100 0'
    num_points = 100
    sort_by = X
    start_point = '100 100 0'
    variable = 'gr0 gr1 contour'
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
  end_time = 120

  automatic_scaling = true
  compute_scaling_once = false

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5 #0.001
    # cutback_factor = 0.9
    # growth_factor = 1.1
    optimal_iterations = 6
    linear_iteration_ratio = 1e5
  []
[]

[Outputs]
  # file_base = MLPaperActualSimulations/Case3
  exodus = true # Exodus file will be outputted
  #nemesis = true
  console = true
  csv = true
  checkpoint = true
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
  # file_base = 20_15gr_aniso${amag}_a${a_tol}_i${i_tol}
  # file_base = 17_bicr_large_withnewKernel
  # file_base = 06_normal
[]
