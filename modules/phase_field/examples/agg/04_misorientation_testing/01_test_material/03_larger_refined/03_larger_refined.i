##############################################################################
# File: 03_larger_refined.i
# File Location: /examples/agg/04_misorientation_testing/01_test_material/03_larger_refined
# Created Date: Tuesday February 10th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Tuesday February 10th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  using a case with elements at 0.5 instead of 1 to get more refinement
#
#
#
##############################################################################

# i_tol = 0
# a_tol = 5
# amag = 0.0

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = '../00_sub/2D_200x200_25gr_fid0.txt'
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
  # parallel_type = DISTRIBUTED # Periodic BCs
  # second_order = false
  # uniform_refine = 1
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 12 #15 # Number of order parameters used
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
    compute_halo_maps = true # Only necessary for displaying HALOS
    execute_on = 'initial timestep_end'
    halo_level = 2
    # use_less_than_threshold_comparison = true
    polycrystal_ic_uo = ebsd
  []
  # [term]
  #   type = Terminator
  #   expression = 'grain_tracker < 60'
  # []
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
  # [sum_inc]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  [aphi1]
    order = CONSTANT
    family = MONOMIAL
  []
  [bPhi]
    order = CONSTANT
    family = MONOMIAL
  []
  [cphi2]
    order = CONSTANT
    family = MONOMIAL
  []
  [ebsd_numbers]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
    # variable_mobility = false
    # order = SECOND
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
  [halos]
    type = FeatureFloodCountAux
    variable = halos
    flood_counter = grain_tracker
    field_display = HALOS
    execute_on = 'initial timestep_end'
  []
  # [sum_inc]
  #   type = MatVectorijSum
  #   property = inclination
  #   variable = sum_inc
  # []
  # The phi will output the Euler angle from EBSD data, and the data structure
  # will change with the guide from grain_tracker
  [aphi1]
    type = OutputEulerAngles
    variable = aphi1
    euler_angle_provider = ebsd_reader
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [bPhi]
    type = OutputEulerAngles
    variable = bPhi
    euler_angle_provider = ebsd_reader
    grain_tracker = grain_tracker
    output_euler_angle = 'Phi'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [cphi2]
    type = OutputEulerAngles
    variable = cphi2
    euler_angle_provider = ebsd_reader
    grain_tracker = grain_tracker
    output_euler_angle = 'phi2'
    execute_on = 'INITIAL TIMESTEP_END'
  []
  # Import the unique grain ID from ebsd data, and the data structure
  # will change with the guide from grain_tracker
  [ebsd_numbers]
    type = EBSDReaderAvgDataAux
    data_name = feature_id
    ebsd_reader = ebsd_reader
    grain_tracker = grain_tracker
    variable = ebsd_numbers
    execute_on = 'initial timestep_end'
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
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L             kappa_op  int_width_iso   mu       sigma' #'gamma_asymm'
    prop_values = '1.15382e-6  2.07337e7        6      5.521269e6  4.60748e6' #'1.5    '
  []
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    grain_ops = 'gr0 gr1 gr2'
    hgb_threshold = 0
    output_properties = 'hgb'
    outputs = 'exodus'
  []
  #
  [GBM]
    type = GBMisorientation
    ebsd_reader = ebsd_reader
    grain_tracker = grain_tracker
    kappa = kappa_op
    output_properties = 'eul_a eul_b eul_c quat_a quat_b quat_c quat_d quat_mag
    misorientation miso_axis_polar miso_axis_azimuth other_out miso_ang_energy miso_ax_energy
    twist_energy tilt_energy f_miso gamma_asymm int_width'
    outputs = exodus
  []
  [gbe]
    type = ParsedMaterial
    property_name = gbe
    # coupled_variables = 'gr0 gr1 gr2'
    material_property_names = 'twist_energy tilt_energy'
    expression = '0.3 + twist_energy'
    outputs = 'exodus'
  []
  # [GB_type]
  #   # The new developed Miso Bnds Aux Kernel
  #   type = ComputeGBMisorientationType
  #   ebsd_reader = ebsd_reader
  #   grain_tracker = grain_tracker
  #   output_properties = 'gb_type'
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
  end_time = 20
  # dtmin = 0.1
  # end_time = 1000000.0
  # num_steps = 5
  dt = 0.4
  automatic_scaling = true
  compute_scaling_once = false
  # dt = 0.05
  # dtmax = 0.5
  # dt = 2e-5
  # [TimeStepper]
  #   type = IterationAdaptiveDT
  #   dt = 1 #0.001
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
  # file_base = test
[]
