##############################################################################
# File: MLPaperCase3.i
# File Location: /examples/agg/01_initial_testing/16_complex_alpha/02_vishal_iso
# Created Date: Friday July 18th 2025
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Friday July 18th 2025
# Modified By: Brandon Battas
# -----
# Description:
#  This simulation predicts GB migration of a 2D copper polycrystal with
#   100 grains represented with 18 order parameters
#  Mesh adaptivity and time step adaptivity are used
#  An AuxVariable is used to calculate the grain boundary locations
#  Postprocessors are used to record time step and the number of grains
#  2 nodes of 32/30 with 8gb/cpu according to vishal (64/512gb)
# Yeah thats not right it needs like 10 nodes of 30 not 1 of 64
# 8 is a comfy 12.5k ish dofs per core?
##############################################################################

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  [gmg]
     type = DistributedRectilinearMeshGenerator
     dim = 2 # Problem dimension
     nx = 512 # Number of elements in the x-direction
     ny = 512 # Number of elements in the y-direction
     xmin = 0    # minimum x-coordinate of the mesh
     xmax = 512 # maximum x-coordinate of the mesh
     ymin = 0    # minimum y-coordinate of the mesh
     ymax = 512 # maximum y-coordinate of the mesh
     elem_type = QUAD4  # Type of elements used in the mesh
     # uniform_refine = 3 # Initial uniform refinement of the mesh
   []
  parallel_type = DISTRIBUTED # Periodic BCs
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 12 # Number of order parameters used
  var_name_base = gr # Base name of grains
  # T = 1400 # Constant temperature of the simulation (for mobility calculation)
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [./PolycrystalVariables]
  [../]
[]

[UserObjects]
  [./voronoi]
    type = PolycrystalVoronoi
    # grain_num = 512 # Number of grains
    # rand_seed = 2460
     # output_adjacency_matrix = true
    coloring_algorithm = bt
    file_name = 'Case3_grains_centers.txt'
    # int_width = 6
  [../]
  [./term]
    type = Terminator
    expression = 'grain_tracker < 25'
  [../]
  [./grain_tracker]
    type = GrainTracker
    threshold = 0.2
    connecting_threshold = 0.08
    compute_halo_maps = true # Only necessary for displaying HALOS
    halo_level = 6
    compute_var_to_feature_map = true

  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    [../]
  [../]
[]

[AuxVariables]
  # Dependent variables
  # [./GBEnergyT]
  # [../]
  [./bnds]
    # Variable used to visualize the grain boundaries in the simulation
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ghost_regions]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./halos]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [./PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
     variable_mobility = false
  [../]
[]

[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  # [./GBEngeryTCal]
  #   type = FunctionAux
  #   variable = GBEnergyT
  #   function = GBEnergyTfun
  # [../]
  [./bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
  [../]
  [./ghosted_entities]
    type = FeatureFloodCountAux
    variable = ghost_regions
    flood_counter = grain_tracker
    field_display = GHOSTED_ENTITIES
    execute_on = 'initial timestep_end'
  [../]
  [./halos]
    type = FeatureFloodCountAux
    variable = halos
    flood_counter = grain_tracker
    field_display = HALOS
    execute_on = 'initial timestep_end'
  [../]
[]


[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y'
    [../]
  [../]
[]



# [Functions]
#   [./GBEnergyTfun]
#     type = ParsedFunction
#     value = '1.56-5.87*T'
#   [../]
# []
[Materials]
  [./UO2]
    # Material properties
    type = GBEvolution
    T = 1400 # Constant temperature of the simulation (for mobility calculation)
    wGB = 6 # Width of the diffuse GB
    GBMobility = 3.2407e-11 #m^4(Js) for copper from Schoenfelder1997
    Q = 3.01 #eV for copper from Schoenfelder1997
    GBenergy = 0.7382 #0.708 #J/m^2 from Schoenfelder1997
    length_scale = 1e-06
    time_scale = 1
    # output_properties = 'M_GB sigma kappa_op gamma_asymm mu L'
    # outputs = 'exodus'
  [../]
  [incl_vec]
    type = GGInclinationVector
    # grain_tracker = grain_tracker
    output_properties = 'inclination_vector ang_dist'
    gb_id_method = HGB
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
  [./dt]
    type = TimestepSize
  [../]
  [./n_elements]
    type = NumElements
    execute_on = 'initial timestep_end'
  [../]
  [./n_nodes]
    type = NumNodes
    execute_on = 'initial timestep_end'
  [../]
  [./DOFs]
    type = NumDOFs
    execute_on = 'initial timestep_end'
  [../]
  [./avg_grain_volumes]
   type = AverageGrainVolume
   feature_counter = grain_tracker
   execute_on = 'initial timestep_end'
   [../]
[]

[VectorPostprocessors]
  [./grain_volumes]
    type = FeatureVolumeVectorPostprocessor
    flood_counter = grain_tracker
    execute_on = 'initial timestep_end'
    # output_centroids = true
  [../]
# [./mem]
#   type = VectorMemoryUsage
#   execute_on = 'INITIAL TIMESTEP_END NONLINEAR LINEAR'
#   report_peak_value = true
#   #mem_units = kilobytes # or bytes, megabytes, gigabytes
#   mem_units = gigabytes
# [../]

[]

[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  petsc_options_value = 'hypre boomeramg 101 ds'

  l_max_its = 30 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 40 # Max number of nonlinear iterations
  nl_rel_tol = 1e-8 #1e-10 # Absolute tolerance for nonlienar solves
  # nl_abs_tol = 1e-10


  start_time = 0.0
  dtmin = 0.1
  end_time = 1000000.0
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    cutback_factor = 0.9
    growth_factor = 1.1
    optimal_iterations = 8
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
  csv = true
  checkpoint = true
  [./console]
    type = Console
    max_rows = 20 # Will print the 20 most recent postprocessor values to the screen
  [../]
[]
