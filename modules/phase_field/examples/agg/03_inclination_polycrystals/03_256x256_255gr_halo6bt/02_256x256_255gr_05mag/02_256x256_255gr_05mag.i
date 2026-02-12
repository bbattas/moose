##############################################################################
# File: 02_256x256_255gr_05mag.i
# File Location: /examples/agg/03_inclination_polycrystals/03_256x256_255gr_halo6bt/02_256x256_255gr_05mag
# Created Date: Thursday February 12th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Thursday February 12th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  Using bt instead of jp, adding more OPs (20->30), and upping halo (2->6)
#  Should resolve the issue of remapping wrong and some grains just instnatly
#   combining when they meet
#  120 or 150 cores is fine, ~1.9m dofs
##############################################################################

i_tol = 0
a_tol = 5
amag = 0.05

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2 # Problem dimension
    nx = 256 # Number of elements in the x-direction
    ny = 256 # Number of elements in the y-direction
    xmin = 0 # minimum x-coordinate of the mesh
    xmax = 256 # maximum x-coordinate of the mesh
    ymin = 0 # minimum y-coordinate of the mesh
    ymax = 256 # maximum y-coordinate of the mesh
    # elem_type = QUAD4 # Type of elements used in the mesh
    # uniform_refine = 3 # Initial uniform refinement of the mesh
  []
  parallel_type = DISTRIBUTED # Periodic BCs
  second_order = false
  # uniform_refine = 1
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 30 #20 # Number of order parameters used
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
      polycrystal_ic_uo = voronoi
    []
  []
[]

[UserObjects]
  [voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = bt #jp #bt
    file_name = '../../02_256x256_255gr/00_sub/2D_256x256_255gr_ctrs.txt'
    int_width = 6
  []
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
  [halos]
    order = CONSTANT
    family = MONOMIAL
  []
  [sum_inc]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
    variable_mobility = false
    # order = SECOND
  []
  [PolycrystalInclinationKernel]
    variable_mobility = false
    hgb_mask = hgb
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
  [incl_vec]
    type = GGInclinationVector
    # grain_tracker = grain_tracker
    output_properties = 'inclination_vector ang_dist'
    gb_id_method = SWITCH
    # ffc = gr_flood_uo
    grain_tracker = grain_tracker
    hgb = hgb
    hgb_threshold = 0.5
    outputs = exodus
  []
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L0         gamma_iso kappa_op  iw_iso gbe_iso'
    prop_values = '1.15382e-6   1.5    2.07337e7   6   4.60748e6'
  []
  # [constants_2]
  #   type = GenericConstantMaterial
  #   prop_names = 'L        '#gamma_asymm'# int_width'#   mu'
  #   prop_values = '1.15382e-6 '#' 1.5 '#'    6        '#5.521269e6'
  # []
  # [GBInc_Base]
  #   type = GBInclinationBase
  #   gb_id_method = ffc
  #   ffc = gr_flood_uo
  #   angular_func = ATAN_2D
  #   intol = ${i_tol} #200 # cut if alpha > intol
  #   altol = ${a_tol} #10 #1.5 # cut if h*alpha > altol
  #   output_properties = 'theta_ij gtnum dtheta_dgradeta d2theta_dgradeta2'
  #   outputs = 'exodus'
  # []
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
    ifunc_a = ${amag}
    ifunc_b = 2
    ifunc_c = 0
    ifunc_d = 0
    # L
    aniso_L = false
    # Other Properties
    gb_energy_iso_name = gbe_iso
    kappa = kappa_op
    free_energy_m = 5.521269e6
    output_properties = 'gtnum int_width gamma_asymm' # testout3'
    # testout1 testout2 testoutgrad testouttens
    outputs = 'exodus'
  []
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    hgb_threshold = 0
    func_type = COMBINED
    output_properties = 'hgb'
    outputs = 'exodus'
  []
  # [gamma_test]
  #   type = ElementalGammaMaterial
  #   gb_energy_iso_name = gbe_iso
  #   kappa = kappa_op
  #   free_energy_m = 5.521269e6
  #   well = false
  #   output_properties = 'gamma_asymm int_noij int_width'
  #   outputs = 'exodus'
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
  [avg_iw]
    type = ElementAverageMaterialProperty
    mat_prop = int_width
  []
  [avg_gamma]
    type = ElementAverageMaterialProperty
    mat_prop = gamma_asymm
  []
  # [gr0]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr0
  # []
  # [gr1]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr1
  # []
  # [gr2]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr2
  # []
  # [gr3]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr3
  # []
  # [gr4]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr4
  # []
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
  # end_time = 20
  # dtmin = 0.1
  # end_time = 1000000.0
  # num_steps = 2
  automatic_scaling = true
  compute_scaling_once = false
  # dt = 0.05
  # dtmax = 0.5
  # dt = 2e-5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1 #0.001
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
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final' # Default is "final"
    level = 2 # Default is 1
    heaviest_branch = true # Default is false
    heaviest_sections = 7 # Default is 0
  []
  # file_base = 20_15gr_aniso${amag}_a${a_tol}_i${i_tol}
  # file_base = 17_bicr_large_withnewKernel
  # file_base = test
[]
