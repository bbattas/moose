##############################################################################
# File: 06_lowgt_bicr.i
# File Location: /examples/agg/01_initial_testing/28_VECTOR_inclination/06_lowgt_bicr
# Created Date: Thursday January 15th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Thursday January 15th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  BICR input corrected using the lower grain tracker threshold to actually
#   produce anisotropy
#
#
##############################################################################

a_tol = 0
i_tol = 0

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2 # Problem dimension
    nx = 160 # Number of elements in the x-direction
    ny = 160 # Number of elements in the y-direction
    xmin = 0 # minimum x-coordinate of the mesh
    xmax = 160 # maximum x-coordinate of the mesh
    ymin = 0 # minimum y-coordinate of the mesh
    ymax = 160 # maximum y-coordinate of the mesh
    # elem_type = QUAD4 # Type of elements used in the mesh
    # uniform_refine = 3 # Initial uniform refinement of the mesh
  []
  parallel_type = DISTRIBUTED # Periodic BCs
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
  [term]
    type = Terminator
    expression = 'grain_tracker < 2'
  []
[]

[ICs]
  [gr0_IC]
    type = SmoothCircleIC
    invalue = 1
    outvalue = 0
    radius = 60
    variable = gr0
    x1 = 80
    y1 = 80
    int_width = 6
  []
  [gr1_IC]
    type = SmoothCircleIC
    invalue = 0
    outvalue = 1
    radius = 60
    variable = gr1
    x1 = 80
    y1 = 80
    int_width = 6
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
  []
  [gr0_ACInc]
    type = ACInterfaceInclinationGamma
    variable = gr0
    coupled_variables = 'gr1'
    debug_kernel = false
    skip_off = false
    variable_L = false
    mask_name = hgb
  []
  [gr1_ACInc]
    type = ACInterfaceInclinationGamma
    variable = gr1
    coupled_variables = 'gr0'
    debug_kernel = false
    skip_off = false
    variable_L = false
    mask_name = hgb
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
    # outputs = exodus
  []
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L0         gamma_iso kappa_op  iw_iso gbe_iso'
    prop_values = '1.15382e-6   1.5    2.07337e7   6   4.60748e6'
  []
  [GBInc_cos]
    type = GBInclination
    gb_id_method = GRAINTRACKER
    # ffc = gr_flood_uo
    grain_tracker = grain_tracker
    angular_func = ATAN_2D
    intol = ${i_tol} #100 #10 #100 #200 # cut if alpha > intol
    altol = ${a_tol} #100 #10 #1.5 # cut if h*alpha > altol
    # Inclination function
    inc_func = COS
    ifunc_a = 0.3
    ifunc_b = 2
    ifunc_c = 0
    ifunc_d = 0
    # L
    aniso_L = false
    # Other Properties
    gb_energy_iso_name = gbe_iso
    kappa = kappa_op
    free_energy_m = 5.521269e6
    # output_properties = 'gtnum no_ij_pairs gamma_a'
    # # testout1 testout2 testoutgrad testouttens
    # outputs = 'exodus'
  []
  [gamma_test]
    type = ElementalGammaMaterial
    gb_energy_iso_name = gbe_iso
    kappa = kappa_op
    free_energy_m = 5.521269e6
    well = false
    output_properties = 'gamma_asymm int_noij int_width'
    outputs = 'out' # h5out nemesis h5nemesis'
  []
  [hgb_a]
    type = ParsedMaterial
    property_name = hgb_a
    coupled_variables = 'gr0 gr1' # gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10' # gr11 gr12 gr13 gr14 gr15'
    expression = 'hb:=gr0*gr0 + gr1*gr1;
                  4 * (1 - hb) * (1 - hb)'
    # hg:=4 * (1 - hb) * (1 - hb);
    # if(hg>1.0,1.0,hg)'
    # outputs = 'exodus'
    #  + gr11*gr11 + gr12*gr12 + gr13*gr13 + gr14*gr14 + gr15*gr15
  []
  [hgb_b]
    type = SwitchingFunctionGBMaterial
    h_name = hgb_b
    grain_ops = 'gr0 gr1' # gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10' # gr11 gr12 gr13 gr14 gr15'
    hgb_threshold = 0
    output_properties = 'hgb_b'
    # outputs = 'exodus'
  []
  [hgb]
    type = ParsedMaterial
    property_name = hgb
    coupled_variables = 'gr0 gr1' # gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10' # gr11 gr12 gr13 gr14 gr15'
    material_property_names = 'hgb_a hgb_b'
    # expression = 'h1:=if(hgb_a>1,1,if(hgb_a<0,0.0,hgb_a));
    #               h2:=if(hgb_b>1,1,if(hgb_b<0,0.0,hgb_b));
    #               h3:=(h1 + h2);
    #               if(h3>1,1,if(h3<0,0.0,h3))'
    expression = 'h3:=(hgb_a + hgb_b)/2;
                  if(h3>1,1,if(h3<0,0.0,h3))'
    # outputs = 'exodus' #h3:=h1+h2;
  []
  [sumgr]
    type = ParsedMaterial
    property_name = sumgr
    coupled_variables = 'gr0 gr1' # gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10' # gr11 gr12 gr13 gr14 gr15'
    expression = 'gr0 + gr1' # + gr2 + gr3 + gr4 + gr5 + gr6 + gr7 + gr8 + gr9 + gr10'
    #  + gr11 + gr12 + gr13 + gr14 + gr15'
    outputs = 'out' # h5out nemesis h5nemesis'
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
  [tot_gr_op]
    type = ElementIntegralMaterialProperty
    mat_prop = sumgr
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
  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -mat_mffd_type'
  # petsc_options_value = 'hypre boomeramg 101 ds'
  petsc_options_iname = '-pc_type -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm        preonly       lu           2'

  l_max_its = 90 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 12 # Max number of nonlinear iterations
  nl_rel_tol = 1e-7 #1e-10 # Relative tolerance for nonlienar solves
  # nl_abs_tol = 1e-10

  start_time = 0.0
  end_time = 30
  # dtmin = 0.1
  # end_time = 1000000.0
  # num_steps = 10
  automatic_scaling = true
  compute_scaling_once = false
  # dt = 0.05
  # dtmax = 0.5
  # dt = 0.01
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    cutback_factor = 0.9
    growth_factor = 1.1
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
  # exodus = false # Exodus file will be outputted
  [out]
    type = Exodus
  []
  # [h5out]
  #   type = Exodus
  #   write_hdf5 = true
  # []
  # [nemesis]
  #   type = Nemesis
  # []
  # [h5nemesis]
  #   type = Nemesis
  #   write_hdf5 = true
  # []
  #nemesis = true
  console = true
  csv = false
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
  # file_base = 12_halt_hypre_nopre_i${i_tol}_a${a_tol}
  # file_base = 17_bicr_large_withnewKernel
  # file_base = test
[]
