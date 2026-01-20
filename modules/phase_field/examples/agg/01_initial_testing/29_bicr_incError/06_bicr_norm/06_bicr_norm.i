##############################################################################
# File: 06_bicr_norm.i
# File Location: /examples/agg/01_initial_testing/29_bicr_incError/06_bicr_norm
# Created Date: Tuesday January 20th 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Tuesday January 20th 2026
# Modified By: Brandon Battas
# -----
# Description:
#  input 01 the .3 bicr with gt at 0.01, but with the u vector normalized in
#   the incliantion material before all the calculations. should be wrong but
#   i need to double check for my sanity that this wont work (its minimizing the
#   torque term so it should be wrong)
##############################################################################

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
  second_order = false
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
  [gr_flood_uo]
    type = FeatureFloodCount
    variable = 'gr0 gr1'
    threshold = 0.01
    connecting_threshold = 0.01
    compute_var_to_feature_map = true
    compute_halo_maps = true # Only necessary for displaying HALOS
    execute_on = 'initial timestep_end'
    # use_less_than_threshold_comparison = true
  []
  [term]
    type = Terminator
    expression = 'gr_flood_uo < 2'
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
  [gr_unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [gr_halos]
    order = CONSTANT
    family = MONOMIAL
  []
  # [theta_00]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  [inc_01]
    order = CONSTANT
    family = MONOMIAL
  []
  # [inc_02]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [inc_12]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [gamma_01_v]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [gamma_02_v]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [gamma_12_v]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [i_0]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [j_0]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  [theta_01_v]
    order = CONSTANT
    family = MONOMIAL
  []
  # [theta_02_v]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [theta_12_v]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [theta_11]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [theta_01_x]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [theta_01_y]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [theta_01_z]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
    variable_mobility = false
    # order = SECOND
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
  # [gr2_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr2
  #   coupled_variables = 'gr0 gr1'
  #   debug_kernel = true
  #   skip_off = false
  #   variable_L = false
  #   mask_name = hgb
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
  [gr_unique_grains]
    type = FeatureFloodCountAux
    variable = gr_unique_grains
    flood_counter = gr_flood_uo
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  []
  [gr_halos]
    type = FeatureFloodCountAux
    variable = gr_halos
    flood_counter = gr_flood_uo
    field_display = HALOS
    execute_on = 'initial timestep_end'
  []
  # [theta_00]
  #   type = MatVectorComponentAux
  #   variable = theta_00
  #   property = theta_ij
  #   i = 0
  #   j = 0
  #   execute_on = 'initial timestep_end'
  # []
  [inc_01]
    type = MatVectorComponentAux
    variable = inc_01
    property = inclination
    i = 0
    j = 1
    execute_on = 'initial timestep_end'
  []
  # [inc_02]
  #   type = MatVectorComponentAux
  #   variable = inc_02
  #   property = inclination
  #   i = 0
  #   j = 2
  #   execute_on = 'initial timestep_end'
  # []
  # [inc_12]
  #   type = MatVectorComponentAux
  #   variable = inc_12
  #   property = inclination
  #   i = 1
  #   j = 2
  #   execute_on = 'initial timestep_end'
  # []
  # [gamma_01_v]
  #   type = MatVectorComponentAux
  #   variable = gamma_01_v
  #   property = gamma_ij
  #   i = 0
  #   j = 1
  #   execute_on = 'initial timestep_end'
  # []
  # [gamma_02_v]
  #   type = MatVectorComponentAux
  #   variable = gamma_02_v
  #   property = gamma_ij
  #   i = 0
  #   j = 2
  #   execute_on = 'initial timestep_end'
  # []
  # [gamma_12_v]
  #   type = MatVectorComponentAux
  #   variable = gamma_12_v
  #   property = gamma_ij
  #   i = 1
  #   j = 2
  #   execute_on = 'initial timestep_end'
  # []
  # [i_0]
  #   type = MatVectorComponentAux
  #   variable = i_0
  #   property = ij_i
  #   is_ij = true
  #   i = 0
  #   j = 1
  #   execute_on = 'initial timestep_end'
  # []
  # [j_0]
  #   type = MatVectorComponentAux
  #   variable = j_0
  #   property = ij_j
  #   is_ij = true
  #   i = 0
  #   j = 1
  #   execute_on = 'initial timestep_end'
  # []
  [theta_01_v]
    type = MatVectorComponentAux
    variable = theta_01_v
    property = theta_ij
    i = 0
    j = 1
    execute_on = 'initial timestep_end'
  []
  # [theta_02_v]
  #   type = MatVectorComponentAux
  #   variable = theta_02_v
  #   property = theta_ij
  #   i = 0
  #   j = 2
  #   execute_on = 'initial timestep_end'
  # []
  # [theta_12_v]
  #   type = MatVectorComponentAux
  #   variable = theta_12_v
  #   property = theta_ij
  #   i = 1
  #   j = 2
  #   execute_on = 'initial timestep_end'
  # []
  # [theta_11]
  #   type = MatVectorComponentAux
  #   variable = theta_11
  #   property = theta_ij
  #   i = 1
  #   j = 1
  #   execute_on = 'initial timestep_end'
  # []
  # [theta_01_x]
  #   type = MatVectorComponentAux
  #   variable = theta_01_x
  #   property = dtheta_dgradeta
  #   gradient = true
  #   component = 0
  #   i = 0
  #   j = 1
  #   execute_on = 'initial timestep_end'
  # []
  # [theta_01_y]
  #   type = MatVectorComponentAux
  #   variable = theta_01_y
  #   property = dtheta_dgradeta
  #   gradient = true
  #   component = 1
  #   i = 0
  #   j = 1
  #   execute_on = 'initial timestep_end'
  # []
  # [theta_01_z]
  #   type = MatVectorComponentAux
  #   variable = theta_01_z
  #   property = dtheta_dgradeta
  #   gradient = true
  #   component = 2
  #   i = 0
  #   j = 1
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
  [incl_vec]
    type = GGInclinationVector
    # grain_tracker = grain_tracker
    output_properties = 'inclination_vector ang_dist'
    gb_id_method = SWITCH
    ffc = gr_flood_uo
    # grain_tracker = gr_flood_uo
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
    gb_id_method = ffc
    ffc = gr_flood_uo
    angular_func = ATAN_2D
    intol = 0 #100 #10 #100 #200 # cut if alpha > intol
    altol = 0 #100 #10 #1.5 # cut if h*alpha > altol
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
    output_properties = 'gtnum no_ij_pairs gamma_a'
    # testout1 testout2 testoutgrad testouttens
    outputs = 'exodus'
  []
  [gamma_test]
    type = ElementalGammaMaterial
    gb_energy_iso_name = gbe_iso
    kappa = kappa_op
    free_energy_m = 5.521269e6
    well = false
    output_properties = 'gamma_asymm int_noij int_width'
    outputs = 'exodus'
  []
  # [aniso_mat]
  #   type = GGInclinationMaterial
  #   gb_case = ffc
  #   # grain_tracker = grain_tracker
  #   ffc = gr_flood_uo
  #   gb_energy_input = gbe_iso
  #   kappa = 2.07337e7 #0.3 #kappa
  #   free_energy_m = 5.521269e6 #4.5e6 #0.9375 #const_m
  #   L0 = L0
  #   gamma0 = gamma_iso
  #   #
  #   moelans_mu = true
  #   aniso_L = false
  #   delta_ij = 0.2
  #   theta_prefactor = 2
  #   inc_ij_0 = 0 #1.57
  #   continuous = false
  #   gt_tol = 0.00
  #   angular_func = ATAN_2D
  #   alphacase = BOTH
  #   intol = 100 #${i_tol} #200 # cut if alpha > intol
  #   altol = 100 #${a_tol} #10 #1.5 # cut if h*alpha > altol
  #   hgb = hgb
  #   # Output Names
  #   inclination_name = inclination_mat
  #   L_name = Lold
  #   gamma_name = gamma_old
  #   mu_name = muold
  #   gb_energy = gbeold
  #   output_properties = 'gamma_old Lold muold gbeold int_width_old' # t2tens' # inclination_mat t2tens'
  #   outputs = 'exodus'
  # []
  [hgb_a]
    type = ParsedMaterial
    property_name = hgb_a
    coupled_variables = 'gr0 gr1' # gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10' # gr11 gr12 gr13 gr14 gr15'
    expression = 'hb:=gr0*gr0 + gr1*gr1;
                  4 * (1 - hb) * (1 - hb)'
    # hg:=4 * (1 - hb) * (1 - hb);
    # if(hg>1.0,1.0,hg)'
    outputs = 'exodus'
    #  + gr11*gr11 + gr12*gr12 + gr13*gr13 + gr14*gr14 + gr15*gr15
  []
  [hgb_b]
    type = SwitchingFunctionGBMaterial
    h_name = hgb_b
    grain_ops = 'gr0 gr1' # gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10' # gr11 gr12 gr13 gr14 gr15'
    hgb_threshold = 0
    output_properties = 'hgb_b'
    outputs = 'exodus'
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
    outputs = 'exodus' #h3:=h1+h2;
  []
  [sumgr]
    type = ParsedMaterial
    property_name = sumgr
    coupled_variables = 'gr0 gr1' # gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10' # gr11 gr12 gr13 gr14 gr15'
    expression = 'gr0 + gr1' # + gr2 + gr3 + gr4 + gr5 + gr6 + gr7 + gr8 + gr9 + gr10'
    #  + gr11 + gr12 + gr13 + gr14 + gr15'
    outputs = 'exodus'
  []
  [sum_inc]
    type = ParsedMaterial
    property_name = sum_inc
    coupled_variables = 'gr0 gr1 inc_01'
    expression = 'num:= inc_01 * gr0 * gr0 * gr1 * gr1;
    den:= gr0 * gr0 * gr1 * gr1;
    num/den'
    outputs = 'exodus'
  []
  [sum_theta]
    type = ParsedMaterial
    property_name = sum_theta
    coupled_variables = 'gr0 gr1 theta_01_v'
    expression = 'num:= theta_01_v * gr0 * gr0 * gr1 * gr1;
    den:= gr0 * gr0 * gr1 * gr1;
    num/den'
    outputs = 'exodus'
  []
  [Lmat]
    type = ParsedMaterial
    property_name = Lmat
    material_property_names = 'L'
    expression = 'L'
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
  # [avg_grain_volumes]
  #   type = AverageGrainVolume
  #   feature_counter = grain_tracker
  #   execute_on = 'initial timestep_end'
  # []
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
  # num_steps = 2
  automatic_scaling = true
  compute_scaling_once = false
  # dt = 0.05
  # dtmax = 0.5
  # dt = 2e-5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.001
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
  exodus = true # Exodus file will be outputted
  #nemesis = true
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
  # [pgraph]
  #   type = PerfGraphOutput
  #   execute_on = 'initial final' # Default is "final"
  #   level = 2 # Default is 1
  #   heaviest_branch = true # Default is false
  #   heaviest_sections = 7 # Default is 0
  # []
  # file_base = 12_halt_hypre_nopre_i${i_tol}_a${a_tol}
  # file_base = 17_bicr_large_withnewKernel
  # file_base = test
[]
