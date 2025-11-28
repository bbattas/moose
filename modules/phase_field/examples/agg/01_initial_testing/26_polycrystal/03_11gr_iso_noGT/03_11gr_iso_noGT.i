##############################################################################
# File: 03_11gr_iso_noGT.i
# File Location: /examples/agg/01_initial_testing/26_polycrystal/03_11gr_iso_noGT
# Created Date: Friday November 28th 2025
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Friday November 28th 2025
# Modified By: Brandon Battas
# -----
# Description:
#  11 grain voronoi polycrystal ic, without using grain tracker
#
#
#
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2 # Problem dimension
    nx = 80 # Number of elements in the x-direction
    ny = 80 # Number of elements in the y-direction
    xmin = 0 # minimum x-coordinate of the mesh
    xmax = 80 # maximum x-coordinate of the mesh
    ymin = 0 # minimum y-coordinate of the mesh
    ymax = 80 # maximum y-coordinate of the mesh
    # elem_type = QUAD4 # Type of elements used in the mesh
    # uniform_refine = 3 # Initial uniform refinement of the mesh
  []
  parallel_type = DISTRIBUTED # Periodic BCs
[]

[GlobalParams]
  op_num = 11
  var_name_base = gr
[]

[Variables]
  [PolycrystalVariables]
  []
[]

[UserObjects]
  [voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = bt
    file_name = '../00_sub/small_11gr.txt'
    int_width = 6
  []
  [gr_flood_uo]
    type = FeatureFloodCount
    variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10'
    threshold = 0.01
    connecting_threshold = 0.01
    compute_var_to_feature_map = true
    compute_halo_maps = true # Only necessary for displaying HALOS
    execute_on = 'initial timestep_end'
    # use_less_than_threshold_comparison = true
  []
  [term]
    type = Terminator
    expression = 'gr_flood_uo < 5'
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
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
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
    variable_mobility = false
    # order = SECOND
  []
  # [gr0_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr0
  #   coupled_variables = 'gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10'
  #   debug_kernel = false
  #   skip_off = false
  #   variable_L = false
  #   mask_name = hgb
  # []
  # [gr1_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr1
  #   coupled_variables = 'gr0 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10'
  #   debug_kernel = false
  #   skip_off = false
  #   variable_L = false
  #   mask_name = hgb
  # []
  # [gr2_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr2
  #   coupled_variables = 'gr0 gr1 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10'
  #   debug_kernel = false
  #   skip_off = false
  #   variable_L = false
  #   mask_name = hgb
  # []
  # [gr3_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr3
  #   coupled_variables = 'gr0 gr1 gr2 gr4 gr5 gr6 gr7 gr8 gr9 gr10'
  #   debug_kernel = false
  #   skip_off = false
  #   variable_L = false
  #   mask_name = hgb
  # []
  # [gr4_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr4
  #   coupled_variables = 'gr0 gr1 gr2 gr3 gr5 gr6 gr7 gr8 gr9 gr10'
  #   debug_kernel = false
  #   skip_off = false
  #   variable_L = false
  #   mask_name = hgb
  # []
  # [gr5_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr5
  #   coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr6 gr7 gr8 gr9 gr10'
  #   debug_kernel = false
  #   skip_off = false
  #   variable_L = false
  #   mask_name = hgb
  # []
  # [gr6_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr6
  #   coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr7 gr8 gr9 gr10'
  #   debug_kernel = false
  #   skip_off = false
  #   variable_L = false
  #   mask_name = hgb
  # []
  # [gr7_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr7
  #   coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr8 gr9 gr10'
  #   debug_kernel = false
  #   skip_off = false
  #   variable_L = false
  #   mask_name = hgb
  # []
  # [gr8_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr8
  #   coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr9 gr10'
  #   debug_kernel = false
  #   skip_off = false
  #   variable_L = false
  #   mask_name = hgb
  # []
  # [gr9_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr9
  #   coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr10'
  #   debug_kernel = false
  #   skip_off = false
  #   variable_L = false
  #   mask_name = hgb
  # []
  # [gr10_ACInc]
  #   type = ACInterfaceInclinationGamma
  #   variable = gr10
  #   coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9'
  #   debug_kernel = false
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
  [GBInc_cos]
    type = GBInclination
    gb_id_method = ffc
    ffc = gr_flood_uo
    angular_func = ATAN_2D
    intol = 0 #100 #10 #100 #200 # cut if alpha > intol
    altol = 0 #100 #10 #1.5 # cut if h*alpha > altol
    # Inclination function
    inc_func = COS
    ifunc_a = 0.0
    ifunc_b = 2
    ifunc_c = 0
    ifunc_d = 0
    # L
    aniso_L = false
    # Other Properties
    gb_energy_iso_name = gbe_iso
    kappa = kappa_op
    free_energy_m = 5.521269e6
    output_properties = 'gtnum no_ij_pairs'
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
  [hgb_a]
    type = ParsedMaterial
    property_name = hgb_a
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10'
    expression = 'hb:=gr0*gr0 + gr1*gr1 + gr2*gr2 + gr3*gr3 + gr4*gr4 + gr5*gr5 + gr6*gr6 + gr7*gr7 + gr8*gr8 + gr9*gr9 + gr10*gr10;
                  4 * (1 - hb) * (1 - hb)'
    outputs = 'exodus'
  []
  [hgb_b]
    type = SwitchingFunctionGBMaterial
    h_name = hgb_b
    grain_ops = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10'
    hgb_threshold = 0
    output_properties = 'hgb_b'
    outputs = 'exodus'
  []
  [hgb]
    type = ParsedMaterial
    property_name = hgb
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10'
    material_property_names = 'hgb_a hgb_b'
    expression = 'h3:=(hgb_a + hgb_b)/2;
                  if(h3>1,1,if(h3<0,0.0,h3))'
    outputs = 'exodus'
  []
  [sumgr]
    type = ParsedMaterial
    property_name = sumgr
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10'
    expression = 'gr0 + gr1 + gr2 + gr3 + gr4 + gr5 + gr6 + gr7 + gr8 + gr9 + gr10'
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

  automatic_scaling = true
  compute_scaling_once = false

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.01
    cutback_factor = 0.9
    growth_factor = 1.1
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
  checkpoint = false
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final' # Default is "final"
    level = 2 # Default is 1
    heaviest_branch = true # Default is false
    heaviest_sections = 7 # Default is 0
  []
[]
