##############################################################################
# File: 07_box_moreAniso_debug.i
# File Location: /examples/agg/01_initial_testing/16_complex_alpha/07_box_moreAniso_debug
# Created Date: Monday July 28th 2025
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Tuesday July 29th 2025
# Modified By: Brandon Battas
# -----
# Description:
#  box ic to look at interfaces in the two extremes of inclination
#
#
#
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 40
    ny = 40
    xmin = 0
    xmax = 60
    ymin = 0
    ymax = 60
    #  elem_type = QUAD4
    #  uniform_refine = 1
  []
  # type = GeneratedMesh
  # dim = 2
  # nx = 40 #80
  # ny = 40 #80
  # nz = 0
  # xmin = 0
  # xmax = 20 #40 #1000
  # ymin = 0
  # ymax = 20 #40 #1000
  # zmin = 0
  # zmax = 0
  # # elem_type = QUAD4
  # parallel_type = DISTRIBUTED
  uniform_refine = 1
  second_order = false
  parallel_type = DISTRIBUTED
[]

[GlobalParams]
  op_num = 6 #8
  var_name_base = 'gr'
[]

# [Variables]
#   [PolycrystalVariables]
#     order = SECOND
#     family = HERMITE
#   []
#   # [gr0]
#   #   order = SECOND
#   #   family = LAGRANGE
#   # []
#   # [gr1]
#   #   order = SECOND
#   #   family = LAGRANGE
#   # []
#   # [gr2]
#   #   order = SECOND
#   #   family = LAGRANGE
#   # []
#   # [gr3]
#   #   order = SECOND
#   #   family = LAGRANGE
#   # []
# []

[UserObjects]
  [voronoi]
    type = PolycrystalVoronoi
    # grain_num = 6#10 #16 # Number of grains
    # rand_seed = 22 #10
    int_width = 6 #1.55
    file_name = box_ctrs.txt
  []
  [grain_tracker]
    type = GrainTracker
    # variable = eta
    threshold = 0.001 #0.2
    connecting_threshold = 0.0008 #0.08
    compute_halo_maps = true
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
  []
  [terminator]
    type = Terminator
    expression = 'grain_tracker < 4'
    execute_on = TIMESTEP_END
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
  []
[]

# [BCs]
#   [Periodic]
#     [All]
#       auto_direction = 'x y'
#     []
#   []
# []

[AuxVariables]
  [bnds]
    order = FIRST
    family = LAGRANGE
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
  [ghost_regions]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Modules]
  [PhaseField]
    [GrainGrowth]
      mobility = L0 #_aniso
      kappa = kappa
      order = FIRST
      family = LAGRANGE
    []
  []
[]

[Kernels]
  # [PolycrystalKernel]
  # []
  [gr0_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr0
    # v = 'gr1 gr2 gr3'
    coupled_variables = 'gr1 gr2 gr3 gr4' # gr5' # gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_0
    d2gamma_dgradop2_name = d2gammadgrad_eta2_0
    mob_name = L0 #_aniso
    # dL_dgradop_name = dLdgrad_eta_0
    # d2L_dgradop2_name = d2Ldgrad_eta2_0
    variable_L = false
    skip_off = false
    mask_name = h_test
  []
  [gr1_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr1
    # v = 'gr0 gr2 gr3'
    coupled_variables = 'gr0 gr2 gr3 gr4 gr5' # gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_1
    d2gamma_dgradop2_name = d2gammadgrad_eta2_1
    mob_name = L0 #_aniso
    # dL_dgradop_name = dLdgrad_eta_1
    # d2L_dgradop2_name = d2Ldgrad_eta2_1
    variable_L = false
    skip_off = false
    mask_name = h_test
  []
  [gr2_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr2
    # v = 'gr0 gr1 gr3'
    coupled_variables = 'gr0 gr1 gr3 gr4 gr5' # gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_2
    d2gamma_dgradop2_name = d2gammadgrad_eta2_2
    mob_name = L0 #_aniso
    # dL_dgradop_name = dLdgrad_eta_2
    # d2L_dgradop2_name = d2Ldgrad_eta2_2
    variable_L = false
    skip_off = false
    mask_name = h_test
  []
  [gr3_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr3
    # v = 'gr0 gr1 gr2'
    coupled_variables = 'gr0 gr1 gr2 gr4 gr5' # gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_3
    d2gamma_dgradop2_name = d2gammadgrad_eta2_3
    mob_name = L0 #_aniso
    # dL_dgradop_name = dLdgrad_eta_3
    # d2L_dgradop2_name = d2Ldgrad_eta2_3
    variable_L = false
    skip_off = false
    mask_name = h_test
  []
  [gr4_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr4
    # v = 'gr0 gr1 gr2'
    coupled_variables = 'gr0 gr1 gr2 gr3 gr5' # gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_4
    d2gamma_dgradop2_name = d2gammadgrad_eta2_4
    mob_name = L0 #_aniso
    # dL_dgradop_name = dLdgrad_eta_4
    # d2L_dgradop2_name = d2Ldgrad_eta2_4
    variable_L = false
    skip_off = false
    mask_name = h_test
  []
  [gr5_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr5
    # v = 'gr0 gr1 gr2'
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4' # gr5' # gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_5
    d2gamma_dgradop2_name = d2gammadgrad_eta2_5
    mob_name = L0 #_aniso
    # dL_dgradop_name = dLdgrad_eta_5
    # d2L_dgradop2_name = d2Ldgrad_eta2_5
    variable_L = false
    skip_off = false
    mask_name = h_test
  []
[]

[AuxKernels]
  [BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = timestep_begin
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = timestep_begin
  []
  [ghosted_entities]
    type = FeatureFloodCountAux
    variable = ghost_regions
    flood_counter = grain_tracker
    field_display = GHOSTED_ENTITIES
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  # [constants]
  #   type = GenericConstantMaterial
  #   prop_names = 'L0 kappa const_m gamma_iso iw_iso gbe_iso'
  #   prop_values = '1.0 0.3 0.9375   1.5       1.55    0.25'
  # []
  # [constants]
  #   type = GenericConstantMaterial
  #   prop_names = 'L0 kappa const_m gamma_iso iw_iso gbe_iso'
  #   prop_values = '1.0 0.6 0.46875   1.5       1.55    0.25'
  # []
  # [iso_constants]
  #   type = GenericConstantMaterial
  #   prop_names = 'L0 kappa const_m gamma_asymm gamma_iso iw_iso gbe_iso mu'
  #   prop_values = '1.0 0.3 0.9375   1.5           1.5     1.55   0.25  0.85'
  # []
  [v_constants]
    type = GenericConstantMaterial
    prop_names = 'L0         gamma_iso kappa  iw_iso gbe_iso'
    prop_values = '1.15382e-6   1.5    2.07337e7   6   4.60748e6'
  []
  [incl_vec]
    type = GGInclinationVector
    # grain_tracker = grain_tracker
    output_properties = 'inclination_vector ang_dist'
    gb_id_method = SWITCH
    hgb = h_test
    hgb_threshold = 0.5
    outputs = exodus
  []
  [incl_test01]
    type = GGInclinationMaterial
    grain_tracker = grain_tracker
    gb_energy_input = gbe_iso
    kappa = 2.07337e7 #0.3 #kappa
    free_energy_m = 4.5e6 #0.46875 #0.9375 #const_m
    L0 = L0
    gamma0 = gamma_iso
    #
    delta_ij = 0.5
    theta_prefactor = 2
    inc_ij_0 = 0 #1.57
    continuous = false
    angular_func = atan
    alphacase = BOTH
    intol = 100#200 # cut if alpha > intol (-2)
    altol = 10#10 #1.5 # cut if h*alpha > altol (-1)
    hgb = h_test
    moelans_mu = false
    aniso_L = false
    # Output Names
    inclination_name = inclination_mat
    L_name = L_aniso
    gamma_name = gamma_asymm
    mu_name = mu
    gb_energy = sigma
    output_properties = 'gamma_asymm L_aniso mu sigma int_width alpha_out
    inclination_mat testout testout2 testout3 dgammadgrad_eta_0 d2gammadgrad_eta2_0 gtnum atens t2tens
    dgammadgrad_eta_1 dgammadgrad_eta_2 dgammadgrad_eta_3'
    outputs = 'exodus'
  []
  # [aniso_mat]
  #   type = GGInclinationMaterial
  #   grain_tracker = grain_tracker
  #   gb_energy_input = gbe_iso
  #   kappa = 2.07337e7 #0.3 #kappa
  #   free_energy_m = 4.5e6 #0.9375 #const_m
  #   L0 = L0
  #   gamma0 = gamma_iso
  #   #
  #   moelans_mu = false
  #   aniso_L = false
  #   theta_prefactor = 2
  #   inc_ij_0 = 0 #1.57
  #   continuous = false
  #   angular_func = atan
  #   alphacase = BOTH
  #   intol = 200 # cut if alpha > intol
  #   altol = 10 #1.5 # cut if h*alpha > altol
  #   hgb = hgb
  #   # Output Names
  #   inclination_name = inclination_mat
  #   L_name = L
  #   gamma_name = gamma_asymm
  #   mu_name = mu
  #   gb_energy = sigma
  #   output_properties = 'gamma_asymm L mu sigma int_width inclination_mat'
  #   outputs = 'exodus'
  # []
  [sumgr]
    type = ParsedMaterial
    property_name = sumgr
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5' # gr6 gr7'
    expression = 'gr0 + gr1 + gr2 + gr3 + gr4 + gr5' # + gr7 + gr7'
    outputs = 'exodus'
  []
  [hgb]
    type = ParsedMaterial
    property_name = hgb
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5'
    expression = 'hb:=gr0*gr0 + gr1*gr1 + gr2*gr2 + gr3*gr3 + gr4*gr4 + gr5*gr5;
                  hg:=4 * (1 - hb) * (1 - hb);
                  if(hg>1.0,hg,hg)'
    outputs = 'exodus'
  []
  [hgb2]
    type = SwitchingFunctionGBMaterial
    h_name = hgb2
    grain_ops = 'gr0 gr1 gr2 gr3 gr4 gr5'
    hgb_threshold = 0
    output_properties = 'hgb2'
    outputs = 'exodus'
  []
  [hgb_thr]
    type = ParsedMaterial
    property_name = hgb_thr
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'hgb'
    expression = 'if(hgb>1,1,hgb)'
    outputs = 'exodus'
  []
  [hgb2_thr]
    type = ParsedMaterial
    property_name = hgb2_thr
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'hgb2'
    expression = 'if(hgb2>1,1,hgb2)'
    outputs = 'exodus'
  []
  [h_test]
    type = ParsedMaterial
    property_name = h_test
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'hgb hgb2'
    expression = 'h1:=if(hgb>1,1,hgb);
                  h2:=if(hgb2>1,1,hgb2);
                  h3:=(hgb + hgb2)/2;
                  if(h3>1,1,h3)'
    outputs = 'exodus' #h3:=h1+h2;
  []
  [h_inc]
    type = ParsedMaterial
    property_name = h_inc
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'hgb h_test'
    expression = 'if(h_test>0.6,1,0)'
    outputs = 'exodus'
  []
  [h_alpha]
    type = ParsedMaterial
    property_name = h_alpha
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'hgb alpha_out h_test'
    expression = 'h_test * alpha_out'
    outputs = 'exodus'
  []
  # Variable Gradients
  [ng_gr0]
    type = VariableGradientMaterial
    prop = ng_gr0
    variable = gr0
    outputs = 'exodus'
  []
  [ng_gr1]
    type = VariableGradientMaterial
    prop = ng_gr1
    variable = gr1
  []
  [ng_gr2]
    type = VariableGradientMaterial
    prop = ng_gr2
    variable = gr2
  []
  [ng_gr3]
    type = VariableGradientMaterial
    prop = ng_gr3
    variable = gr3
  []
  [ng_gr4]
    type = VariableGradientMaterial
    prop = ng_gr4
    variable = gr4
  []
  [ng_gr5]
    type = VariableGradientMaterial
    prop = ng_gr5
    variable = gr5
  []
  # [ng_gr6]
  #   type = VariableGradientMaterial
  #   prop = ng_gr6
  #   variable = gr6
  # []
  # [ng_gr7]
  #   type = VariableGradientMaterial
  #   prop = ng_gr7
  #   variable = gr7
  # []
  # Free energies 6 OP
  [fe_bulk_manual]
    type = ParsedMaterial
    property_name = fe_bulk_manual
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'gamma_asymm'
    expression = 'gmeta:= gr0^2 * gr1^2 + gr0^2 * gr2^2 + gr0^2 * gr3^2 + gr0^2 * gr4^2 + gr0^2 * gr5^2 +
     gr1^2 * gr2^2 + gr1^2 * gr3^2 + gr1^2 * gr4^2 + gr1^2 * gr5^2 +
     gr2^2 * gr3^2 + gr2^2 * gr4^2 + gr2^2 * gr5^2 +
     gr3^2 * gr4^2 + gr3^2 * gr5^2 +
     gr4^2 * gr5^2;
    etaover:= 0.25*(gr0^4 + gr1^4 + gr2^4 + gr3^4 + gr4^4 + gr5^4) - 0.5*(gr0^2 + gr1^2 + gr2^2 + gr3^2 + gr4^2 + gr5^2);
    etaover + gamma_asymm * gmeta + 0.25'
    outputs = 'exodus'
  []
  [fe_grad_manual]
    type = ParsedMaterial
    property_name = fe_grad_manual
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'kappa ng_gr0 ng_gr1 ng_gr2 ng_gr3 ng_gr4 ng_gr5'
    expression = '0.5 * kappa * (ng_gr0^2 + ng_gr1^2 + ng_gr2^2 + ng_gr3^2 + ng_gr4^2 + ng_gr5^2)'
    outputs = 'exodus'
  []
  [fe_tot_manual]
    type = ParsedMaterial
    property_name = fe_tot_manual
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5'
    material_property_names = 'fe_bulk_manual fe_grad_manual'
    expression = 'fe_bulk_manual + fe_grad_manual'
    outputs = 'exodus'
  []
[]

[Postprocessors]
  [runtime]
    type = PerfGraphData
    section_name = "Root"
    data_type = TOTAL
  []
  [timestep]
    type = TimestepSize
  []
  [nl_its]
    type = NumNonlinearIterations
  []
  [l_its]
    type = NumLinearIterations
  []
  [max_mpi_memory]
    type = MemoryUsage
    value_type = max_process
    report_peak_value = True
    mem_units = MEGABYTES
    execute_on = 'NONLINEAR TIMESTEP_END'
  []
  [ctr_grain_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
    execute_on = 'initial timestep_end'
  []
  [tot_grain_area]
    type = ElementIntegralMaterialProperty
    mat_prop = sumgr
  []
  [gbe_total]
    type = ElementIntegralMaterialProperty
    mat_prop = sigma
  []
  [mu_total]
    type = ElementIntegralMaterialProperty
    mat_prop = mu
  []
  [fe_bulk]
    type = ElementIntegralMaterialProperty
    mat_prop = fe_bulk_manual
  []
  [fe_grad]
    type = ElementIntegralMaterialProperty
    mat_prop = fe_grad_manual
  []
  [fe_tot]
    type = ElementIntegralMaterialProperty
    mat_prop = fe_tot_manual
  []
  # Console value checks
  [avg_iw]
    type = ElementAverageMaterialProperty
    mat_prop = int_width
  []
  [avg_gamma]
    type = ElementAverageMaterialProperty
    mat_prop = gamma_asymm
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = 'PJFNK'

  # petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  # petsc_options_value = 'hypre boomeramg 31'
  petsc_options_iname = '-pc_type -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm        preonly       lu           2'

  nl_max_its = 12
  l_max_its = 60
  l_tol = 1e-05
  nl_rel_tol = 1e-8 #default is 1e-8
  # nl_abs_tol = 1e-14 #only needed when near equilibrium or veeeery small dt
  # nl_abs_tol = 1e-14

  start_time = 0.0
  end_time = 30
  # num_steps = 6

  # dt = 1
  automatic_scaling = true
  compute_scaling_once = false
  # line_search = none
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    linear_iteration_ratio = 1e5
    dt = 0.01
    # dt = 1
  []
[]

[Outputs]
  exodus = true
  console = true
  csv = true
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'initial final' # Default is "final"
    level = 2 # Default is 1
    heaviest_branch = true # Default is false
    heaviest_sections = 7 # Default is 0
  []
  [dbg]
    type = VariableResidualNormsDebugOutput
  []
[]

# [Debug]
#   [dbg]
#     type = VariableResidualNormsDebugOutput
#   []
# []
