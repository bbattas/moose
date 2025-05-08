##############################################################################
# File: 03_voronoi_early_anisoInverse.i
# File Location: /examples/agg/01_initial_testing/10_voronoi_early_large/03_voronoi_early_anisoInverse
# Created Date: Thursday May 8th 2025
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Thursday May 8th 2025
# Modified By: Brandon Battas
# -----
# Description:
#  Building a polycrystalline input so i can do aniso vs iso on a longer gg problem
#  Using a +pi/4 in the inclination function to invert the anisotropy
#
#
##############################################################################

[Mesh]
  # [gmg]
  #   type = DistributedRectilinearMeshGenerator
  #   dim = 2
  #   nx = 40
  #   ny = 40
  #   nz = 0
  #   xmin = 0
  #   xmax = 16 #1000
  #   ymin = 0
  #   ymax = 16 #1000
  #   zmin = 0
  #   zmax = 0
  #   # elem_type = QUAD4
  #   parallel_type = DISTRIBUTED
  #   uniform_refine = 0
  #   second_order = true
  # []
  type = GeneratedMesh
  dim = 2
  nx = 80
  ny = 80
  nz = 0
  xmin = 0
  xmax = 40 #1000
  ymin = 0
  ymax = 40 #1000
  zmin = 0
  zmax = 0
  # elem_type = QUAD4
  parallel_type = DISTRIBUTED
  uniform_refine = 0
  second_order = true
[]

[GlobalParams]
  op_num = 8
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
    grain_num = 16 # Number of grains
    rand_seed = 10
    int_width = 1.55
  []
  [grain_tracker]
    type = GrainTracker
    # variable = eta
    threshold = 0.2
    connecting_threshold = 0.08
    compute_halo_maps = true
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
  []
  [terminator]
    type = Terminator
    expression = 'grain_tracker < 3'
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
  # [ebsd_ic]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [ebsd_grains]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # # Halos
  # [halos]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [halo0]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [halo1]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [halo2]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
  # [halo3]
  #   order = CONSTANT
  #   family = MONOMIAL
  # []
[]

[Modules]
  [PhaseField]
    [GrainGrowth]
      mobility = L_aniso
      kappa = kappa
      order = SECOND
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
    coupled_variables = 'gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_0
    d2gamma_dgradop2_name = d2gammadgrad_eta2_0
    mob_name = L_aniso
    dL_dgradop_name = dLdgrad_eta_0
    d2L_dgradop2_name = d2Ldgrad_eta2_0
    variable_L = true
  []
  [gr1_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr1
    # v = 'gr0 gr2 gr3'
    coupled_variables = 'gr0 gr2 gr3 gr4 gr5 gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_1
    d2gamma_dgradop2_name = d2gammadgrad_eta2_1
    mob_name = L_aniso
    dL_dgradop_name = dLdgrad_eta_1
    d2L_dgradop2_name = d2Ldgrad_eta2_1
    variable_L = true
  []
  [gr2_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr2
    # v = 'gr0 gr1 gr3'
    coupled_variables = 'gr0 gr1 gr3 gr4 gr5 gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_2
    d2gamma_dgradop2_name = d2gammadgrad_eta2_2
    mob_name = L_aniso
    dL_dgradop_name = dLdgrad_eta_2
    d2L_dgradop2_name = d2Ldgrad_eta2_2
    variable_L = true
  []
  [gr3_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr3
    # v = 'gr0 gr1 gr2'
    coupled_variables = 'gr0 gr1 gr2 gr4 gr5 gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_3
    d2gamma_dgradop2_name = d2gammadgrad_eta2_3
    mob_name = L_aniso
    dL_dgradop_name = dLdgrad_eta_3
    d2L_dgradop2_name = d2Ldgrad_eta2_3
    variable_L = true
  []
  [gr4_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr4
    # v = 'gr0 gr1 gr2'
    coupled_variables = 'gr0 gr1 gr2 gr3 gr5 gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_4
    d2gamma_dgradop2_name = d2gammadgrad_eta2_4
    mob_name = L_aniso
    dL_dgradop_name = dLdgrad_eta_4
    d2L_dgradop2_name = d2Ldgrad_eta2_4
    variable_L = true
  []
  [gr5_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr5
    # v = 'gr0 gr1 gr2'
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr6 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_5
    d2gamma_dgradop2_name = d2gammadgrad_eta2_5
    mob_name = L_aniso
    dL_dgradop_name = dLdgrad_eta_5
    d2L_dgradop2_name = d2Ldgrad_eta2_5
    variable_L = true
  []
  [gr6_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr6
    # v = 'gr0 gr1 gr2'
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr7'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_6
    d2gamma_dgradop2_name = d2gammadgrad_eta2_6
    mob_name = L_aniso
    dL_dgradop_name = dLdgrad_eta_6
    d2L_dgradop2_name = d2Ldgrad_eta2_6
    variable_L = true
  []
  [gr7_ACIaniso]
    type = ACInterfaceAnisoGamma
    variable = gr7
    # v = 'gr0 gr1 gr2'
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6'
    gamma_name = gamma_asymm
    dgamma_dgradop_name = dgammadgrad_eta_7
    d2gamma_dgradop2_name = d2gammadgrad_eta2_7
    mob_name = L_aniso
    dL_dgradop_name = dLdgrad_eta_7
    d2L_dgradop2_name = d2Ldgrad_eta2_7
    variable_L = true
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
  # [ebsd_ic_aux]
  #   type = EBSDReaderPointDataAux
  #   variable = ebsd_ic
  #   ebsd_reader = ebsd_reader
  #   data_name = 'feature_id'
  #   execute_on = 'initial timestep_end'
  # []
  # [ebsd_grains]
  #   type = EBSDReaderAvgDataAux
  #   data_name = 'feature_id'
  #   ebsd_reader = ebsd_reader
  #   grain_tracker = grain_tracker
  #   variable = ebsd_grains
  #   execute_on = 'initial timestep_end'
  # []
  # [halos]
  #   type = FeatureFloodCountAux
  #   variable = halos
  #   flood_counter = grain_tracker
  #   field_display = HALOS
  #   execute_on = 'initial timestep_end'
  # []
  # HALOS
  # [halo0]
  #   type = FeatureFloodCountAux
  #   variable = halo0
  #   map_index = 0
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # []
  # [halo1]
  #   type = FeatureFloodCountAux
  #   variable = halo1
  #   map_index = 1
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # []
  # [halo2]
  #   type = FeatureFloodCountAux
  #   variable = halo2
  #   map_index = 2
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # []
  # [halo3]
  #   type = FeatureFloodCountAux
  #   variable = halo3
  #   map_index = 3
  #   field_display = HALOS
  #   flood_counter = grain_tracker
  # []
[]

# [BCs]
#   [Periodic]
#     [All]
#       auto_direction = 'x y'
#     []
#   []
# []

[Materials]
  # [Moly_GB]
  #   type = GBEvolution
  #   time_scale = 1.0
  #   GBmob0 = 3.986e-6
  #   T = 500 # K
  #   wGB = 60 # nm
  #   Q = 1.0307
  #   GBenergy = 2.4
  # []
  [constants]
    type = GenericConstantMaterial
    prop_names = 'L0 kappa const_m gamma_iso iw_iso gbe_iso'
    prop_values = '1.0 0.3 0.9375   1.5       1.55    0.25'
  []
  # [iso_constants]
  #   type = GenericConstantMaterial
  #   prop_names = 'L0 kappa const_m gamma_asymm gamma_iso iw_iso gbe_iso mu'
  #   prop_values = '1.0 0.3 0.9375   1.5           1.5     1.55   0.25  0.85'
  # []
  [incl_test01]
    type = GGInclinationMaterial
    grain_tracker = grain_tracker
    gb_energy_input = gbe_iso
    kappa = 0.3 #kappa
    free_energy_m = 0.9375 #const_m
    L0 = L0
    gamma0 = gamma_iso
    inc_ij_0 = 0.785
    # Output Names
    inclination_name = inclination_mat
    L_name = L_aniso
    gamma_name = gamma_asymm
    mu_name = mu
    output_properties = 'gamma_asymm L_aniso mu gb_energy int_width'
    outputs = 'exodus'
  []
  [sumgr]
    type = ParsedMaterial
    property_name = sumgr
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    expression = 'gr0 + gr1 + gr2 + gr3 + gr4 + gr5 + gr7 + gr7'
    outputs = 'exodus'
  []
  # Variable Gradients
  [ng_gr0]
    type = VariableGradientMaterial
    prop = ng_gr0
    variable = gr0
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
  [ng_gr6]
    type = VariableGradientMaterial
    prop = ng_gr6
    variable = gr6
  []
  [ng_gr7]
    type = VariableGradientMaterial
    prop = ng_gr7
    variable = gr7
  []
  # Free energies
  [fe_bulk_manual]
    type = ParsedMaterial
    property_name = fe_bulk_manual
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    material_property_names = 'gamma_asymm'
    expression = 'gmeta:= gr0^2 * gr1^2 + gr0^2 * gr2^2 + gr0^2 * gr3^2 + gr0^2 * gr4^2 + gr0^2 * gr5^2 + gr0^2 * gr6^2 + gr0^2 * gr7^2 +
     gr1^2 * gr2^2 + gr1^2 * gr3^2 + gr1^2 * gr4^2 + gr1^2 * gr5^2 + gr1^2 * gr6^2 + gr1^2 * gr7^2 +
     gr2^2 * gr3^2 + gr2^2 * gr4^2 + gr2^2 * gr5^2 + gr2^2 * gr6^2 + gr2^2 * gr7^2 +
     gr3^2 * gr4^2 + gr3^2 * gr5^2 + gr3^2 * gr6^2 + gr3^2 * gr7^2 +
     gr4^2 * gr5^2 + gr4^2 * gr6^2 + gr4^2 * gr7^2 + gr5^2 * gr6^2 + gr5^2 * gr7^2 + gr6^2 * gr7^2;
    etaover:= 0.25*(gr0^4 + gr1^4 + gr2^4 + gr3^4 + gr4^4 + gr5^4 + gr6^4 + gr7^4) - 0.5*(gr0^2 + gr1^2 + gr2^2 + gr3^2 + gr4^2 + gr5^2 + gr6^2 + gr7^2);
    etaover + gamma_asymm * gmeta + 0.25'
    outputs = 'exodus'
  []
  [fe_grad_manual]
    type = ParsedMaterial
    property_name = fe_grad_manual
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    material_property_names = 'kappa ng_gr0 ng_gr1 ng_gr2 ng_gr3 ng_gr4 ng_gr5 ng_gr6 ng_gr7'
    expression = '0.5 * kappa * (ng_gr0^2 + ng_gr1^2 + ng_gr2^2 + ng_gr3^2 + ng_gr4^2 + ng_gr5^2 + ng_gr6^2 + ng_gr7^2)'
    outputs = 'exodus'
  []
  [fe_tot_manual]
    type = ParsedMaterial
    property_name = fe_tot_manual
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    material_property_names = 'fe_bulk_manual fe_grad_manual'
    expression = 'fe_bulk_manual + fe_grad_manual'
    outputs = 'exodus'
  []
[]

[Postprocessors]
  # [ctr_grain_area]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr2
  #   execute_on = 'initial timestep_end'
  # []
  [tot_grain_area]
    type = ElementIntegralMaterialProperty
    mat_prop = sumgr
  []
  [gbe_total]
    type = ElementIntegralMaterialProperty
    mat_prop = gb_energy
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
  # Variable Residuals
  # [R_gr0]
  #   type = VariableResidual
  #   variable = gr0
  # []
  # [R_gr1]
  #   type = VariableResidual
  #   variable = gr1
  # []
  # [R_gr2]
  #   type = VariableResidual
  #   variable = gr2
  # []
  # [R_gr3]
  #   type = VariableResidual
  #   variable = gr3
  # []
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

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  nl_max_its = 30
  l_max_its = 60
  l_tol = 1e-05
  nl_rel_tol = 1e-8 #default is 1e-8
  # nl_abs_tol = 1e-14 #only needed when near equilibrium or veeeery small dt

  start_time = 0.0
  # end_time = 20
  # dt = 1
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    linear_iteration_ratio = 1e5
    dt = 1
  []
[]

[Outputs]
  exodus = true
  console = true
  csv = true
[]
