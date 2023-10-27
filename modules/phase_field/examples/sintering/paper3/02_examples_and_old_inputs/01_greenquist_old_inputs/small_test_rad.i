##############################################################################
# File: small_test_rad.i
# File Location: /examples/sintering/paper3/02_examples_and_old_inputs/01_greenquist_old_inputs
# Created Date: Thursday October 26th 2023
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Thursday October 26th 2023
# Modified By: Brandon Battas
# -----
# Description:
#  An old input from Ian's Marmot Fork (4 years old) so see some example of
#  the needed kernels/materials and whatnot for irradiation
#
#
##############################################################################

[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 10000
  nx = 125
  ymin = 0
  ymax = 10000
  ny = 125
  uniform_refine = 0
[]

[GlobalParams]
  op_num = 9
  var_name_base = gr
  int_width = 320
[]

[Functions]
  [f_T]
    type = ConstantFunction
    value = 1200
  []
  [f_rhomv_avg]
    type = ParsedFunction
    vars = 'rhomv_tot'
    vals = 'rhomv_tot'
    value = 'rhomv_tot * 1e-8'
  []
  [f_rhoi_avg]
    type = ParsedFunction
    vars = 'rhomi_tot'
    vals = 'rhomi_tot'
    value = 'rhomi_tot * 1e-8'
  []
  [f_phi]
    type = SolutionFunction
    solution = soln_function
    from_variable = phi
  []
  [f_gr0]
    type = SolutionFunction
    solution = soln_function
    from_variable = gr0
  []
  [f_gr1]
    type = SolutionFunction
    solution = soln_function
    from_variable = gr1
  []
  [f_gr2]
    type = SolutionFunction
    solution = soln_function
    from_variable = gr2
  []
  [f_gr3]
    type = SolutionFunction
    solution = soln_function
    from_variable = gr3
  []
  [f_gr4]
    type = SolutionFunction
    solution = soln_function
    from_variable = gr4
  []
  [f_gr5]
    type = SolutionFunction
    solution = soln_function
    from_variable = gr5
  []
  [f_gr6]
    type = SolutionFunction
    solution = soln_function
    from_variable = gr6
  []
  [f_gr7]
    type = SolutionFunction
    solution = soln_function
    from_variable = gr7
  []
  [f_gr8]
    type = SolutionFunction
    solution = soln_function
    from_variable = gr8
  []
[]

[Variables]
  [PolycrystalVariables]
    scaling = 10
  []
  [phi]
    scaling = 10
  []
  [wv]
  []
  [wi]
    scaling = 1e4
  []
[]

[Modules]
  [PhaseField]
    [GrandPotential]
      switching_function_names = 'hm hb'

      chemical_potentials = 'wv    wi'
      susceptibilities = 'chiv  chii'
      mobilities = 'chiDv chiDi'
      anisotropy = '1     0'
      free_energies_w = 'rhomv rhobv rhomi rhobi'

      kappa_gr = kappa
      gamma_gr = gamma
      free_energies_gr = 'omegam omegab'
      mobility_name_gr = Lm

      additional_ops = 'phi'
      kappa_op = kappa
      gamma_grxop = gamma
      free_energies_op = 'omegam omegab'
      mobility_name_op = Lb
    []
  []
[]

[Kernels]
  [barrier_phi]
    type = ACBarrierFunction
    variable = phi
    v = 'gr0 gr1 gr2 gr3'
    gamma = gamma
    mob_name = Lb
  []
  [kappa_phi]
    type = ACKappaFunction
    variable = phi
    mob_name = Lb
    kappa_name = kappa
  []
  #sources
  [source_wv]
    type = MaskedBodyForce
    variable = wv
    mask = source
  []
  [source_wi]
    type = MaskedBodyForce
    variable = wi
    mask = source
  []
  #sinks
  [recom_wv]
    type = GrandPotentialRecombination
    variable = wv
    rho = rhomv
    rho_r = rhomi
    D = Di
    omega = 0.04092
    a0 = 0.25
    Z = 250
    hm = 1
    args = 'wv wi phi gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8'
  []
  [recom_wi]
    type = GrandPotentialRecombination
    variable = wi
    rho = rhomi
    rho_r = rhomv
    D = Di
    omega = 0.04092
    a0 = 0.25
    Z = 250
    hm = 1
    args = 'wv wi phi gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8'
  []
  #damage
  [damage_wv]
    type = MaskedBodyForce
    variable = wv
    mask = rad_damage_v
  []
  [damage_wi]
    type = MaskedBodyForce
    variable = wi
    mask = rad_damage_i
  []
  #noise
  #  [./noise_wv]
  #    type = ConservedLangevinNoise
  #    variable = wv
  #    noise = noise
  #    amplitude = 1e-10
  #    multiplier = hm
  #  [../]
  #  [./noise_wi]
  #    type = ConservedLangevinNoise
  #    variable = wi
  #    noise = noise
  #    amplitude = 1e-10
  #    multiplier = hm
  #  [../]
[]

[AuxVariables]
  [bnds]
  []
  [T]
    order = CONSTANT
    family = MONOMIAL
  []
  [halos]
    order = CONSTANT
    family = MONOMIAL
  []
  [grns]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  [T_aux]
    type = FunctionAux
    variable = T
    function = f_T
    execute_on = 'timestep_begin'
  []
  [halos_aux]
    type = FeatureFloodCountAux
    variable = halos
    flood_counter = grain_tracker
    field_display = HALOS
  []
  [grns_aux]
    type = FeatureFloodCountAux
    variable = grns
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  []
[]

[ICs]
  [IC_phi]
    type = FunctionIC
    variable = phi
    function = f_phi
  []
  [IC_gr0]
    type = FunctionIC
    variable = gr0
    function = f_gr0
  []
  [IC_gr1]
    type = FunctionIC
    variable = gr1
    function = f_gr1
  []
  [IC_gr2]
    type = FunctionIC
    variable = gr2
    function = f_gr2
  []
  [IC_gr3]
    type = FunctionIC
    variable = gr3
    function = f_gr3
  []
  [IC_gr4]
    type = FunctionIC
    variable = gr4
    function = f_gr4
  []
  [IC_gr5]
    type = FunctionIC
    variable = gr5
    function = f_gr5
  []
  [IC_gr6]
    type = FunctionIC
    variable = gr6
    function = f_gr6
  []
  [IC_gr7]
    type = FunctionIC
    variable = gr7
    function = f_gr7
  []
  [IC_gr8]
    type = FunctionIC
    variable = gr8
    function = f_gr8
  []
  [IC_T]
    type = FunctionIC
    variable = T
    function = f_T
  []
[]

[Materials]
  [constants]
    type = GenericConstantMaterial
    prop_names = 'Va      kmv  kbv  kmi    kbi    ' #f_dot'
    prop_values = '0.04092 26.8 26.8 3.6e21 3.6e21 ' #1e-8'
  []
  [f_dot]
    type = RandomMaterial
    f_name = f_dot
    noise = noise
    outputs = exodus
    min_value = 0.8e-8
    max_value = 1.2e-8
  []
  [chiDv]
    type = GrandPotentialTensorMaterial
    f_name = chiDv
    solid_mobility = Lb
    void_mobility = Lm
    chi = chiv
    surface_energy = 19.7
    c = phi
    T = T
    D0 = 8.33e9
    GBmob0 = 1.4759e9
    Q = 2.77
    Em = 3.608
    bulkindex = 1
    gbindex = 1e6
    surfindex = 1e9
  []
  [chiDi]
    type = GrandPotentialScalarMaterial
    D_name = Di
    M_name = chiDi
    chemical_potential = wi
    void_variable = phi
    Temperature = T
    susceptibility = chii
    D0 = 1e13
    migration_energy = 2.0
    bulk_index = 1
    GB_index = 1e6
    surface_index = 1e9
  []
  [cmv_eq]
    type = GBSegregation
    c_name = cmv_eq
    formation_energy = -0.21343
    c_prefactor = 0.00011854
    GB_energy = 1.5
    Temperature = T
  []
  [cmi_eq]
    type = GBSegregation
    c_name = cmi_eq
    formation_energy = 6.9
    GB_energy = 1.5
    Temperature = T
  []
  [radiation_damage]
    type = GrandPotentialRadiationDamage
    chemical_potentials = 'wv wi'
    void_op = phi
    Temperature = T
    void_energy_coefficient = 'kbv kbi'
    solid_energy_coefficient = 'kmv kmi'
    solid_energy_model = PARABOLIC
    surface_energy = 19.7
    GB_energy = 9.86
    atomic_volume = 0.04092
    solid_equilibrium_concentrations = 'cmv_eq cmi_eq'
    void_equilibrium_concentrations = '1.0 0.0'
    interface_switch_value = 0.3
    switching_base = h
    potential_base = omega
    void_subscript = b
    solid_subscript = m
    potential_subscripts = 'v i'
    susceptibility_base = chi
    kappa_name = kappa
    gamma_name = gamma
    output_properties = 'hm hb rhomi rhomv rhobi rhobv'
    outputs = exodus
  []
  [source]
    type = DerivativeParsedMaterial
    f_name = source
    args = 'phi'
    material_property_names = 'hm(phi) f_dot'
    constant_names = 'X' #defects generated per cascade
    constant_expressions = '3'
    function = '2*hm*f_dot*X'
    derivative_order = 1
    outputs = exodus
  []
  #RADIATION DAMAGE
  [avg_functions]
    type = GenericFunctionMaterial
    prop_names = 'rhomv_avg rhomi_avg'
    prop_values = 'f_rhomv_avg f_rhomi_avg'
  []
  [hmrhomv]
    type = DerivativeParsedMaterial
    f_name = hmrhomv
    args = 'phi wv'
    material_property_names = 'hm(phi) rhomv(phi,wv)'
    function = 'hm*rhomv'
    derivative_order = 1
  []
  [hmrhomi]
    type = DerivativeParsedMaterial
    f_name = hmrhomi
    args = 'phi wi'
    material_property_names = 'hm(phi) rhomi(phi,wi)'
    function = 'hm*rhomi'
    derivative_order = 1
  []
  [rad_damage_v]
    type = DerivativeParsedMaterial
    f_name = rad_damage_v
    args = 'phi wv'
    material_property_names = 'hm(phi) rhomv(wv) rhomv_avg f_dot'
    constant_names = 'Vc' #volume of a collision cascade
    constant_expressions = '268'
    function = '-f_dot * 2 * Vc * (rhomv - rhomv_avg)'
    derivative_order = 1
    outputs = exodus
  []
  [rad_damage_i]
    type = DerivativeParsedMaterial
    f_name = rad_damage_i
    args = 'phi wi'
    material_property_names = 'hm(phi) rhomi(wi) rhomi_avg f_dot'
    constant_names = 'Vc' #volume of a collision cascade
    constant_expressions = '268'
    function = '-f_dot * 2 * Vc * (rhomi - rhomi_avg)'
    derivative_order = 1
    outputs = exodus
  []
  #OUTPUT VISUALIZATIONS
  [rhov]
    type = DerivativeParsedMaterial
    f_name = rhov
    args = 'wv phi'
    material_property_names = 'hm(phi) hb(phi) rhomv(wv) rhobv(wv)'
    function = 'hm * rhomv + hb * rhobv'
    outputs = exodus
    derivative_order = 1
  []
  [rhoi]
    type = DerivativeParsedMaterial
    f_name = rhoi
    args = 'wi phi'
    material_property_names = 'hm(phi) hb(phi) rhomi(wi) rhobi(wi)'
    function = 'hm * rhomi + hb * rhobi'
    outputs = exodus
    derivative_order = 1
  []
  [rhov_dev]
    type = ParsedMaterial
    f_name = rhov_dev
    material_property_names = 'rhov hm hb cmv_eq Va'
    constant_names = 'cbv_eq'
    constant_expressions = '1.0'
    function = 'rhov - (hm * cmv_eq + hb * cbv_eq) / Va'
    outputs = exodus
  []
  [rhoi_dev]
    type = ParsedMaterial
    f_name = rhoi_dev
    material_property_names = 'rhoi hm hb cmi_eq Va'
    constant_names = 'cbi_eq'
    constant_expressions = '0.0'
    function = 'rhoi - (hm * cmi_eq + hb * cbi_eq) / Va'
    outputs = exodus
  []
  [f_phase]
    type = ParsedMaterial
    f_name = f_phase
    args = 'wv wi'
    material_property_names = 'hm hb rhomv rhobv rhomi rhobi omegam omegab'
    function = 'hm * (omegam + rhomv * wv + rhomi * wi) + hb * (omegab + rhobv * wv + rhobi * wi)'
    outputs = exodus
  []
[]

[UserObjects]
  [grain_tracker]
    type = GrainTracker
    compute_halo_maps = true
    compute_var_to_feature_map = true
    remap_grains = true
    enable_var_coloring = true
    threshold = 0.2
    connecting_threshold = 0.1
    halo_level = 4
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [void_tracker]
    type = FeatureFloodCount
    variable = phi
    execute_on = 'INITIAL TIMESTEP_END'
    compute_var_to_feature_map = true
    threshold = 0.4
    connecting_threshold = 0.38
  []
  [terminator]
    type = Terminator
    expression = 'void_tracker = 1'
  []
  [soln_function]
    type = SolutionUserObject
    mesh = 'small_test_IC_out.e'
    timestep = LATEST
    execute_on = 'INITIAL'
  []
  [noise]
    type = ConservedUniformNoise
    execute_on = 'INITIAL TIMESTEP_BEGIN'
    allow_duplicate_execution_on_initial = true
  []
[]

[Postprocessors]
  [memory]
    type = MemoryUsage
    mem_type = physical_memory
  []
  [rhov_tot]
    type = ElementIntegralMaterialProperty
    mat_prop = rhov
    execute_on = TIMESTEP_BEGIN
  []
  [rhoi_tot]
    type = ElementIntegralMaterialProperty
    mat_prop = rhoi
    execute_on = TIMESTEP_BEGIN
  []
  [hm_tot]
    type = ElementIntegralMaterialProperty
    mat_prop = hm
    execute_on = TIMESTEP_BEGIN
  []
  [rhomv_tot]
    type = ElementIntegralMaterialProperty
    mat_prop = rhomv
    execute_on = TIMESTEP_BEGIN
  []
  [rhomi_tot]
    type = ElementIntegralMaterialProperty
    mat_prop = rhomi
    execute_on = TIMESTEP_BEGIN
  []
  [f_phase_tot]
    type = ElementIntegralMaterialProperty
    mat_prop = f_phase
    execute_on = TIMESTEP_END
  []
[]

[VectorPostprocessors]
  [voids]
    type = FeatureVolumeVectorPostprocessor
    output_centroids = true
    flood_counter = void_tracker
    execute_on = 'INITIAL TIMESTEP_END'
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
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -ksp_gmres_restart -sub_ksp_type -pc_factor_levels'
  petsc_options_value = ' asm      ilu          2               31                 preonly       4'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_package'
  #petsc_options_value = 'lu NONZERO superlu_dist'
  nl_max_its = 25
  l_max_its = 30
  l_tol = 1e-5
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-4
  start_time = 0
  end_time = 3e7
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 360
    optimal_iterations = 12
    iteration_window = 3
    growth_factor = 1.5
    cutback_factor = 0.5
  []
  line_search = none
  abort_on_solve_fail = false
[]

[Outputs]
  exodus = true
  csv = true
  checkpoint = true
  print_linear_residuals = false
[]

[Debug]
  show_var_residual_norms = true
[]
