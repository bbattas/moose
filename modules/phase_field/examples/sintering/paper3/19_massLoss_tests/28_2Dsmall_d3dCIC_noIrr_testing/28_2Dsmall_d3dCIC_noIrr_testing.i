##############################################################################
# File: 28_2Dsmall_d3dCIC_noIrr_testing.i
# File Location: /examples/sintering/paper3/19_massLoss_tests/28_2Dsmall_d3dCIC_noIrr_testing
# Created Date: Wednesday June 19th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Monday June 24th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Testing this new approach, trying various things to improve the solve
#
#
#
##############################################################################

f_dot = 1e-8
# # Thermal 5% Volume
# ks_vac = 1.302e4
# ks_int = 1.092e9
# # Irradiation 5% Volume
# ks_vac = 5.753e3
# ks_int = 1.116e7
# Thermal GB 1600K
ks_vac = 6.569e2
ks_int = 5.461e7
# # Irradiation GB 1600K
# ks_vac = 7.751e2
# ks_int = 5.711e5

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = ../00_d3d_txt/2D_20x20um_8umavg_cThEq.txt
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = ebsd_mesh
    combinatorial_geometry = 'x > 19000'
    block_id = 1
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 0
[]

[GlobalParams]
  profile = TANH
  int_width = 2000
  op_num = 3
  var_name_base = gr
[]

[Variables]
  [wvac]
    initial_condition = 0
  []
  [wint]
    initial_condition = 0
  []
  # [cvac_var]
  # []
  # [cint_var]
  # []
  [PolycrystalVariables]
  []
  [phi]
  []
[]

[AuxVariables]
  [bnds]
  []
  [T]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1600
  []
  [voids]
    order = CONSTANT
    family = MONOMIAL
  []
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
    []
  []
  [Void_external_IC]
    type = BoundingBoxIC
    variable = phi
    block = 1
    inside = 1
    outside = 0.0 #0.01
    x1 = 20000
    x2 = 35000
    y1 = -2000
    y2 = 22000
  []
  [Void_internal_IC]
    type = SmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0.0 #0.01
    radius = 5046.26504
    x1 = 8893.09397
    y1 = 8144.4245
    z1 = 0
    block = 0
  []
  # Concentrations from EBSD
  # [c_ic]
  #   type = EBSDConsVarIC
  #   variable = cvac_var
  #   ebsd_reader = ebsd_reader
  #   column = 0
  # []
  # [ci_ic]
  #   type = EBSDConsVarIC
  #   variable = cint_var
  #   ebsd_reader = ebsd_reader
  #   column = 1
  # []
[]

[BCs]
  [phi]
    type = NeumannBC
    variable = phi
    value = 0
    boundary = 'left right top bottom'
  []
  [gr0]
    type = NeumannBC
    variable = gr0
    value = 0
    boundary = 'left right top bottom'
  []
  [gr1]
    type = NeumannBC
    variable = gr1
    value = 0
    boundary = 'left right top bottom'
  []
  [gr2]
    type = NeumannBC
    variable = gr2
    value = 0
    boundary = 'left right top bottom'
  []
[]

[Modules]
  [PhaseField]
    [GrandPotential]
      switching_function_names = 'hv hs'
      # Chempot
      chemical_potentials = 'wvac wint'
      mobilities = 'chiuD chiiD' #cons_mob
      anisotropic = 'false false'
      susceptibilities = 'chiu chii'
      free_energies_w = 'rhovu rhosu rhovi rhosi'
      # Grains
      mobility_name_gr = L_mat
      kappa_gr = kappa
      gamma_gr = gamma
      free_energies_gr = 'omegav omegas'
      # Other OPs (Phi)
      additional_ops = 'phi'
      mobility_name_op = L_mat
      kappa_op = kappa
      gamma_grxop = gamma
      free_energies_op = 'omegav omegas' #empty when no phi'omegaa omegab'
      # Mass Conservation
      mass_conservation = false
      concentrations = 'cvac_var cint_var'
      hj_over_kVa = 'hoverk_vu hoverk_su hoverk_vi hoverk_si' #'hv_over_kVa hs_over_kVa' #
      hj_c_min = 'hv csueq_mask cvieq_mask csieq_mask' #cvueq_mask=hv*1 'hv_c_min hs_c_min' #
    []
  []
[]

# [Kernels]
#   # Irradiation
#   # Source/Generation
#   [source_vac]
#     type = MaskedBodyForce
#     variable = cvac_var #wvac
#     mask = rho_gen #change_vac #rho_gen
#     coupled_variables = 'phi wvac' # wint'
#   []
#   [source_int]
#     type = MaskedBodyForce
#     variable = cint_var #wint
#     mask = rho_gen #change_vac #rho_gen
#     coupled_variables = 'phi wint' # wvac'
#   []
#   # Sink/Recombination
#   [recombination_vac]
#     type = MatReaction
#     variable = cvac_var #wvac
#     mob_name = rho_recomb
#     args = 'wvac wint phi gr0 gr1' # gr0 gr1 gr2'
#     #coupled_variables = 'phi gr0 gr1 gr2'
#   []
#   [recombination_int]
#     type = MatReaction
#     variable = cint_var #wint
#     mob_name = rho_recomb
#     args = 'wvac wint phi gr0 gr1' # gr0 gr1 gr2'
#     #coupled_variables = 'wvac phi gr0 gr1 gr2'
#   []
#   # Damage/Mixing
#   [ballistic_mix_vac]
#     type = MatDiffusion
#     variable = cvac_var #wvac
#     diffusivity = rho_mixing_vac
#   []
#   [ballistic_mix_int]
#     type = MatDiffusion
#     variable = cint_var #wint
#     diffusivity = rho_mixing_int
#   []
# []

[AuxKernels]
  [bnds_aux]
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  [voids_aux]
    type = FeatureFloodCountAux
    variable = voids
    flood_counter = void_tracker
    field_display = UNIQUE_REGION
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  []
[]

[Materials]
  [consts]
    type = GenericConstantMaterial
    prop_names = 'Va cvieq_mask'
    prop_values = '0.04092 0.0'
  []
  [k_constants]
    type = GenericConstantMaterial
    prop_names = 'ksu kvu ksi kvi' # Using the GB based values (lowest of mine)
    prop_values = '${ks_vac} ${ks_vac} ${ks_int} ${ks_int}'
  []
  [gb_e_mat] # eV/nm^2
    type = ParsedMaterial
    property_name = gb_e_mat
    coupled_variables = 'T'
    constant_names = 'a b c'
    constant_expressions = '-5.87e-4 1.56 6.24151'
    expression = 'c*(a*T + b)'
    outputs = none
  []
  [hgb]
    type = DerivativeParsedMaterial
    property_name = hgb
    # derivative_order = 2
    coupled_variables = 'gr0 gr1 gr2'
    expression = '16 * ( (gr0 * gr1)^2 + (gr0 * gr2)^2  + (gr1 * gr2)^2 )'
    # expression = '4*(1 - (gr0^2 + gr1^2 + etad0^2))^2'
    outputs = none
  []
  [L_mat]
    type = DerivativeParsedMaterial
    property_name = L_mat
    material_property_names = 'L Lv'
    coupled_variables = 'phi'
    constant_names = 'p0'
    constant_expressions = '0.3'
    expression = 'hv:=if(phi<=0.0,0.0,if(phi>=p0,1.0,6*(phi/p0)^5 - 15*(phi/p0)^4 + 10*(phi/p0)^3));
                hv*Lv + (1-hv)*L'
    outputs = none #'nemesis'
  []
  [chiuD]
    type = GrandPotentialIsoMaterial
    f_name = chiuD
    solid_mobility = L #CHANGED FROM L
    void_mobility = Lv
    chi = chiu
    c = phi
    T = T
    D0 = 4.2488e11 #8.33e9
    Em = 4.23317 #3.608
    GBmob0 = 3.42828e10 # nm4/eVs #1.4759e9 # new value from Tonks/PC/Jake GG Paper
    Q = 3.01 #2.77 # new value from Tonks/PC/Jake GG Paper
    bulkindex = 1
    gbindex = -1 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = -1 #1e11
    GBwidth = 1.0
    surf_thickness = 1.0 #0.5
    iw_scaling = true
    D_out_name = vac_diffus
  []
  [chiiD]
    type = GrandPotentialIsoIntMaterial
    f_name = chiiD
    chi = chii
    c = phi
    T = T
    D0 = 4.0767e11 #1e13 #Ian's irradiation paper (Matzke 1987)
    Em = 4.08453089 #2 #Ian's irradiation paper (Matzke 1987)
    bulkindex = 1
    gbindex = -1 #10 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = -1 #100 #1e11
    GBwidth = 1.0
    surf_thickness = 1.0 #0.5
    iw_scaling = true
    D_out_name = int_diffus
  []
  # [cons_mob]
  #   type = DerivativeParsedMaterial
  #   property_name = cons_mob
  #   # derivative_order = 2
  #   material_property_names = 'chiD(w,phi,gr0,gr1) Va'
  #   expression = 'Va * chiD'
  # []
  [densificaiton]
    type = GrandPotentialMultiSinteringMaterial
    chemical_potential_vac = wvac
    chemical_potential_int = wint
    void_op = phi
    Temperature = T
    surface_energy = gb_e_mat #surf_e_mat #19.7
    grainboundary_energy = gb_e_mat #9.86
    vac_solid_energy_coefficient = ksu
    int_solid_energy_coefficient = ksi
    vac_void_energy_coefficient = kvu
    int_void_energy_coefficient = kvi
    equilibrium_vacancy_concentration = cv_eq
    equilibrium_interstitial_concentration = ci_eq
    solid_energy_model = PARABOLIC
    mass_conservation = false
  []
  [cv_eq]
    type = DerivativeParsedMaterial
    property_name = cv_eq
    # derivative_order = 2
    coupled_variables = 'gr0 gr1 gr2 phi wvac'
    material_property_names = 'hgb(phi,gr0,gr1,gr2)' #'rhovi(vac) rhosi hv(phi)'
    constant_names = 'cb cgb'
    # constant_expressions = '3.877e-04 4.347e-03' #Irradiation
    constant_expressions = '2.424e-06 5.130e-03' #No Irradiation- LANL
    expression = 'cgb * hgb + (1 - hgb)*cb'
    outputs = nemesis #'nemesis' # + phi^2
  []
  [ci_eq]
    type = DerivativeParsedMaterial
    property_name = ci_eq
    # derivative_order = 2
    coupled_variables = 'gr0 gr1 gr2 phi wint'
    material_property_names = 'hgb(phi,gr0,gr1,gr2)' # 'rhovi(wint) rhosi(wint) hv(phi)'
    constant_names = 'cb cgb'
    # constant_expressions = '7.258e-09 5.900e-06' #Irradiation
    constant_expressions = '1.667e-32 6.170e-08' #No Irradiation- LANL
    expression = 'cgb * hgb + (1 - hgb)*cb'
    outputs = nemesis #'nemesis' #+ phi^2
  []
  [cvac]
    type = ParsedMaterial
    property_name = cvac
    material_property_names = 'hs rhosu hv rhovu Va'
    expression = 'Va*(hs*rhosu + hv*rhovu)'
    outputs = nemesis
  []
  [cint]
    type = ParsedMaterial
    property_name = cint
    material_property_names = 'hs rhosi hv rhovi Va'
    expression = 'Va*(hs*rhosi + hv*rhovi)'
    outputs = nemesis
  []
  [dvac]
    type = ParsedMaterial
    property_name = dvac
    material_property_names = 'vac_diffus'
    expression = 'vac_diffus'
    outputs = nemesis
  []
  [dint]
    type = ParsedMaterial
    property_name = dint
    material_property_names = 'int_diffus'
    expression = 'int_diffus'
    outputs = nemesis
  []
  [hvoid]
    type = ParsedMaterial
    property_name = hvoid
    material_property_names = 'hv'
    expression = 'hv'
    outputs = nemesis
  []
  [hsolid]
    type = ParsedMaterial
    property_name = hsolid
    material_property_names = 'hs'
    expression = 'hs'
    outputs = nemesis
  []
  [rhosi_out]
    type = ParsedMaterial
    property_name = rhosi_out
    material_property_names = 'hs rhosi hv rhovi Va'
    expression = 'rhosi'
    outputs = nemesis
  []
  [rhovi_out]
    type = ParsedMaterial
    property_name = rhovi_out
    material_property_names = 'hs rhosi hv rhovi Va'
    expression = 'rhovi'
    outputs = nemesis
  []
  [rhosu_out]
    type = ParsedMaterial
    property_name = rhosu_out
    material_property_names = 'hs rhosu hv rhovu Va'
    expression = 'rhosu'
    outputs = nemesis
  []
  [rhovu_out]
    type = ParsedMaterial
    property_name = rhovu_out
    material_property_names = 'hs rhosu hv rhovu Va'
    expression = 'rhovu'
    outputs = nemesis
  []
  # CONSERVATION
  # h / (V*k)
  [hoverk_vu]
    type = DerivativeParsedMaterial
    property_name = hoverk_vu
    coupled_variables = 'phi gr0 gr1 gr2'
    derivative_order = 2
    material_property_names = 'hv(phi) Va kvu'
    expression = 'hv / (Va * kvu)'
  []
  [hoverk_su]
    type = DerivativeParsedMaterial
    property_name = hoverk_su
    coupled_variables = 'phi gr0 gr1 gr2'
    derivative_order = 2
    material_property_names = 'hs(phi) Va ksu'
    expression = 'hs / (Va * ksu)'
  []
  [hoverk_vi]
    type = DerivativeParsedMaterial
    property_name = hoverk_vi
    coupled_variables = 'phi gr0 gr1 gr2'
    derivative_order = 2
    material_property_names = 'hv(phi) Va kvi'
    expression = 'hv / (Va * kvi)'
  []
  [hoverk_si]
    type = DerivativeParsedMaterial
    property_name = hoverk_si
    coupled_variables = 'phi gr0 gr1 gr2'
    derivative_order = 2
    material_property_names = 'hs(phi) Va ksi'
    expression = 'hs / (Va * ksi)'
  []
  # h*ceq Masks
  # [cvueq_mask] # cvueq_mask = hv*1
  #   type = DerivativeParsedMaterial
  #   property_name = cvueq_mask
  #   coupled_variables = 'phi gr0 gr1 gr2'
  #   derivative_order = 2
  #   material_property_names = 'hv(phi)'
  #   expression = 'hv * 1'
  # []
  [csueq_mask]
    type = DerivativeParsedMaterial
    property_name = csueq_mask
    coupled_variables = 'phi gr0 gr1 gr2'
    derivative_order = 2
    material_property_names = 'hs(phi) cv_eq(phi,gr0,gr1,gr2)'
    expression = 'hs * cv_eq'
  []
  # [cvieq_mask]  # cvieq_mask = hv*0
  #   type = DerivativeParsedMaterial
  #   property_name = cvieq_mask
  #   coupled_variables = 'phi gr0 gr1'
  #   derivative_order = 2
  #   material_property_names = 'hv(phi)'
  #   expression = 'hv * 0.0'
  # []
  [csieq_mask]
    type = DerivativeParsedMaterial
    property_name = csieq_mask
    coupled_variables = 'phi gr0 gr1 gr2'
    derivative_order = 2
    material_property_names = 'hs(phi) ci_eq(phi,gr0,gr1,gr2)'
    expression = 'hs * ci_eq'
  []
  # IRRADIATION
  [combined_rho_vac]
    type = DerivativeParsedMaterial
    property_name = combined_rho_vac
    coupled_variables = 'wvac phi gr0 gr1 gr2'
    derivative_order = 2
    material_property_names = 'rhovu(wvac,gr0,gr1,gr2) rhosu(wvac,gr0,gr1,gr2) hv(phi) hs(phi)'
    expression = 'hv*rhovu + hs*rhosu'
    outputs = none #'nemesis'
  []
  [combined_rho_int]
    type = DerivativeParsedMaterial
    property_name = combined_rho_int
    coupled_variables = 'wint phi gr0 gr1 gr2'
    derivative_order = 2
    material_property_names = 'rhovi(wint,gr0,gr1,gr2) rhosi(wint,gr0,gr1,gr2) hv(phi) hs(phi)'
    expression = 'hv*rhovi + hs*rhosi' #'(1-hv)*rhos' #
    outputs = none #'nemesis'
  []
  [a_r]
    type = ParsedMaterial
    property_name = a_r
    coupled_variables = 'phi T gr0 gr1 gr2'
    constant_names = 'Va Z a_0 kB Di_0 Ei_B' # Di_0 Ei_B'
    constant_expressions = '0.04092 250 0.25 8.617343e-5 4.0767e11 4.08453089' #1e13 2
    material_property_names = 'hs(phi)'
    expression = 'dint:=Di_0 * exp(-Ei_B / (kB * T));
                  hs * Va * Z * dint / (a_0^2)'
    outputs = none
  []
  [rho_gen]
    type = DerivativeParsedMaterial
    property_name = rho_gen
    coupled_variables = 'phi'
    derivative_order = 1
    constant_names = 'Nc Nd noise f_dot'
    constant_expressions = '2 5 1 ${f_dot}'
    material_property_names = 'hs(phi) Va'
    expression = 'f_dot * noise * Nc * Nd * hs * Va'
    outputs = none
  []
  [rho_recomb]
    type = DerivativeParsedMaterial
    property_name = rho_recomb
    coupled_variables = 'wvac wint phi gr0 gr1 gr2'
    derivative_order = 2
    # additional_derivative_symbols = w # combined_rho_vac combined_rho_int
    material_property_names = 'Va a_r(phi,gr0,gr1,gr2)
      combined_rho_vac(wvac,phi,gr0,gr1,gr2) combined_rho_int(wint,phi,gr0,gr1,gr2)'
    expression = 'out:=a_r * combined_rho_vac * combined_rho_int * Va;
                  if(out>0.0,0.0-out,0.0)'
    outputs = none
  []
  [rho_mixing_vac]
    type = DerivativeParsedMaterial
    property_name = rho_mixing_vac
    coupled_variables = 'wvac phi'
    derivative_order = 2
    constant_names = 'Nc Vc noise tc Dc f_dot'
    constant_expressions = '2 268 1 1e-11 1e12 ${f_dot}'
    material_property_names = 'chiu(phi,wvac) Va'
    expression = 'f_dot * noise * Nc * tc * Vc * Dc * chiu * Va' # * hs
    outputs = none
  []
  [rho_mixing_int]
    type = DerivativeParsedMaterial
    property_name = rho_mixing_int
    coupled_variables = 'wint phi'
    derivative_order = 2
    constant_names = 'Nc Vc noise tc Dc f_dot'
    constant_expressions = '2 268 1 1e-11 1e12 ${f_dot}'
    material_property_names = 'chii(phi,wint) Va'
    expression = 'f_dot * noise * Nc * tc * Vc * Dc * chii * Va' # * hs
    outputs = none
  []
[]

[Postprocessors]
  # [memoryAll]
  #   type = MemoryUsage
  #   mem_units = megabytes
  #   outputs = csv
  #   execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
  #   report_peak_value = false
  # []
  # [memoryPeak]
  #   type = MemoryUsage
  #   mem_units = megabytes
  #   outputs = csv
  #   execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
  #   report_peak_value = true
  # []
  [memory1CPU]
    type = MemoryUsage
    mem_units = megabytes
    outputs = csv
    execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
    value_type = max_process
  []
  [n_DOFs]
    type = NumDOFs
    outputs = csv
  []
  # [cv_var_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = cvac_var
  # []
  # [ci_var_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = cint_var
  # []
  [cv_total]
    type = ElementIntegralMaterialProperty
    mat_prop = cvac
    # outputs = csv
  []
  [ci_total]
    type = ElementIntegralMaterialProperty
    mat_prop = cint
    # outputs = csv
  []
  [cv_avg]
    type = ElementAverageMaterialProperty
    mat_prop = cvac
    outputs = csv
  []
  [ci_avg]
    type = ElementAverageMaterialProperty
    mat_prop = cint
    outputs = csv
  []
  [hgb_total]
    type = ElementIntegralMaterialProperty
    mat_prop = hgb
    outputs = csv
  []
  [nonlinear]
    type = NumNonlinearIterations
    outputs = csv
  []
  [linear]
    type = NumLinearIterations
    outputs = csv
  []
  [residuals]
    type = NumResidualEvaluations
    outputs = csv
  []
  [runtime]
    type = PerfGraphData
    section_name = "Root"
    data_type = TOTAL
  []
  [timestep]
    type = TimestepSize
    outputs = csv
  []
  # FFC
  [void_tracker]
    type = FeatureFloodCount
    variable = phi
    threshold = 0.5 # 0.1 #0.2
    connecting_threshold = 0.5 #0.09 #0.08
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
    outputs = csv
  []
  [grain_tracker_fc]
    type = FeatureFloodCount
    variable = 'gr0 gr1 gr2' # gr3 gr4 gr5 gr6 gr7 gr8 gr9' #'gr0 gr1 gr2'
    threshold = 0.5 #0.2
    connecting_threshold = 0.5 #0.08
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
    outputs = csv
  []
  # Other irradiation based PPs
  [total_phi]
    type = ElementIntegralVariablePostprocessor
    variable = phi
    outputs = csv
  []
  [avg_rhov]
    type = ElementAverageMaterialProperty
    mat_prop = combined_rho_vac
    outputs = csv
  []
  [avg_rhoi]
    type = ElementAverageMaterialProperty
    mat_prop = combined_rho_int
    outputs = csv
  []
  [total_rhov]
    type = ElementIntegralMaterialProperty
    mat_prop = combined_rho_vac
    outputs = csv
  []
  [total_rhoi]
    type = ElementIntegralMaterialProperty
    mat_prop = combined_rho_int
    outputs = csv
  []
[]

[VectorPostprocessors]
  [voids]
    type = FeatureVolumeVectorPostprocessor
    flood_counter = void_tracker
    execute_on = 'initial timestep_end final'
    output_centroids = true #false #was true
    outputs = csv
  []
  [grain_sizes]
    type = FeatureVolumeVectorPostprocessor
    flood_counter = grain_tracker_fc
    execute_on = 'initial timestep_end final'
    output_centroids = false #was true
    outputs = csv
  []
  # [vectorMemory]
  #   type = VectorMemoryUsage
  #   mem_units = gigabytes
  #   outputs = csv
  #   execute_on = 'NONLINEAR LINEAR TIMESTEP_END'
  # []
[]

[UserObjects]
  [ebsd_reader]
    type = EBSDReader
    custom_columns = 2
  []
  [ebsd]
    type = PolycrystalEBSD
    # coloring_algorithm = bt
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    output_adjacency_matrix = true
    phase = 1
  []
  [grain_tracker]
    type = GrainTracker
    threshold = 0.1 #0.2
    connecting_threshold = 0.09 #0.08
    compute_halo_maps = true #true#false
    verbosity_level = 1
  []
  [terminator_2void]
    type = Terminator
    expression = 'void_tracker < 2'
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
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type'
  petsc_options_value = ' asm      lu           2                nonzero'
  nl_max_its = 30
  l_max_its = 60
  l_tol = 1e-06 #4
  nl_rel_tol = 1e-6 #6 #default is 1e-8
  # nl_abs_tol = 1e-14 #only needed when near equilibrium or veeeery small dt
  start_time = 0
  end_time = 1e5
  # num_steps = 50
  # steady_state_detection = true
  # # From tonks ode input
  automatic_scaling = true
  compute_scaling_once = false
  # line_search = none
  # dt = 0.001
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 0.001
  []
[]

[Outputs]
  csv = true
  exodus = false
  checkpoint = false
  print_linear_residuals = false
  [nemesis]
    type = Nemesis
  []
  file_base = NMC_kvol
[]

# [Debug]
#   show_var_residual_norms = true
# []
