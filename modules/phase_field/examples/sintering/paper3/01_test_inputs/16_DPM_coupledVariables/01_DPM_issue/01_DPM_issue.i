##############################################################################
# File: 01_DPM_issue.i
# File Location: /examples/sintering/paper3/01_test_inputs/16_DPM_coupledVariables/01_DPM_issue
# Created Date: Monday August 26th 2024
# Author: Battas,Brandon Scott (bbattas@ufl.edu)
# -----
# Last Modified: Tuesday August 27th 2024
# Modified By: Battas,Brandon Scott
# -----
# Description:
#  A version of the input that hangs based on hgb
#
#
#
##############################################################################



[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = ../../../21_multiApps_IC/05_2D_irr_tests/00_sub/2D_40x40um_8umavg_8pore.txt
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = ebsd_mesh
    combinatorial_geometry = 'x > 38000'
    block_id = 1
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 0
[]

[GlobalParams]
  # profile = TANH
  int_width = 2000
  op_num = 8
  var_name_base = gr
[]

[Variables]
  [wvac]
    initial_condition = 0
  []
  [wint]
    initial_condition = 0
  []
  [PolycrystalVariables]
  []
  [phi]
  []
[]

[AuxVariables]
  [T]
    order = CONSTANT
    family = MONOMIAL
    initial_condition = 1600
  []
  # Vacancy mat to variable
  [cvac_aux]
    family = LAGRANGE
    order = FIRST
  []
  [cvac_aux_elem]
    family = MONOMIAL
    order = CONSTANT
  []
  # Interstitial mat to variable
  [cint_aux]
    family = LAGRANGE
    order = FIRST
  []
  [cint_aux_elem]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  # Vacancy
  [cvac_aux_elem]
    type = MaterialRealAux
    variable = 'cvac_aux_elem'
    property = 'cvac'
  []
  [cvac_aux]
    type = ProjectionAux
    variable = cvac_aux
    v = cvac_aux_elem
  []
  # Interstitial
  [cint_aux_elem]
    type = MaterialRealAux
    variable = 'cint_aux_elem'
    property = 'cint'
  []
  [cint_aux]
    type = ProjectionAux
    variable = cint_aux
    v = cint_aux_elem
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = ebsd
    []
  []
  [VoidIC] #External
    type = BoundingBoxIC
    variable = phi
    block = 1
    inside = 1
    outside = 0.0 #0.01
    x1 = 40000
    x2 = 80000
    y1 = -5000
    y2 = 55000
  []
  [voidIC2] #Internal
    type = SpecifiedSmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0 #0.01
    radii = '2157.03 1863.05 2081.   2036.04 1826.13 1772.42 1798.55 2060.82'
    x_positions = '30095.6  29046.79 20283.12  9370.67 19428.33 11649.88 24936.2   7106.3'
    y_positions = '7653.64 18136.86 20731.16 27581.73  8101.88 11721.34 28519.04 18280.74'
    z_positions = '    0.       0.       0.       0.       0.       0.       0.       0.  '
    # radii = '2157.03 1863.05 2081.  2036.04 1826.13 1772.42 1798.55 2060.82'
    # x_positions = '30095.6  29046.79 20283.12  9370.67 19428.33 11649.88 24936.2   7106.3'
    # y_positions = '7653.64 18136.86 20731.16 27581.73  8101.88 11721.34 28519.04 18280.74'
    # z_positions = '    0.       0.       0.       0.       0.      0.        0.       0.'
    block = 0
  []
[]

[BCs]
  [phi_bc]
    type = NeumannBC
    variable = 'phi'
    boundary = 'left right top bottom'
    value = 0
  []
  [gr0_bc]
    type = NeumannBC
    variable = 'gr0'
    boundary = 'left right top bottom'
    value = 0
  []
  [gr1_bc]
    type = NeumannBC
    variable = 'gr1'
    boundary = 'left right top bottom'
    value = 0
  []
  [gr2_bc]
    type = NeumannBC
    variable = 'gr2'
    boundary = 'left right top bottom'
    value = 0
  []
  [gr3_bc]
    type = NeumannBC
    variable = 'gr3'
    boundary = 'left right top bottom'
    value = 0
  []
  [gr4_bc]
    type = NeumannBC
    variable = 'gr4'
    boundary = 'left right top bottom'
    value = 0
  []
  [gr5_bc]
    type = NeumannBC
    variable = 'gr5'
    boundary = 'left right top bottom'
    value = 0
  []
  [gr6_bc]
    type = NeumannBC
    variable = 'gr6'
    boundary = 'left right top bottom'
    value = 0
  []
  [gr7_bc]
    type = NeumannBC
    variable = 'gr7'
    boundary = 'left right top bottom'
    value = 0
  []
  [wv_bc]
    type = NeumannBC
    variable = 'wvac'
    boundary = 'left right top bottom'
    value = 0
  []
  [wi_bc]
    type = NeumannBC
    variable = 'wint'
    boundary = 'left right top bottom'
    value = 0
  []
[]

[Modules]
  [PhaseField]
    [GrandPotentialAlt]
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
      # concentrations = 'cvac_var cint_var'
      # hj_over_kVa = 'hoverk_vu hoverk_su hoverk_vi hoverk_si' #'hv_over_kVa hs_over_kVa' #
      # hj_c_min = 'cvueq_mask csueq_mask cvieq_mask csieq_mask' #cvueq_mask=hv*1 'hv_c_min hs_c_min' #
    []
  []
[]

# [Kernels]
#   # Sintering Terms for gb/surface energy
#   [barrier_phi]
#     type = ACBarrierFunction
#     variable = phi
#     v = 'gr0 gr1 gr2' #gr2 gr3'# gr4 gr5'# gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13 gr14 gr15'
#     gamma = gamma
#     mob_name = L_mat
#   []
#   [kappa_phi]
#     type = ACKappaFunction
#     variable = phi
#     mob_name = L_mat
#     kappa_name = kappa
#   []
# []

[Materials]
  [consts]
    type = GenericConstantMaterial
    prop_names = 'Va negOverVa cvieq_mask'
    prop_values = '0.04092 -24.4379 0.0'
  []
  [k_constants]
    type = GenericConstantMaterial
    prop_names = 'ksu kvu ksi kvi' # Using the GB based values (lowest of mine)
    # prop_values = '7.751e2 7.751e2 5.711e5 5.711e5' # Irradiation
    prop_values = '6.569e2 6.569e2  5.461e7 5.461e7' # No Irradiation
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
  # [surf_e_mat] # eV/nm^2
  #   type = ParsedMaterial
  #   property_name = surf_e_mat
  #   material_property_names = gb_e_mat
  #   expression = '2 * gb_e_mat'
  #   outputs = none
  # []
  [hgb_out]
    type = SwitchingFunctionGBMaterial
    h_name = hgb_out
    grain_ops = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
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
    f_name = chiuD #chiuD
    solid_mobility = L #CHANGED FROM L
    void_mobility = Lv
    chi = chiu
    c = phi
    T = T
    D0 = 4.2488e11 #8.33e9
    Em = 4.23317 #3.608
    GBmob0 = 3.42828e10 # nm4/eVs #1.4759e9 # new value from Tonks/PC/Jake GG Paper
    Q = 3.01 #2.77 # new value from Tonks/PC/Jake GG Paper
    vaporindex = 1
    bulkindex = 1
    gbindex = -1 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = -1 #1e11
    GBwidth = 3.5 # based on avg of two lanl values
    surf_thickness = 3.5 # keeping equal to gb for simplicity
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
    vaporindex = 1
    bulkindex = 1
    gbindex = -1 #10 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = -1 #100 #1e11
    GBwidth = 3.5 # based on avg of two lanl values
    surf_thickness = 3.5 # keeping equal to gb for simplicity
    iw_scaling = true
    D_out_name = int_diffus
  []
  [densificaiton]
    type = GrandPotentialMultiSinteringMaterial
    chemical_potential_vac = wvac
    chemical_potential_int = wint
    void_op = phi
    Temperature = T
    surface_energy = gb_e_mat #gb_e_mat #surf_e_mat #19.7
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
    derivative_order = 2
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 phi wvac' #gr6 gr7
    material_property_names = 'hgb_out(gr0,gr1,gr2,gr3,gr4,gr5,gr6,gr7)'#,gr6,gr7 #'rhovi(vac) rhosi hv(phi)'
    constant_names = 'cb cgb'
    # constant_expressions = '3.877e-04 4.347e-03' #Irradiation
    constant_expressions = '2.424e-06 5.130e-03' #No Irradiation- LANL
    # expression = 'hg:=16 * ( (gr0 * gr1)^2 + (gr0 * gr2)^2 + (gr0 * gr3)^2 + (gr0 * gr4)^2 +
    #                       (gr0 * gr5)^2 + (gr0 * gr6)^2 + (gr0 * gr7)^2 +
    #                       (gr1 * gr2)^2 + (gr1 * gr3)^2 + (gr1 * gr4)^2 +
    #                       (gr1 * gr5)^2 + (gr1 * gr6)^2 + (gr1 * gr7)^2 +
    #                       (gr2 * gr3)^2 + (gr2 * gr4)^2 + (gr2 * gr5)^2 + (gr2 * gr6)^2 + (gr2 * gr7)^2 +
    #                       (gr3 * gr4)^2 + (gr3 * gr5)^2 + (gr3 * gr6)^2 + (gr3 * gr7)^2 +
    #                       (gr4 * gr5)^2 + (gr4 * gr6)^2 + (gr4 * gr7)^2 +
    #                       (gr5 * gr6)^2 + (gr5 * gr7)^2 + (gr6 * gr7)^2 );
    #               hgb:=if(hg>1e-4,hg,0.0);
    #               cgb * hgb + (1 - hgb)*cb'
    expression = 'cgb * hgb_out + (1 - hgb_out)*cb'
    outputs = none #'nemesis' # + phi^2
  []
  [ci_eq]
    type = DerivativeParsedMaterial
    property_name = ci_eq
    derivative_order = 2
    coupled_variables = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 phi wint' #gr6 gr7
    # material_property_names = 'hgb'#(phi,gr0,gr1,gr2,gr3,gr4,gr5)'#,gr6,gr7 # 'rhovi(wint) rhosi(wint) hv(phi)'
    constant_names = 'cb cgb'
    # constant_expressions = '7.258e-09 5.900e-06' #Irradiation
    constant_expressions = '1.667e-32 6.170e-08' #'1.667e-32 6.170e-08' #No Irradiation- LANL
    # constant_expressions = '2.424e-06 5.130e-03' #No Irradiation VACANCY- LANL
    expression = 'hg:=16 * ( (gr0 * gr1)^2 + (gr0 * gr2)^2 + (gr0 * gr3)^2 + (gr0 * gr4)^2 +
                          (gr0 * gr5)^2 + (gr0 * gr6)^2 + (gr0 * gr7)^2 +
                          (gr1 * gr2)^2 + (gr1 * gr3)^2 + (gr1 * gr4)^2 +
                          (gr1 * gr5)^2 + (gr1 * gr6)^2 + (gr1 * gr7)^2 +
                          (gr2 * gr3)^2 + (gr2 * gr4)^2 + (gr2 * gr5)^2 + (gr2 * gr6)^2 + (gr2 * gr7)^2 +
                          (gr3 * gr4)^2 + (gr3 * gr5)^2 + (gr3 * gr6)^2 + (gr3 * gr7)^2 +
                          (gr4 * gr5)^2 + (gr4 * gr6)^2 + (gr4 * gr7)^2 +
                          (gr5 * gr6)^2 + (gr5 * gr7)^2 + (gr6 * gr7)^2 );
                  hgb:=if(hg>1e-4,hg,0.0);
                  cgb * hgb + (1 - hgb)*cb'
    outputs = none #'nemesis' #+ phi^2
  []
  [cvac]
    type = ParsedMaterial
    property_name = cvac
    material_property_names = 'hs rhosu hv rhovu Va'
    expression = 'Va*(hs*rhosu + hv*rhovu)'
    outputs = none
  []
  [cint]
    type = ParsedMaterial
    property_name = cint
    material_property_names = 'hs rhosi hv rhovi Va'
    expression = 'Va*(hs*rhosi + hv*rhovi)'
    outputs = none
  []
[]

[Postprocessors]
  # [cv_var_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = cvac_var
  # []
  [wvac_total]
    type = ElementIntegralVariablePostprocessor
    variable = wvac
  []
  # [ci_var_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = cint_var
  # []
  [wint_total]
    type = ElementIntegralVariablePostprocessor
    variable = wint
  []
  [phi_total]
    type = ElementIntegralVariablePostprocessor
    variable = phi
  []
  # [gr0_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr0
  # []
  # [gr1_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr1
  # []
  # [gr2_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = gr2
  # []
  [hgb_total]
    type = ElementIntegralMaterialProperty
    mat_prop = hgb_out
  []
  [cvac_total]
    type = ElementIntegralMaterialProperty
    mat_prop = cvac
  []
  [cint_total]
    type = ElementIntegralMaterialProperty
    mat_prop = cint
  []
  [runtime]
    type = PerfGraphData
    section_name = "Root"
    data_type = TOTAL
  []
  # [void_tracker]
  #   type = FeatureFloodCount
  #   variable = phi
  #   threshold = 0.5 # 0.1 #0.2
  #   connecting_threshold = 0.5 #0.09 #0.08
  #   compute_var_to_feature_map = true
  #   execute_on = 'initial timestep_end'
  # []
[]

# [VectorPostprocessors]
#   [voids]
#     type = FeatureVolumeVectorPostprocessor
#     flood_counter = void_tracker
#     execute_on = 'initial timestep_end final'
#     output_centroids = true #false #was true
#     outputs = csv
#   []
# []

[UserObjects]
  [ebsd_reader]
    type = EBSDReader
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
    compute_halo_maps = false #true#false
    verbosity_level = 1
  []
  #   [terminator]
  #     type = Terminator
  #     expression = 'void_tracker < 2'
  #     execute_on = TIMESTEP_END
  #   []
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
  # start_time = 0
  # end_time = 10000
  num_steps = 1
  # steady_state_detection = true
  # # From tonks ode input
  automatic_scaling = true
  compute_scaling_once = false
  # line_search = none
  # dt = 1.0
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 0.01
  []
[]

[Outputs]
  csv = false
  exodus = false
  nemesis = false
  checkpoint = false
  # file_base = fullList
[]

# [Debug]
#   show_var_residual_norms = true
# []
