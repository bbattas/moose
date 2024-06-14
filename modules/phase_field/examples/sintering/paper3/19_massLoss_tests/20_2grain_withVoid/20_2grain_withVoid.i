##############################################################################
# File: 20_2grain_withVoid.i
# File Location: /examples/sintering/paper3/19_massLoss_tests/20_2grain_withVoid
# Created Date: Tuesday June 11th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Wednesday June 12th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Testing the mass conservation bit by bit adding all my stuff into it with
#   a manually defined IC so i can specify C
#
#
##############################################################################

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 60
    ny = 30
    xmin = 0
    xmax = 40000
    ymin = 0
    ymax = 15000
  []
  [subdomain_right]
    type = ParsedSubdomainMeshGenerator
    input = gmg
    combinatorial_geometry = 'x > 20000'
    block_id = 1
  []
  uniform_refine = 1
[]

[GlobalParams]
  profile = TANH
  int_width = 2000
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [w]
    initial_condition = 0
  []
  [c]
  []
  [gr0]
  []
  [gr1]
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
[]

[ICs]
  [c_IC_L]
    type = FunctionIC
    variable = c
    function = ic_func_cGB
    block = 0
  []
  [c_IC_R]
    type = SmoothCircleIC
    variable = c
    x1 = 30000
    y1 = 7500
    radius = 1500
    invalue = 1
    outvalue = 2.424e-06
    block = 1
  []
  [gr0_IC]
    type = FunctionIC
    variable = gr0
    function = ic_func_gr0
    block = '0 1'
  []
  [gr1_IC_L]
    type = FunctionIC
    variable = gr1
    function = ic_func_gr1
    block = '0 1'
  []
  [gr1_IC_R]
    type = SmoothCircleIC
    variable = gr1
    x1 = 30000
    y1 = 7500
    radius = 1500
    invalue = 0
    outvalue = 1
    block = 1
  []
  [phi_IC]
    type = SmoothCircleIC
    variable = phi
    x1 = 30000
    y1 = 7500
    radius = 1500
    invalue = 1
    outvalue = 0.0
    block = '0 1'
  []
[]

[Functions]
  [ic_func_gr0] # Left grain
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r'
    symbol_values = '2000 -5000 7500 15000'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);1-0.5*(1.0-tanh((r-d)/iw))'
  []
  [ic_func_gr1] # Right grain
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r'
    symbol_values = '2000 -5000 7500 15000'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);0.5*(1.0-tanh((r-d)/iw))'
  []
  [ic_func_cGB]
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r cb cgb'
    symbol_values = '2000 -5000 7500 15000 2.424e-06 5.130e-03'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);
    cb + (cgb - cb)*16*((1-0.5*(1.0-tanh((r-d)/iw))) * (0.5*(1.0-tanh((r-d)/iw))))^2'
  []
[]

[Modules]
  [PhaseField]
    [GrandPotential]
      switching_function_names = 'hv hs'
      # Chempot
      chemical_potentials = 'w'
      mobilities = 'cons_mob' #cons_mob
      anisotropic = 'false'
      susceptibilities = 'chi'
      free_energies_w = 'rhov rhos'
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
      mass_conservation = true
      concentrations = 'c'
      hj_over_kVa = 'hoverk_vu hoverk_su' #'hv_over_kVa hs_over_kVa' #
      hj_c_min = 'cvueq_mask csueq_mask' #cvueq_mask=hv*1 'hv_c_min hs_c_min' #
    []
  []
[]

[Materials]
  [consts]
    type = GenericConstantMaterial
    prop_names = 'Va'
    prop_values = '0.04092'
  []
  [k_constants]
    type = GenericConstantMaterial
    prop_names = 'ksu kvu ksi kvi' # Using the GB based values (lowest of mine)
    prop_values = '6.569e2 6.569e2 5.461e7 5.461e7' #'1.302e4 1.302e4 1.092e9 1.092e9' #'26.8 26.8 3.6e21 3.6e21'#154.16 154.16
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
    coupled_variables = 'gr0 gr1'
    expression = '16 * ( (gr0 * gr1)^2 )' #+ (gr0 * etad0)^2 + (gr1 * etad0)^2)
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
  [chiD]
    type = GrandPotentialIsoMaterial
    f_name = chiD
    solid_mobility = L #CHANGED FROM L
    void_mobility = Lv
    chi = chi
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
  [cons_mob]
    type = DerivativeParsedMaterial
    property_name = cons_mob
    # derivative_order = 2
    material_property_names = 'chiD(w,phi,gr0,gr1) Va'
    expression = 'Va * chiD'
  []
  [sintering]
    type = GrandPotentialSinteringMaterial
    chemical_potential = w
    void_op = phi
    Temperature = T
    surface_energy = gb_e_mat
    grainboundary_energy = gb_e_mat
    void_energy_coefficient = kvu
    solid_energy_coefficient = ksu
    solid_energy_model = PARABOLIC
    equilibrium_vacancy_concentration = cv_eq
  []
  [cv_eq]
    type = DerivativeParsedMaterial
    property_name = cv_eq
    # derivative_order = 2
    coupled_variables = 'gr0 gr1 phi w'
    material_property_names = 'hgb(phi,gr0,gr1)' #'rhovi(vac) rhosi hv(phi)'
    constant_names = 'cb cgb'
    # constant_expressions = '3.877e-04 4.347e-03' #Irradiation
    constant_expressions = '2.424e-06 5.130e-03' #No Irradiation- LANL
    expression = 'cgb * hgb + (1 - hgb)*cb'
    outputs = none #'nemesis' # + phi^2
  []
  [cvac]
    type = ParsedMaterial
    property_name = cvac
    material_property_names = 'hs rhos hv rhov Va'
    expression = 'Va*(hs*rhos + hv*rhov)'
    outputs = exodus
  []
  [dvac]
    type = ParsedMaterial
    property_name = dvac
    material_property_names = 'vac_diffus'
    expression = 'vac_diffus'
    outputs = exodus
  []
  # CONSERVATION
  [hoverk_vu]
    type = DerivativeParsedMaterial
    property_name = hoverk_vu
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hv(phi) Va kvu'
    expression = 'hv / (Va * kvu)'
  []
  [hoverk_su]
    type = DerivativeParsedMaterial
    property_name = hoverk_su
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hs(phi) Va ksu'
    expression = 'hs / (Va * ksu)'
  []
  # cvueq_mask = hv*1
  [cvueq_mask]
    type = DerivativeParsedMaterial
    property_name = cvueq_mask
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hv(phi)'
    expression = 'hv * 1'
  []
  [csueq_mask]
    type = DerivativeParsedMaterial
    property_name = csueq_mask
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hs(phi) cv_eq(phi,gr0,gr1)'
    expression = 'hs * cv_eq'
  []
[]

[Postprocessors]
  # [c_total]
  #   type = ElementIntegralVariablePostprocessor
  #   variable = c
  # []
  [phi_total]
    type = ElementIntegralVariablePostprocessor
    variable = phi
  []
  [gr0_total]
    type = ElementIntegralVariablePostprocessor
    variable = gr0
  []
  [gr1_total]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  []
  [hgb_total]
    type = ElementIntegralMaterialProperty
    mat_prop = hgb
  []
  [cvac_total]
    type = ElementIntegralMaterialProperty
    mat_prop = cvac
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
  end_time = 1e6
  # num_steps = 100
  # steady_state_detection = true
  # # From tonks ode input
  automatic_scaling = true
  compute_scaling_once = false
  # line_search = none
  # dt = 1.0
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 0.001
  []
[]

[Outputs]
  csv = true
  exodus = true
  checkpoint = false
[]

# [Debug]
#   show_var_residual_norms = true
# []
