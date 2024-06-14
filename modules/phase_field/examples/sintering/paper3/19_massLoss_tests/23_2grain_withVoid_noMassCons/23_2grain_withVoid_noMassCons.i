##############################################################################
# File: 23_2grain_withVoid_noMassCons.i
# File Location: /examples/sintering/paper3/19_massLoss_tests/23_2grain_withVoid_noMassCons
# Created Date: Wednesday June 12th 2024
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Wednesday June 12th 2024
# Modified By: Brandon Battas
# -----
# Description:
#  Input 22 (using wvac and wint as variables, both with mass conservation
#   BUT no interaction between them (no GAM)) WITHOUT THE MASS CONSERVATION
#  For comparison to 22 mass loss/growth
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
  [wvac]
    initial_condition = 0
  []
  # [cvac_var]
  # []
  [wint]
    initial_condition = 0
  []
  # [cint_var]
  # []
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
  # # C Vac
  # [cv_IC_L]
  #   type = FunctionIC
  #   variable = cvac_var
  #   function = ic_func_cvGB
  #   block = 0
  # []
  # [cv_IC_R]
  #   type = SmoothCircleIC
  #   variable = cvac_var
  #   x1 = 30000
  #   y1 = 7500
  #   radius = 1500
  #   invalue = 1
  #   outvalue = 2.424e-06
  #   block = 1
  # []
  # # C Int
  # [ci_IC_L]
  #   type = FunctionIC
  #   variable = cint_var
  #   function = ic_func_ciGB
  #   block = 0
  # []
  # [ci_IC_R]
  #   type = SmoothCircleIC
  #   variable = cint_var
  #   x1 = 30000
  #   y1 = 7500
  #   radius = 1500
  #   invalue = 0
  #   outvalue = 1.667e-32
  #   block = 1
  # []
  # Grains
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
  [ic_func_cvGB]
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r cb cgb'
    symbol_values = '2000 -5000 7500 15000 2.424e-06 5.130e-03'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);
    cb + (cgb - cb)*16*((1-0.5*(1.0-tanh((r-d)/iw))) * (0.5*(1.0-tanh((r-d)/iw))))^2'
  []
  [ic_func_ciGB]
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r cb cgb'
    symbol_values = '2000 -5000 7500 15000 1.667e-32 6.170e-08'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);
    cb + (cgb - cb)*16*((1-0.5*(1.0-tanh((r-d)/iw))) * (0.5*(1.0-tanh((r-d)/iw))))^2'
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
      # mass_conservation = true
      # concentrations = 'cvac_var cint_var'
      # hj_over_kVa = 'hoverk_vu hoverk_su hoverk_vi hoverk_si' #'hv_over_kVa hs_over_kVa' #
      # hj_c_min = 'cvueq_mask csueq_mask cvieq_mask csieq_mask' #cvueq_mask=hv*1 'hv_c_min hs_c_min' #
    []
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
    mass_conservation = true
  []
  [cv_eq]
    type = DerivativeParsedMaterial
    property_name = cv_eq
    # derivative_order = 2
    coupled_variables = 'gr0 gr1 phi wvac'
    material_property_names = 'hgb(phi,gr0,gr1)' #'rhovi(vac) rhosi hv(phi)'
    constant_names = 'cb cgb'
    # constant_expressions = '3.877e-04 4.347e-03' #Irradiation
    constant_expressions = '2.424e-06 5.130e-03' #No Irradiation- LANL
    expression = 'cgb * hgb + (1 - hgb)*cb'
    outputs = none #'nemesis' # + phi^2
  []
  [ci_eq]
    type = DerivativeParsedMaterial
    property_name = ci_eq
    # derivative_order = 2
    coupled_variables = 'gr0 gr1 phi wint'
    material_property_names = 'hgb(phi,gr0,gr1)' # 'rhovi(wint) rhosi(wint) hv(phi)'
    constant_names = 'cb cgb'
    # constant_expressions = '7.258e-09 5.900e-06' #Irradiation
    constant_expressions = '1.667e-32 6.170e-08' #No Irradiation- LANL
    expression = 'cgb * hgb + (1 - hgb)*cb'
    outputs = none #'nemesis' #+ phi^2
  []
  [cvac]
    type = ParsedMaterial
    property_name = cvac
    material_property_names = 'hs rhosu hv rhovu Va'
    expression = 'Va*(hs*rhosu + hv*rhovu)'
    outputs = exodus
  []
  [cint]
    type = ParsedMaterial
    property_name = cint
    material_property_names = 'hs rhosi hv rhovi Va'
    expression = 'Va*(hs*rhosi + hv*rhovi)'
    outputs = exodus
  []
  # CONSERVATION
  # h / (V*k)
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
  [hoverk_vi]
    type = DerivativeParsedMaterial
    property_name = hoverk_vi
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hv(phi) Va kvi'
    expression = 'hv / (Va * kvi)'
  []
  [hoverk_si]
    type = DerivativeParsedMaterial
    property_name = hoverk_si
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hs(phi) Va ksi'
    expression = 'hs / (Va * ksi)'
  []
  # h*ceq Masks
  [cvueq_mask] # cvueq_mask = hv*1
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
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hs(phi) ci_eq(phi,gr0,gr1)'
    expression = 'hs * ci_eq'
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
  [cint_total]
    type = ElementIntegralMaterialProperty
    mat_prop = cint
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
