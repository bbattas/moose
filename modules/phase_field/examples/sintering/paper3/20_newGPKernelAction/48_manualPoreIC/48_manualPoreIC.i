##############################################################################
# File: 48_manualPoreIC.i
# File Location: /examples/sintering/paper3/20_newGPKernelAction/48_manualPoreIC
# Created Date: Tuesday July 30th 2024
# Author: Battas,Brandon Scott (bbattas@ufl.edu)
# -----
# Last Modified: Wednesday July 31st 2024
# Modified By: Battas,Brandon Scott
# -----
# Description:
#  Manually defining the pore with a function IC so i can make cv/ci use the
#   switching function function of that function?
#
#
##############################################################################

# f_dot = 1e-8
# cvic = 500

[Mesh]
  [gmg]
    type = DistributedRectilinearMeshGenerator
    dim = 2
    nx = 80
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
  # [subdomain_botright]
  #   type = ParsedSubdomainMeshGenerator
  #   input = subdomain_right
  #   combinatorial_geometry = 'y < 7500'
  #   block_id = 2
  #   excluded_subdomains = 0
  # []
  uniform_refine = 2
  # second_order = true
[]

[GlobalParams]
  # profile = TANH #Not needed since its all function ICs
  int_width = 1000
  op_num = 2
  var_name_base = gr
[]

[Variables]
  [wvac]
    initial_condition = 0
  []
  [wint]
    initial_condition = 0
  []
  [cvac_var]
  []
  [cint_var]
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
  # C Vac
  [cv_IC_L]
    type = FunctionIC
    variable = cvac_var
    function = ic_func_cvGB
    block = 0
  []
  [cv_IC_R]
    type = FunctionIC
    variable = cvac_var
    function = ic_func_cvPore
    block = 1
  []
  # C INT
  [ci_IC_L]
    type = FunctionIC
    variable = cint_var
    function = ic_func_ciGB
    block = 0
  []
  [ci_IC_R]
    type = FunctionIC
    variable = cint_var
    function = ic_func_ciPore
    block = 1
  []
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
    block = '0' #1
  []
  [gr1_IC_R]
    type = FunctionIC
    variable = gr1
    function = ic_func_GR1pore
    block = '1'
  []
  [phi_IC]
    type = FunctionIC
    variable = phi
    function = ic_func_pore
    block = '0 1'
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
  [cv_bc]
    type = NeumannBC
    variable = 'cvac_var'
    boundary = 'left right top bottom'
    value = 0
  []
  [ci_bc]
    type = NeumannBC
    variable = 'cint_var'
    boundary = 'left right top bottom'
    value = 0
  []
[]

[Functions]
  [ic_func_gr0] # Left grain
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r'
    symbol_values = '1000 -5000 7500 15000'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);1-0.5*(1.0-tanh((r-d)/iw))'
  []
  [ic_func_gr1] # Right grain
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r'
    symbol_values = '1000 -5000 7500 15000'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);0.5*(1.0-tanh((r-d)/iw))'
  []
  [ic_func_pore] # Right side pore
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r'
    symbol_values = '1000 30000 7500 5000'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);
    1-0.5*(1.0-tanh((r-d)/iw))'
  []
  [ic_func_GR1pore] # Right side pore grain value (1-ic_func_pore)
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r'
    symbol_values = '1000 30000 7500 5000'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);
    0.5*(1.0-tanh((r-d)/iw))'
  []
  # C Functions
  [ic_func_cvGB]
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r cb cgb'
    symbol_values = '1000 -5000 7500 15000 2.424e-06 5.130e-03' #3.877e-04 4.347e-03' # irr
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);
    hgbex:=16*((1-0.5*(1.0-tanh((r-d)/iw))) * (0.5*(1.0-tanh((r-d)/iw))))^2;
    hgb:=if(hgbex>1e-6,hgbex,0.0);
    cb + (cgb - cb)*hgb'
  []
  [ic_func_ciGB]
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r cb cgb'
    symbol_values = '1000 -5000 7500 15000 1.667e-32 6.170e-08' #7.258e-09 5.900e-06' # irr
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);
    hgbex:=16*((1-0.5*(1.0-tanh((r-d)/iw))) * (0.5*(1.0-tanh((r-d)/iw))))^2;
    hgb:=if(hgbex>1e-6,hgbex,0.0);
    cb + (cgb - cb)*hgb'
  []
  # C pore
  [ic_func_cvPore] # Right side pore
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r cb cp'
    symbol_values = '1000 30000 7500 5000 2.424e-06 1'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);
    circ:=1-0.5*(1.0-tanh((r-d)/iw));
    hv:=circ*circ*circ*(10 + circ * (-15 + circ * 6));
    hv*cp + (1-hv)*cb'
  []
  [ic_func_ciPore] # Right side pore
    type = ParsedFunction
    symbol_names = 'iw x0 y0 r cb cp'
    symbol_values = '1000 30000 7500 5000 1.667e-32 0'
    expression = 'd:=sqrt((x-x0)^2+(y-y0)^2);
    circ:=1-0.5*(1.0-tanh((r-d)/iw));
    hv:=circ*circ*circ*(10 + circ * (-15 + circ * 6));
    hv*cp + (1-hv)*cb'
  []
[]

# [Modules]
#   [PhaseField]
#     [GrandPotential]
#       switching_function_names = 'hv hs'
#       # Chempot
#       chemical_potentials = 'wvac wint'
#       mobilities = 'chiuD chiiD' #cons_mob
#       anisotropic = 'false false'
#       susceptibilities = 'chiu chii'
#       free_energies_w = 'rhovu rhosu rhovi rhosi'
#       # Grains
#       mobility_name_gr = L_mat
#       kappa_gr = kappa
#       gamma_gr = gamma
#       free_energies_gr = 'omegav omegas'
#       # Other OPs (Phi)
#       additional_ops = 'phi'
#       mobility_name_op = L_mat
#       kappa_op = kappa
#       gamma_grxop = gamma
#       free_energies_op = 'omegav omegas' #empty when no phi'omegaa omegab'
#       # Mass Conservation
#       mass_conservation = true
#       concentrations = 'cvac_var cint_var'
#       hj_over_kVa = 'hoverk_vu hoverk_su' # hoverk_vi hoverk_si' #'hv_over_kVa hs_over_kVa' #
#       hj_c_min = 'cvueq_mask csueq_mask cvieq_mask csieq_mask' #cvueq_mask=hv*1 'hv_c_min hs_c_min' #
#     []
#   []
# []

[Kernels]
  # Order parameter phi
  [DT_phi]
    type = TimeDerivative
    variable = phi
  []
  [ACInt_phi]
    type = ACInterface
    variable = phi
    kappa_name = kappa
    mob_name = L_mat
    coupled_variables = 'gr0 gr1'
  []
  [ACSwitch_phi]
    type = ACSwitching
    variable = phi
    mob_name = L_mat
    Fj_names = 'omegav omegas'
    hj_names = 'hv     hs'
    coupled_variables = 'gr0 gr1 wvac wint' # action technically includes all vars here including phi
  []
  [AcGrGr_phi]
    type = ACGrGrMulti
    variable = phi
    mob_name = L_mat
    v = 'gr0 gr1'
    gamma_names = 'gamma gamma'
  []
  # Order parameter gr0
  [DT_gr0]
    type = TimeDerivative
    variable = gr0
  []
  [ACInt_gr0]
    type = ACInterface
    variable = gr0
    kappa_name = kappa
    mob_name = L_mat
    coupled_variables = 'phi gr1'
  []
  [ACSwitch_gr0]
    type = ACSwitching
    variable = gr0
    mob_name = L_mat
    Fj_names = 'omegav omegas'
    hj_names = 'hv     hs'
    coupled_variables = 'phi gr1 wvac wint' # action technically includes all vars here including gr0
  []
  [AcGrGr_gr0]
    type = ACGrGrMulti
    variable = gr0
    mob_name = L_mat
    v = 'phi gr1'
    gamma_names = 'gamma gamma'
  []
  # Order parameter gr1
  [DT_gr1]
    type = TimeDerivative
    variable = gr1
  []
  [ACInt_gr1]
    type = ACInterface
    variable = gr1
    kappa_name = kappa
    mob_name = L_mat
    coupled_variables = 'phi gr0'
  []
  [ACSwitch_gr1]
    type = ACSwitching
    variable = gr1
    mob_name = L_mat
    Fj_names = 'omegav omegas'
    hj_names = 'hv     hs'
    coupled_variables = 'phi gr0 wvac wint' # action technically includes all vars here including gr1
  []
  [AcGrGr_gr1]
    type = ACGrGrMulti
    variable = gr1
    mob_name = L_mat
    v = 'phi gr0'
    gamma_names = 'gamma gamma'
  []
  # MASS CONSERVATION
  # Concentration Vacancies
  [DT_cv]
    type = TimeDerivative
    variable = cvac_var
  []
  [MatDif_cv]
    type = MatDiffusion
    variable = cvac_var
    v = wvac # in action its 'wvac wint' ?
    diffusivity = chiuD
    args = 'phi gr0 gr1 wvac wint' #NOT SURE IF IT SHOULD HAVE CROSS TERMS?
  []
  # Concentration Interstitials
  [DT_ci]
    type = TimeDerivative
    variable = cint_var
  []
  [MatDif_ci]
    type = MatDiffusion
    variable = cint_var
    v = wint # in action its 'wvac wint' ?
    diffusivity = chiiD
    args = 'phi gr0 gr1 wvac wint'
  []
  # Chemical Potential Vacancies
  [MR_c_wvac]
    type = MatReaction
    variable = wvac
    v = cvac_var # in action its 'cvac_var cint_var' ?
    mob_name = '-1'
  []
  [MR_wvac_hoverk_vu]
    type = MatReaction
    variable = wvac
    args = 'phi gr0 gr1'
    mob_name = hoverk_vu
  []
  [MBD_wvac_cvueq_mask]
    type = MaskedBodyForce
    variable = wvac
    mask = cvueq_mask
    coupled_variables = 'phi gr0 gr1'
  []
  [MR_wvac_hoverk_su]
    type = MatReaction
    variable = wvac
    args = 'phi gr0 gr1'
    mob_name = hoverk_su
  []
  [MBD_wvac_csueq_mask]
    type = MaskedBodyForce
    variable = wvac
    mask = csueq_mask
    coupled_variables = 'phi gr0 gr1'
  []
  # Chemical Potential Interstitials
  [MR_c_wint]
    type = MatReaction
    variable = wint
    v = cint_var # in action its 'cvac_var cint_var' ?
    mob_name = '-1'
  []
  [MR_wint_hoverk_vi]
    type = MatReaction
    variable = wint
    args = 'phi gr0 gr1'
    mob_name = hoverk_vi
  []
  [MBD_wint_cvieq_mask]
    type = MaskedBodyForce
    variable = wint
    mask = cvieq_mask
    coupled_variables = 'phi gr0 gr1'
  []
  [MR_wint_hoverk_si]
    type = MatReaction
    variable = wint
    args = 'phi gr0 gr1'
    mob_name = hoverk_si
  []
  [MBD_wint_csieq_mask]
    type = MaskedBodyForce
    variable = wint
    mask = csieq_mask
    coupled_variables = 'phi gr0 gr1'
  []
  # Sintering Terms for gb/surface energy
  [barrier_phi]
    type = ACBarrierFunction
    variable = phi
    v = 'gr0 gr1' #gr2 gr3'# gr4 gr5'# gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13 gr14 gr15'
    gamma = gamma
    mob_name = L_mat
  []
  [kappa_phi]
    type = ACKappaFunction
    variable = phi
    mob_name = L_mat
    kappa_name = kappa
  []
[]

[Materials]
  [consts]
    type = GenericConstantMaterial
    prop_names = 'Va negOverVa negone' #cvieq_mask cvueq_mask csieq_mask csueq_mask'
    prop_values = '0.04092 -24.4379 -1' #0.0 0.0 0.0 0.0'
  []
  [k_constants]
    type = GenericConstantMaterial
    prop_names = 'ksu kvu ksi kvi' # Using the GB based values (lowest of mine)
    # prop_values = '7.751e2 7.751e2 5.711e5 5.711e5' # Irradiation
    prop_values = '6.569e2 6.569e2 5.461e7 5.461e7' # No Irradiation
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
  [surf_e_mat] # eV/nm^2
    type = ParsedMaterial
    property_name = surf_e_mat
    material_property_names = 'gb_e_mat'
    expression = '2 * gb_e_mat'
    outputs = none
  []
  [hgb]
    type = DerivativeParsedMaterial
    property_name = hgb
    # derivative_order = 2
    coupled_variables = 'gr0 gr1'
    # expression = '16 * ( (gr0 * gr1)^2 )'
    expression = 'hg:=16 * ( (gr0 * gr1)^2 );
                  if(hg>1e-6,hg,0.0)' #+ (gr0 * etad0)^2 + (gr1 * etad0)^2)
    # expression = '4*(1 - (gr0^2 + gr1^2 + etad0^2))^2'
    # outputs = exodus
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
    bulkindex = 1
    gbindex = -1 # -1 sets the GB D to the LANL MD Value in GPIsoMat
    surfindex = -1 #1e11
    vaporindex = 1
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
    vaporindex = 1
    GBwidth = 1.0
    surf_thickness = 1.0 #0.5
    iw_scaling = true
    D_out_name = int_diffus
  []
  # [chiu]
  #   type = DerivativeParsedMaterial
  #   property_name = chiu
  #   # derivative_order = 2
  #   material_property_names = 'chi(wvac,phi,gr0,gr1) Va'
  #   expression = 'Va * chi'
  # []
  # [sintering]
  #   type = GrandPotentialSinteringMaterial
  #   chemical_potential = wvac
  #   void_op = phi
  #   Temperature = T
  #   surface_energy = surf_e_mat #gb_e_mat
  #   grainboundary_energy = gb_e_mat
  #   void_energy_coefficient = kvu
  #   solid_energy_coefficient = ksu
  #   solid_energy_model = PARABOLIC
  #   equilibrium_vacancy_concentration = cv_eq
  # []
  [densificaiton]
    type = GrandPotentialMultiSinteringMaterial
    chemical_potential_vac = wvac
    chemical_potential_int = wint
    void_op = phi
    Temperature = T
    surface_energy = surf_e_mat #gb_e_mat #surf_e_mat #19.7
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
    constant_expressions = '1.667e-32 6.170e-08' #'1.667e-32 6.170e-08' #No Irradiation- LANL
    # constant_expressions = '2.424e-06 5.130e-03' #No Irradiation VACANCY- LANL
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
    # outputs = exodus
  []
  [hoverk_su]
    type = DerivativeParsedMaterial
    property_name = hoverk_su
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hs(phi) Va ksu'
    expression = 'hs / (Va * ksu)'
    # outputs = exodus
  []
  [hoverk_vi]
    type = DerivativeParsedMaterial
    property_name = hoverk_vi
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hv(phi) Va kvi'
    expression = 'hv / (Va * kvi)'
    # outputs = exodus
  []
  [hoverk_si]
    type = DerivativeParsedMaterial
    property_name = hoverk_si
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hs(phi) Va ksi'
    expression = 'hs / (Va * ksi)'
    # outputs = exodus
  []
  # h*ceq Masks
  [cvueq_mask] # cvueq_mask = hv*1
    type = DerivativeParsedMaterial
    property_name = cvueq_mask
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hv(phi)'
    expression = 'hv'
    # outputs = exodus
  []
  [csueq_mask]
    type = DerivativeParsedMaterial
    property_name = csueq_mask
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hs(phi) cv_eq(phi,gr0,gr1)'
    expression = 'hs * cv_eq'
    # outputs = exodus
  []
  [cvieq_mask] # cvieq_mask = hv*0
    type = DerivativeParsedMaterial
    property_name = cvieq_mask
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hv(phi)'
    expression = '0' #'hv' #0 not 1 normally
    # outputs = exodus
  []
  [csieq_mask]
    type = DerivativeParsedMaterial
    property_name = csieq_mask
    coupled_variables = 'phi gr0 gr1'
    derivative_order = 2
    material_property_names = 'hs(phi) ci_eq(phi,gr0,gr1)'
    expression = 'hs * ci_eq'
    # outputs = exodus
  []
  # Extras
  [hv_out]
    type = ParsedMaterial
    property_name = hv_out
    material_property_names = 'hv'
    expression = 'hv'
    outputs = exodus
  []
[]

[Postprocessors]
  [cv_var_total]
    type = ElementIntegralVariablePostprocessor
    variable = cvac_var
  []
  [ci_var_total]
    type = ElementIntegralVariablePostprocessor
    variable = cint_var
  []
  [wvac_total]
    type = ElementIntegralVariablePostprocessor
    variable = wvac
  []
  [wint_total]
    type = ElementIntegralVariablePostprocessor
    variable = wint
  []
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
  # end_time = 1e6
  dt = 0.001
  num_steps = 1
  # steady_state_detection = true
  # # From tonks ode input
  automatic_scaling = true
  compute_scaling_once = false
  # line_search = none
  # dt = 1.0
  # [TimeStepper]
  #   type = IterationAdaptiveDT
  #   optimal_iterations = 6
  #   dt = 0.001
  # []
[]

[Outputs]
  csv = false
  exodus = true
  checkpoint = false
  # file_base = 44_hgbif_forBlockBoundary
[]

# [Debug]
#   show_var_residual_norms = true
# []
