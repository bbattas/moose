##############################################################################
# File: 02_3D_80x80x50_5p.i
# File Location: /examples/sintering/paper3/21_multiApps_IC/11_3D_new_tests/02_3D_80x80x50_5p
# Created Date: Wednesday July 1st 2026
# Author: Brandon Battas (bbattas@ufl.edu)
# -----
# Last Modified: Wednesday July 1st 2026
# Modified By: Brandon Battas
# -----
# Description:
#  80x80x50 instead of 100x100x50 since that didnt even start
#
#
#
##############################################################################

[Mesh]
  [ebsd_mesh]
    type = EBSDMeshGenerator
    filename = 3D_80x80x50nm_8nmavg_test_plusVoid.txt
  []
  [subdomain_external]
    type = ParsedSubdomainMeshGenerator
    input = ebsd_mesh
    combinatorial_geometry = 'x > 98000'
    block_id = 1
  []
  parallel_type = DISTRIBUTED
  uniform_refine = 0
[]

[GlobalParams]
  # profile = TANH
  int_width = 2000
  op_num = 20
  var_name_base = gr
[]

[UserObjects]
  [ebsd_reader]
    type = EBSDReader
  []
  [ebsd]
    type = PolycrystalEBSD
    coloring_algorithm = bt
    ebsd_reader = ebsd_reader
    enable_var_coloring = true
    output_adjacency_matrix = true
    phase = 1
  []
  [grain_tracker]
    type = GrainTracker
    threshold = 0.1 #0.2
    connecting_threshold = 0.09 #0.08
    # compute_halo_maps = false #true#false
    verbosity_level = 1
  []
  #   [terminator]
  #     type = Terminator
  #     expression = 'void_tracker < 2'
  #     execute_on = TIMESTEP_END
  #   []
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
    x1 = 100000
    x2 = 150000
    y1 = -5000
    y2 = 105000
    z1 = -5000
    z2 = 55000
  []
  [voidIC2] #Internal
    type = SpecifiedSmoothCircleIC
    variable = phi
    invalue = 1
    outvalue = 0 #0.01
    radii = '2666.0   2602.0   2414.0   2318.0   2449.0   2459.0   2368.0   2194.0   2191.0   2609.0   2344.0   2680.0   2294.0   2473.0   2695.0   2242.0   2276.0   2460.0   2336.0   2302.0   2664.0   2331.0   2632.0   2494.0   2575.0   2426.0   2325.0   2698.0   2295.0   2607.0   2529.0   2317.0   2724.0   2603.0   2512.0   2481.0   2297.0   2416.0   2626.0   2504.0   2389.0   2420.0   2207.0   2631.0   2681.0   2648.0   2298.0   2338.0   2303.0   2544.0   2566.0   2617.0   2614.0   2307.0   2477.0   2378.0   2609.0   2518.0   2320.0   2431.0   2706.0   2274.0   2308.0   2328.0   2555.0   2720.0   2643.0   2524.0   2299.0   2371.0   2530.0   2183.0   2578.0   2479.0   2658.0   2704.0   2232.0   2718.0   2293.0   2570.0   2650.0   2580.0   2195.0   2610.0   2520.0   2545.0   2628.0   2606.0   2613.0   2470.0   2533.0   2581.0   2473.0   2361.0   2397.0   2440.0   2634.0   2336.0   2360.0   2471.0   2485.0   2468.0   2614.0   2393.0   2713.0   2670.0   2540.0   2314.0   2434.0   2700.0   2478.0   2408.0   2202.0   2206.0   2719.0   2358.0   2548.0   2185.0   2572.0   2495.0   2414.0   2631.0   2379.0   2457.0   2208.0   2718.0   2372.0   2653.0   2391.0   2412.0   2569.0   2479.0   2312.0   2515.0   2673.0   2574.0   2470.0   2218.0   2429.0   2691.0   2523.0   2545.0   2312.0   2425.0   2723.0   2498.0   2310.0   2501.0   2578.0   2579.0   2691.0   2386.0   2630.0   2641.0   2685.0   2587.0   2235.0   2367.0   2721.0   2201.0   2293.0   2637.0   2595.0   2586.0   2579.0   2353.0   2502.0   2530.0   2512.0   2593.0   2420.0   2502.0   2460.0   2703.0   2674.0   2270.0   2630.0   2379.0   2542.0   2589.0   2284.0   2481.0   2469.0   2701.0   2474.0   2405.0   2588.0   2609.0   2690.0   2697.0   2268.0   2708.0   2680.0   2650.0   2373.0   2578.0   2560.0   2614.0   2408.0   2713.0   2631.0   2443.0   2545.0   2516.0   2510.0   2321.0   2679.0   2245.0   2418.0   2724.0   2549.0   2332.0   2245.0   2337.0   2551.0   2642.0   2631.0   2643.0   2646.0   2643.0   2531.0   2282.0   2535.0   2483.0   2483.0   2596.0   2415.0   2412.0   2607.0   2547.0   2289.0   2450.0   2710.0   2344.0   2475.0   2479.0   2318.0   2689.0   2227.0   2355.0   2495.0   2618.0   2519.0   2532.0   2392.0'
    x_positions = '40242.0   48096.0   73530.0   43578.0   46543.0   49314.0   70472.0   68498.0   39078.0   57340.0   37154.0   50634.0   6889.0   14762.0   56484.0   62474.0   65900.0   12526.0   66723.0   20097.0   39906.0   72155.0   12373.0   37176.0   57689.0   16746.0   73310.0   10563.0   36209.0   73573.0   46333.0   24257.0   26572.0   18562.0   31347.0   61839.0   24021.0   36461.0   56387.0   16028.0   5978.0   15604.0   69149.0   22684.0   58815.0   35407.0   5694.0   39502.0   32461.0   53615.0   67980.0   71190.0   43372.0   31263.0   34768.0   63364.0   17586.0   51293.0   7269.0   5893.0   57048.0   68343.0   28283.0   5924.0   45657.0   15217.0   28318.0   53029.0   14837.0   25643.0   47690.0   38649.0   66102.0   14531.0   73246.0   51073.0   62998.0   48290.0   17155.0   27537.0   10021.0   30884.0   63995.0   12666.0   65391.0   24576.0   57873.0   37561.0   47746.0   66402.0   38571.0   66526.0   39618.0   14053.0   43263.0   71215.0   27267.0   39483.0   53208.0   18576.0   65523.0   60149.0   25187.0   19222.0   18074.0   38150.0   31751.0   14694.0   20892.0   25592.0   66479.0   25162.0   50606.0   8914.0   58126.0   22806.0   71848.0   67074.0   16250.0   31846.0   26691.0   49422.0   7223.0   44815.0   52667.0   13301.0   24348.0   24887.0   43307.0   68787.0   27433.0   47070.0   17952.0   28871.0   48268.0   5962.0   70151.0   72210.0   61720.0   44117.0   43758.0   23599.0   22569.0   70601.0   32694.0   52559.0   54253.0   51318.0   33425.0   62419.0   34849.0   42971.0   22639.0   39374.0   51924.0   72595.0   62849.0   31627.0   14823.0   42496.0   17617.0   32829.0   15212.0   15110.0   35223.0   30304.0   62165.0   35952.0   61495.0   38576.0   6643.0   30648.0   21133.0   17804.0   16672.0   29769.0   56621.0   31885.0   8232.0   41384.0   56832.0   73096.0   36938.0   31935.0   29738.0   66715.0   29728.0   25671.0   10417.0   10621.0   29015.0   43784.0   20741.0   43957.0   35222.0   7610.0   64327.0   46968.0   53265.0   65973.0   52252.0   12018.0   8780.0   37177.0   19625.0   20376.0   44629.0   70789.0   49381.0   37737.0   19871.0   51450.0   60362.0   7470.0   14755.0   50854.0   23978.0   33745.0   23729.0   20167.0   72998.0   13702.0   8829.0   48866.0   32178.0   50282.0   51846.0   18273.0   6412.0   55671.0   46674.0   8379.0   44058.0   72733.0   8812.0   21445.0   6361.0   38638.0   27642.0   53175.0   59405.0   59293.0   7186.0   45069.0   59069.0'
    y_positions = '15070.0   55541.0   74063.0   70056.0   42274.0   36923.0   12253.0   46891.0   71421.0   53538.0   57177.0   22348.0   6618.0   6119.0   25742.0   33004.0   59872.0   11058.0   40063.0   67255.0   26885.0   38942.0   73919.0   7086.0   11429.0   25432.0   38383.0   55706.0   71491.0   10736.0   27358.0   64775.0   9311.0   50687.0   34038.0   26323.0   70222.0   25679.0   36121.0   69069.0   24361.0   73905.0   64344.0   61607.0   10592.0   48658.0   15550.0   54845.0   69614.0   30512.0   19555.0   11465.0   49354.0   71726.0   41187.0   58824.0   63991.0   11183.0   57408.0   38491.0   29913.0   50819.0   40680.0   28938.0   62876.0   32039.0   6649.0   19215.0   50732.0   49679.0   37792.0   16526.0   27512.0   33765.0   22221.0   69869.0   5664.0   67078.0   48075.0   42466.0   52061.0   66585.0   51405.0   30547.0   74291.0   60410.0   24646.0   6944.0   40652.0   37793.0   73092.0   40891.0   21595.0   45945.0   52299.0   23243.0   51272.0   41408.0   7307.0   73541.0   28946.0   36871.0   25803.0   54619.0   18230.0   6346.0   11857.0   10120.0   38353.0   14141.0   26412.0   51476.0   74423.0   65243.0   8705.0   15548.0   28929.0   59398.0   10693.0   28660.0   69890.0   40180.0   51586.0   65370.0   61764.0   65605.0   6959.0   36852.0   15042.0   6027.0   38468.0   5694.0   50711.0   55720.0   13987.0   69305.0   17698.0   45501.0   20118.0   30331.0   25218.0   23495.0   12355.0   33970.0   21502.0   13507.0   42862.0   33882.0   23155.0   20489.0   36188.0   67751.0   7055.0   28795.0   19555.0   26034.0   50360.0   9390.0   23539.0   46882.0   13429.0   17723.0   41383.0   39905.0   30811.0   64227.0   34362.0   62006.0   8847.0   58287.0   46913.0   61295.0   42754.0   65610.0   58950.0   57507.0   18123.0   12834.0   18932.0   32834.0   68374.0   31126.0   45242.0   18054.0   57957.0   68733.0   31095.0   21244.0   64701.0   67578.0   48090.0   34064.0   34257.0   11707.0   73631.0   42558.0   16122.0   21805.0   57628.0   11837.0   48363.0   33801.0   65150.0   54255.0   23875.0   42929.0   45314.0   47172.0   27035.0   50394.0   18174.0   47727.0   69732.0   33541.0   22364.0   55337.0   28103.0   36646.0   73144.0   34745.0   7564.0   74242.0   18123.0   47006.0   42138.0   15410.0   27164.0   57434.0   20083.0   68176.0   73812.0   60386.0   9424.0   39992.0   39774.0   54310.0   5718.0   15821.0   74537.0   64866.0   47870.0   39329.0   24585.0   58615.0   48655.0'
    z_positions = '25411.0   28724.0   39207.0   43198.0   13898.0   31435.0   42664.0   13860.0   8247.0   34994.0   16987.0   17887.0   44083.0   21517.0   25950.0   27793.0   7328.0   30469.0   26174.0   35099.0   22466.0   12062.0   7128.0   9479.0   6656.0   23701.0   38996.0   10664.0   15298.0   9152.0   29279.0   27380.0   23665.0   12904.0   40026.0   32558.0   21878.0   15166.0   6469.0   27794.0   22645.0   35541.0   39025.0   15180.0   23527.0   35519.0   22834.0   28725.0   31260.0   17083.0   33824.0   30324.0   15379.0   43698.0   8317.0   14769.0   7047.0   6931.0   43827.0   32497.0   39009.0   27087.0   42099.0   33148.0   42969.0   37438.0   40537.0   42061.0   6098.0   33126.0   5839.0   36882.0   18061.0   18150.0   14440.0   34799.0   9140.0   8278.0   21185.0   16426.0   19290.0   22161.0   36741.0   29948.0   37681.0   36087.0   13086.0   26469.0   38779.0   34108.0   35752.0   6519.0   43110.0   43927.0   39032.0   43467.0   40207.0   33658.0   44220.0   12603.0   39774.0   39577.0   21635.0   42884.0   39559.0   34942.0   16169.0   38322.0   43264.0   36571.0   9488.0   18954.0   42303.0   44501.0   32405.0   30013.0   24998.0   31190.0   13955.0   30390.0   13170.0   22727.0   30502.0   29887.0   38812.0   13465.0   8293.0   34112.0   5711.0   38767.0   6355.0   26637.0   34588.0   28855.0   34283.0   14041.0   24822.0   31559.0   40961.0   13439.0   5847.0   29779.0   44407.0   5600.0   23780.0   15227.0   34253.0   43084.0   7955.0   22515.0   29924.0   15344.0   33088.0   36078.0   9403.0   35402.0   17819.0   31958.0   15714.0   7412.0   6416.0   30923.0   14564.0   33926.0   7367.0   6329.0   11819.0   33111.0   41347.0   44235.0   37120.0   42840.0   26776.0   42827.0   31497.0   6033.0   32354.0   5462.0   15414.0   43212.0   6137.0   16539.0   43303.0   43627.0   16013.0   31543.0   14803.0   14519.0   33483.0   21271.0   7501.0   21524.0   23693.0   42103.0   23847.0   10324.0   8599.0   36556.0   13574.0   17827.0   17998.0   9454.0   5858.0   9516.0   6082.0   5604.0   28663.0   40900.0   43776.0   21882.0   20351.0   42487.0   14312.0   39946.0   31461.0   5836.0   39415.0   19270.0   30187.0   11356.0   22917.0   43253.0   37822.0   5914.0   25447.0   25155.0   5811.0   22948.0   6981.0   27088.0   22387.0   22769.0   16948.0   18224.0   21851.0   5667.0   17108.0   14137.0   5418.0   17363.0   10416.0   18489.0   43636.0   20657.0   25681.0'
    block = 0
  []
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
  # THEQ
  # Vacancy mat to variable
  [cvac_aux_th]
    family = LAGRANGE
    order = FIRST
  []
  [cvac_aux_elem_th]
    family = MONOMIAL
    order = CONSTANT
  []
  # Interstitial mat to variable
  [cint_aux_th]
    family = LAGRANGE
    order = FIRST
  []
  [cint_aux_elem_th]
    family = MONOMIAL
    order = CONSTANT
  []
  # IRR EQ
  # Vacancy mat to variable
  [cvac_aux_ir]
    family = LAGRANGE
    order = FIRST
  []
  [cvac_aux_elem_ir]
    family = MONOMIAL
    order = CONSTANT
  []
  # Interstitial mat to variable
  [cint_aux_ir]
    family = LAGRANGE
    order = FIRST
  []
  [cint_aux_elem_ir]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  # THEQ
  # Vacancy
  [cvac_aux_elem_th]
    type = MaterialRealAux
    variable = 'cvac_aux_elem_th'
    property = 'cvac_th'
  []
  [cvac_aux_th]
    type = ProjectionAux
    variable = cvac_aux_th
    v = cvac_aux_elem_th
  []
  # Interstitial
  [cint_aux_elem_th]
    type = MaterialRealAux
    variable = 'cint_aux_elem_th'
    property = 'cint_th'
  []
  [cint_aux_th]
    type = ProjectionAux
    variable = cint_aux_th
    v = cint_aux_elem_th
  []
  # THEQ
  # Vacancy
  [cvac_aux_elem_ir]
    type = MaterialRealAux
    variable = 'cvac_aux_elem_ir'
    property = 'cvac_ir'
  []
  [cvac_aux_ir]
    type = ProjectionAux
    variable = cvac_aux_ir
    v = cvac_aux_elem_ir
  []
  # Interstitial
  [cint_aux_elem_ir]
    type = MaterialRealAux
    variable = 'cint_aux_elem_ir'
    property = 'cint_ir'
  []
  [cint_aux_ir]
    type = ProjectionAux
    variable = cint_aux_ir
    v = cint_aux_elem_ir
  []
[]

[BCs]
  [phi_bc]
    type = NeumannBC
    variable = 'phi'
    boundary = 'left right top bottom'
    value = 0
  []
  [grain_bc]
    type = NeumannBC
    variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8 gr9 gr10 gr11 gr12 gr13 gr14 gr15 gr16 gr17 gr18 gr19'
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
      # hj_over_kVa = 'hv_over_kvVa hs_over_kvVa hv_over_kiVa hs_over_kiVa' #'hoverk_vu hoverk_su hoverk_vi hoverk_si'
      # hj_c_min = 'hv_cv_min hs_cv_min hv_ci_min hs_ci_min' #'cvueq_mask csueq_mask cvieq_mask csieq_mask'
    []
  []
[]

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
  [hgb]
    type = SwitchingFunctionGBMaterial
    h_name = hgb
    grain_ops = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7'
    hgb_threshold = 0.0001
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
    surfindex = 5.8422e10 #-1 #1e11
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
    surfindex = 1.3989e11 #-1 #100 #1e11
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
    type = UO2CeqMaterial
    ceq_name = cv_eq
    hgb = hgb
    vi = VAC_TH
  []
  [ci_eq]
    type = UO2CeqMaterial
    ceq_name = ci_eq
    hgb = hgb
    vi = INT_TH
  []
  [cv_eq_ir]
    type = UO2CeqMaterial
    ceq_name = cv_eq_ir
    hgb = hgb
    vi = VAC_IRR
  []
  [ci_eq_ir]
    type = UO2CeqMaterial
    ceq_name = ci_eq_ir
    hgb = hgb
    vi = INT_IRR
  []
  # Outputs for visualization
  [cvac_th]
    type = ParsedMaterial
    property_name = cvac_th
    material_property_names = 'hs cv_eq hv Va'
    expression = 'cv:=(hs*cv_eq ) + hv*1.0;
                  if(cv<0.0, 0.0, cv)'
    # outputs = nemesis
  []
  [cint_th]
    type = ParsedMaterial
    property_name = cint_th
    material_property_names = 'hs ci_eq hv Va'
    expression = 'ci:=(hs*ci_eq) + hv*0.0;
                  if(ci<0.0, 0.0, ci)'
    # outputs = nemesis
  []
  [cvac_ir]
    type = ParsedMaterial
    property_name = cvac_ir
    material_property_names = 'hs cv_eq_ir hv Va'
    expression = 'cv:=(hs*cv_eq_ir ) + hv*1.0;
                  if(cv<0.0, 0.0, cv)'
    # outputs = nemesis
  []
  [cint_ir]
    type = ParsedMaterial
    property_name = cint_ir
    material_property_names = 'hs ci_eq_ir hv Va'
    expression = 'ci:=(hs*ci_eq_ir) + hv*0.0;
                  if(ci<0.0, 0.0, ci)'
    # outputs = nemesis
  []
  # # hgb check for interface width
  # [hgb_out]
  #   type = ParsedMaterial
  #   property_name = hgb_out
  #   material_property_names = 'hgb'
  #   expression = 'hgb'
  #   outputs = nemesis
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
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -sub_pc_type -pc_asm_overlap -sub_pc_factor_shift_type'
  petsc_options_value = ' asm      lu           2                nonzero'
  nl_max_its = 30
  l_max_its = 60
  l_tol = 1e-04 #4
  nl_rel_tol = 1e-6 #6 #default is 1e-8
  # nl_abs_tol = 1e-14 #only needed when near equilibrium or veeeery small dt
  # start_time = 0
  end_time = 10000 #1e8
  # num_steps = 0
  # steady_state_detection = true
  # # From tonks ode input
  automatic_scaling = true
  compute_scaling_once = false
  # line_search = none
  # dt = 1.0
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 6
    dt = 0.1
    linear_iteration_ratio = 1e5
  []
[]

[Outputs]
  csv = true
  exodus = false
  nemesis = true
  checkpoint = false
  # file_base = MC_GifRif_highDs
[]

# [Debug]
#   show_var_residual_norms = true
# []
