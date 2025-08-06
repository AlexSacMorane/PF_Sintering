# fine no flux, working simulation with 4 grains grainTracker included
# from https://github.com/idaholab/moose/discussions/20805
[Debug]
  show_var_residual_norms = true
[]

[Mesh]
# make mesh size 1nm*1nm
    type = GeneratedMesh
    dim = 2
    nx = 400
    ny = 200
    xmin = 0
    xmax =400
    ymin = 0
    ymax = 200

[]
[GlobalParams]
    op_num = 4
    var_name_base = gr
  # let's output all material properties for demonstration purposes
  #outputs = exodus
[]
#
# These AuxVariables hold the directly calculated free energy density in the
# simulation cell. They are provided for visualization purposes.
#
[AuxVariables]
  [./local_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./cross_energy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./bnds]
  # Auxivariables that define the location of grain boundaries in a polycrystalline sample
  [../]

[]

[AuxKernels]
  [./local_free_energy]
    type = TotalFreeEnergy
    variable = local_energy
    interfacial_vars = 'c'
    kappa_names = 'kappa_gr'
    # additional_free_energy = cross_energy
  [../]
  [./bnds]
  # Calculate location of grain boundaries in a polycrystalline sample
    type = BndsCalcAux
    variable = bnds
    var_name_base = gr
    op_num = 4
    v = 'gr0 gr1 gr2 gr3'
  [../]



[]

[Variables]
  [./c]
    order = FIRST
    family = LAGRANGE
  [../]

  # Automatically creates order parameter variables
  [./PolycrystalVariables]
  [../]

  [./w]               # functional derivative of the total free energy
    order = FIRST
    family = LAGRANGE
  [../]

[]

[ICs]

# TWO more here
[./gr0_IC]
    type = SpecifiedSmoothCircleIC
    x_positions = '80  240 '
    y_positions = '60  60   '
    z_positions = '0  0'
    radii =       '40 40'
    variable = gr0
    invalue = 1.0
    outvalue = 0.0
    int_width = 6
[../]
[./gr1_IC]
    type = SpecifiedSmoothCircleIC
    x_positions = '160  320 '
    y_positions = '60  60   '
    z_positions = '0  0'
    radii =       '40 40'
    variable = gr1
    invalue = 1.0
    outvalue = 0.0
    int_width = 6
[../]
[./gr2_IC]
    type = SpecifiedSmoothCircleIC
    x_positions = '80  240 '
    y_positions = '140  140   '
    z_positions = '0  0'
    radii =       '40 40'
    variable = gr2
    invalue = 1.0
    outvalue = 0.0
    int_width = 6
[../]
[./gr3_IC]
    type = SpecifiedSmoothCircleIC
    x_positions = '160  320'
    y_positions = '140  140   '
    z_positions = '0  0'
    radii =       '40 40'
    variable = gr3
    invalue = 1.0
    outvalue = 0.0
    int_width = 6
[../]




 [./c] #fmol
 type = SpecifiedSmoothCircleIC
 x_positions = '80 160 80 160 240 240 320 320'
 y_positions = '60 60 140 140 60   140   60  140'
 z_positions = '0   0   0 0   0    0   0'
 radii =       '40 40 40 40 40 40 40 40'
 variable = c
 invalue = 1.0
 outvalue = 0.0
 int_width = 6
  [../]

[]


[Functions]
  [./load_x]
    # Defines the force on the grains in the x-direction
    type = ParsedFunction
    value = 0.00*cos(x*pi/600)
  [../]
  [./load_y]
    # Defines the force on the grains in the y-direction
    type = ConstantFunction
    value = 0.0
  [../]

[]

[./NoRBMMultiKernel]
# Creates all of the necessary Allen Cahn kernels automatically
c = c
f_name = F # uses customized Free energy
mob_name = L
kappa_name = kappa_gr
[../]

[Kernels]




    # split
    [./c_res]
      type = SplitCHParsed
      variable = c
      f_name = F
      kappa_name = kappa_c
      args = 'gr0 gr1 gr2 gr3'
      w = w
    [../]
    [./w_res]
      type = SplitCHWRes
      variable = w
      mob_name = M
    [../]
    [./time]
      type = CoupledTimeDerivative
      variable = w
      v = c
    [../]



[]

[BCs]
  [./bcs]
    #zero flux BC
    type = NeumannBC
    value = 0
    variable = c
    boundary = '0 1 2 3'
  [../]


[]

[Materials]
  # here we declare some of the model parameters: the mobilities and interface
  # gradient prefactors. For this example we use arbitrary numbers. In an actual simulation
  # physical mobilities would be used, and the interface gradient prefactors would
  # be readjusted to the free energy magnitudes.
  [./consts]
    type = GenericConstantMaterial
    prop_names  = '   L  '
    prop_values = '   1   '
  [../]

  [./constsDiff]
    type = GenericConstantMaterial
    #              nm^2/ns
    prop_names  = 'M_gb   M_s   M_b'
    prop_values = '0.37e1    3.98e1   0.013e1  '
  [../]

  [./consts2]
    type = GenericConstantMaterial
    #              g/mol       g/nm^3   K         ev/K-mol
     prop_names  = 'MM       density    T         R  '
    prop_values = '97.474     4.1e-21  1000  5.1896e19  '
  [../]

  [./freeEnergyConsts2]
    type = GenericConstantMaterial
    # ev/nm^3
    prop_names = 'B_f                C_f'
    prop_values = '27.46      2.4968'
  [../]

  [./constants]
    type = GenericConstantMaterial
    # ev/nm
    prop_names = 'kappa_c   kappa_gr    '

    prop_values = '14.9808     7.4904   '

  [../]

  # h(eta)
  [./h_eta]
    type = SwitchingFunctionMaterial
    h_order = HIGH
    eta = c
  [../]

  # 1- h(eta), putting in function explicitly
  [./one_minus_h_eta_explicit]
    type = DerivativeParsedMaterial
    f_name = one_minus_h_explicit
    args = c
    function = 1-c^3*(6*c^2-15*c+10)
    outputs = exodus
  [../]


  [./CHmobility]    # Cahn and Hilliard Mobility incluing all diffusin path (m^5/J/s)
    type = DerivativeParsedMaterial
    f_name = M
    args = 'c gr0 gr1 gr2 gr3'
    function = 'M_b*c^3*(10 - 15*c + 6*c^2)
                + M_s*c*(1 - c)
                + M_gb*(gr0*gr1+gr0*gr2+gr0*gr3+gr1*gr2+gr1*gr3+gr2*gr3)'
    material_property_names = 'M_b M_s M_gb '
    derivative_order = 2
  [../]


  [./chemical_free_energy]
      type = DerivativeParsedMaterial
      f_name = F # to be added with mechanical
      args = 'c gr0 gr1 gr2 gr3 ' #Must be changed as op_num changes. Copy/paste from line 4
      function = 'B_f*c^2*(1-c)^2+C_f*(c^2+6*(1-c)*(gr0^2+gr1^2+gr2^2+gr3^2)
                  -4*(2-c)*(gr0^3+gr1^3+gr2^3+gr3^3)
                  +3*(gr0^2+gr1^2+gr2^2+gr3^2)^2)'
                                   #Copy/paste from lines 5-6
      derivative_order = 2
      material_property_names='B_f C_f'
  [../]

  [./force_density]
    type = ExternalForceDensityMaterial
    c = c
    k = 10.0
    force_x = load_x
    force_y = load_y
  [../]
[]





[Preconditioning]

  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
    type = Transient
    scheme = bdf2
    solve_type = PJFNK
    automatic_scaling=true
    petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_ksp_type
                           -sub_pc_type -pc_asm_overlap'
    petsc_options_value = 'asm      31                  preonly
                           ilu          2'
    l_tol = 1e-05
    nl_max_its = 30
    l_max_its = 30
    nl_rel_tol = 1e-07
    nl_abs_tol = 1e-09
    start_time = 0.0
    end_time = 1000
    dt = 1e-1
[]

[Debug]
  # show_var_residual_norms = true
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
  checkpoint=true
[]
