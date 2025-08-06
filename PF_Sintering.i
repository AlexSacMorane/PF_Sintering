[Mesh]
    type = GeneratedMesh
    dim = 2
    nx = 301
    ny = 301
    nz = 0
    xmin = -30.0
    xmax = 29.999999999999787
    ymin = -30.0
    ymax = 29.999999999999787
    elem_type = QUAD4
[]

[Variables]
	[./eta0]
		order = FIRST
		family = LAGRANGE
		outputs = exodus
		[./InitialCondition]
			type = FunctionIC
			function = eta0_txt
		[../]
	[../]
	[./eta1]
		order = FIRST
		family = LAGRANGE
		outputs = exodus
		[./InitialCondition]
			type = FunctionIC
			function = eta1_txt
		[../]
	[../]
	[./eta2]
		order = FIRST
		family = LAGRANGE
		outputs = exodus
		[./InitialCondition]
			type = FunctionIC
			function = eta2_txt
		[../]
	[../]
	[./eta3]
		order = FIRST
		family = LAGRANGE
		outputs = exodus
		[./InitialCondition]
			type = FunctionIC
			function = eta3_txt
		[../]
	[../]
	[./eta4]
		order = FIRST
		family = LAGRANGE
		outputs = exodus
		[./InitialCondition]
			type = FunctionIC
			function = eta4_txt
		[../]
	[../]
	[./eta5]
		order = FIRST
		family = LAGRANGE
		outputs = exodus
		[./InitialCondition]
			type = FunctionIC
			function = eta5_txt
		[../]
	[../]
	[./eta6]
		order = FIRST
		family = LAGRANGE
		outputs = exodus
		[./InitialCondition]
			type = FunctionIC
			function = eta6_txt
		[../]
	[../]
	[./eta7]
		order = FIRST
		family = LAGRANGE
		outputs = exodus
		[./InitialCondition]
			type = FunctionIC
			function = eta7_txt
		[../]
	[../]
	[./eta8]
		order = FIRST
		family = LAGRANGE
		outputs = exodus
		[./InitialCondition]
			type = FunctionIC
			function = eta8_txt
		[../]
	[../]
    [./c]
        outputs = exodus
        [./InitialCondition]
            type = FunctionIC
            function = c_txt
        [../]
    [../]
    [./w] # functional derivative of the total free energy  
        order = FIRST
        family = LAGRANGE
    [../]
[]

[Kernels]
	# Order parameter eta0
	[./deta0dt]
		type = TimeDerivative
		variable = eta0
	[../]
	[./ACBulk_eta0]
		type = AllenCahn
		variable = eta0
		coupled_variables = 'eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 c'
		mob_name = L
		f_name = f_tot
	[../]
	[./ACInterface_eta0]
		type = ACInterface
		variable = eta0
		mob_name = L
		kappa_name = kappa_eta
	[../]
	# Order parameter eta1
	[./deta1dt]
		type = TimeDerivative
		variable = eta1
	[../]
	[./ACBulk_eta1]
		type = AllenCahn
		variable = eta1
		coupled_variables = 'eta0 eta2 eta3 eta4 eta5 eta6 eta7 eta8 c'
		mob_name = L
		f_name = f_tot
	[../]
	[./ACInterface_eta1]
		type = ACInterface
		variable = eta1
		mob_name = L
		kappa_name = kappa_eta
	[../]
	# Order parameter eta2
	[./deta2dt]
		type = TimeDerivative
		variable = eta2
	[../]
	[./ACBulk_eta2]
		type = AllenCahn
		variable = eta2
		coupled_variables = 'eta0 eta1 eta3 eta4 eta5 eta6 eta7 eta8 c'
		mob_name = L
		f_name = f_tot
	[../]
	[./ACInterface_eta2]
		type = ACInterface
		variable = eta2
		mob_name = L
		kappa_name = kappa_eta
	[../]
	# Order parameter eta3
	[./deta3dt]
		type = TimeDerivative
		variable = eta3
	[../]
	[./ACBulk_eta3]
		type = AllenCahn
		variable = eta3
		coupled_variables = 'eta0 eta1 eta2 eta4 eta5 eta6 eta7 eta8 c'
		mob_name = L
		f_name = f_tot
	[../]
	[./ACInterface_eta3]
		type = ACInterface
		variable = eta3
		mob_name = L
		kappa_name = kappa_eta
	[../]
	# Order parameter eta4
	[./deta4dt]
		type = TimeDerivative
		variable = eta4
	[../]
	[./ACBulk_eta4]
		type = AllenCahn
		variable = eta4
		coupled_variables = 'eta0 eta1 eta2 eta3 eta5 eta6 eta7 eta8 c'
		mob_name = L
		f_name = f_tot
	[../]
	[./ACInterface_eta4]
		type = ACInterface
		variable = eta4
		mob_name = L
		kappa_name = kappa_eta
	[../]
	# Order parameter eta5
	[./deta5dt]
		type = TimeDerivative
		variable = eta5
	[../]
	[./ACBulk_eta5]
		type = AllenCahn
		variable = eta5
		coupled_variables = 'eta0 eta1 eta2 eta3 eta4 eta6 eta7 eta8 c'
		mob_name = L
		f_name = f_tot
	[../]
	[./ACInterface_eta5]
		type = ACInterface
		variable = eta5
		mob_name = L
		kappa_name = kappa_eta
	[../]
	# Order parameter eta6
	[./deta6dt]
		type = TimeDerivative
		variable = eta6
	[../]
	[./ACBulk_eta6]
		type = AllenCahn
		variable = eta6
		coupled_variables = 'eta0 eta1 eta2 eta3 eta4 eta5 eta7 eta8 c'
		mob_name = L
		f_name = f_tot
	[../]
	[./ACInterface_eta6]
		type = ACInterface
		variable = eta6
		mob_name = L
		kappa_name = kappa_eta
	[../]
	# Order parameter eta7
	[./deta7dt]
		type = TimeDerivative
		variable = eta7
	[../]
	[./ACBulk_eta7]
		type = AllenCahn
		variable = eta7
		coupled_variables = 'eta0 eta1 eta2 eta3 eta4 eta5 eta6 eta8 c'
		mob_name = L
		f_name = f_tot
	[../]
	[./ACInterface_eta7]
		type = ACInterface
		variable = eta7
		mob_name = L
		kappa_name = kappa_eta
	[../]
	# Order parameter eta8
	[./deta8dt]
		type = TimeDerivative
		variable = eta8
	[../]
	[./ACBulk_eta8]
		type = AllenCahn
		variable = eta8
		coupled_variables = 'eta0 eta1 eta2 eta3 eta4 eta5 eta6 eta7 c'
		mob_name = L
		f_name = f_tot
	[../]
	[./ACInterface_eta8]
		type = ACInterface
		variable = eta8
		mob_name = L
		kappa_name = kappa_eta
	[../]
    # Order parameter c
    [./c_res]
        type = SplitCHParsed
        variable = c
        f_name = f_tot
        kappa_name = kappa_c
        coupled_variables =  'eta0 eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 '
        w = w
    [../]
    # variable w (derivative of the total free energy)
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

[]

[Materials]
    [./consts]
        type = GenericConstantMaterial
        prop_names  = 'L kappa_eta M kappa_c'
        prop_values ='1 0.5 1 1'
    [../]
    [./free_energy]
        type = DerivativeParsedMaterial
        block = 0
        property_name = f_tot
        coupled_variables = 'eta0 eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8 c'
        constant_names = 'A B'
        constant_expressions = '16 1'
        expression =  'A*(c^2)*((1-c)^2) + B*(c^2 + 6*(1-c)*(eta0^2+eta1^2+eta2^2+eta3^2+eta4^2+eta5^2+eta6^2+eta7^2+eta8^2) - 4*(2-c)*(eta0^3+eta1^3+eta2^3+eta3^3+eta4^3+eta5^3+eta6^3+eta7^3+eta8^3) + 3*(eta0^2+eta1^2+eta2^2+eta3^2+eta4^2+eta5^2+eta6^2+eta7^2+eta8^2)^2)'
        enable_jit = true
        derivative_order = 2
        #outputs = exodus
    [../]
[]

[Functions]
	[eta0_txt]
		type = PiecewiseMultilinear
		data_file = data/eta0.txt
	[../]
	[eta1_txt]
		type = PiecewiseMultilinear
		data_file = data/eta1.txt
	[../]
	[eta2_txt]
		type = PiecewiseMultilinear
		data_file = data/eta2.txt
	[../]
	[eta3_txt]
		type = PiecewiseMultilinear
		data_file = data/eta3.txt
	[../]
	[eta4_txt]
		type = PiecewiseMultilinear
		data_file = data/eta4.txt
	[../]
	[eta5_txt]
		type = PiecewiseMultilinear
		data_file = data/eta5.txt
	[../]
	[eta6_txt]
		type = PiecewiseMultilinear
		data_file = data/eta6.txt
	[../]
	[eta7_txt]
		type = PiecewiseMultilinear
		data_file = data/eta7.txt
	[../]
	[eta8_txt]
		type = PiecewiseMultilinear
		data_file = data/eta8.txt
	[../]
    [c_txt]
        type = PiecewiseMultilinear
        data_file = data/c.txt
    []
[]

[Preconditioning]
    # This preconditioner makes sure the Jacobian Matrix is fully populated. Our
    # kernels compute all Jacobian matrix entries.
    # This allows us to use the Newton solver below.
    [./SMP]
        type = SMP
        full = true
    [../]
[]

[Executioner]
    type = Transient
    scheme = 'bdf2'

    # Automatic differentiation provides a _full_ Jacobian in this example
    # so we can safely use NEWTON for a fast solve
    solve_type = 'NEWTON'

    l_max_its = 20
    l_tol = 0.0001
    l_abs_tol = 0.0001

    nl_max_its = 10
    nl_rel_tol = 0.0001
    nl_abs_tol = 0.0001

    start_time = 0.0
    num_steps = 100

    [./TimeStepper]
        type = SolutionTimeAdaptiveDT
        dt = 0.01
    [../]
[]

[Postprocessors]
	[eta0_pp]
		type = ElementAverageValue
		variable = eta0
	[../]
	[eta1_pp]
		type = ElementAverageValue
		variable = eta1
	[../]
	[eta2_pp]
		type = ElementAverageValue
		variable = eta2
	[../]
	[eta3_pp]
		type = ElementAverageValue
		variable = eta3
	[../]
	[eta4_pp]
		type = ElementAverageValue
		variable = eta4
	[../]
	[eta5_pp]
		type = ElementAverageValue
		variable = eta5
	[../]
	[eta6_pp]
		type = ElementAverageValue
		variable = eta6
	[../]
	[eta7_pp]
		type = ElementAverageValue
		variable = eta7
	[../]
	[eta8_pp]
		type = ElementAverageValue
		variable = eta8
	[../]
    [c_pp]
      type = ElementAverageValue
      variable = c
    []
[]

[UserObjects]
    [./grain_tracker]
        type = GrainTracker
        variable =  'eta0 eta1 eta2 eta3 eta4 eta5 eta6 eta7 eta8'
        outputs = none
        compute_var_to_feature_map = true
        execute_on = 'initial timestep_begin'
    [../]
[]

[Outputs]
    execute_on = 'initial timestep_end'
    exodus = true
    #[./other]
    #    type = VTK
    #    execute_on = 'TIMESTEP_END'
    #[../]
    [console]
        type = Console
        execute_on = 'nonlinear'
        all_variable_norms = true
        max_rows = 3
        show =  'c_pp eta0_pp eta1_pp eta2_pp eta3_pp eta4_pp eta5_pp eta6_pp eta7_pp eta8_pp'
    []
    [./csv]
        type = CSV
        execute_on = 'TIMESTEP_END'
    [../]
[]
