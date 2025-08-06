#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def write_input_pf(dict_user):
    '''
    Write an input file from a template.

    This file is used for the microstructure evolution simulation.
    '''
    print('\nwrite input')
    file_to_write = open('PF_Sintering.i','w')
    file_to_read = open('PF_Sintering_template.i','r')
    lines = file_to_read.readlines()
    file_to_read.close()

    j = 0
    for line in lines:
        j = j + 1
        if j == 4:
            line = line[:-1] + ' ' + str(len(dict_user['L_x'])) + '\n'
        if j == 5:
            line = line[:-1] + ' ' + str(len(dict_user['L_y'])) + '\n'
        if j == 7:
            line = line[:-1] + ' ' + str(dict_user['L_x'][0]) + '\n'
        if j == 8:
            line = line[:-1] + ' ' + str(dict_user['L_x'][-1]) + '\n'
        if j == 9:
            line = line[:-1] + ' ' + str(dict_user['L_y'][0]) + '\n'
        if j == 10:
            line = line[:-1] + ' ' + str(dict_user['L_y'][-1]) + '\n'
        if j == 15:
            line = ''
            for i_grain in range(len(dict_user['L_M_eta'])):
                line = line + '\t[./eta'+str(i_grain)+']\n'+\
                              '\t\torder = FIRST\n'+\
                              '\t\tfamily = LAGRANGE\n'+\
                              '\t\toutputs = exodus\n'+\
                              '\t\t[./InitialCondition]\n'+\
                              '\t\t\ttype = FunctionIC\n'+\
                              '\t\t\tfunction = eta'+str(i_grain)+'_txt\n'+\
                              '\t\t[../]\n'+\
                              '\t[../]\n'
        if j == 30:
            line = ''
            for i_grain in range(len(dict_user['L_M_eta'])):
                # prepare coupled_variables str
                coupled_variables = "'"
                for j_grain in range(len(dict_user['L_M_eta'])):
                    if j_grain != i_grain:
                        coupled_variables = coupled_variables + 'eta'+str(j_grain) + ' '
                coupled_variables = coupled_variables + "c'"
                # write line
                line = line + '\t# Order parameter eta'+str(i_grain)+'\n'+\
                              '\t[./deta'+str(i_grain)+'dt]\n'+\
                              '\t\ttype = TimeDerivative\n'+\
                              '\t\tvariable = eta'+str(i_grain)+'\n'+\
                              '\t[../]\n'+\
                              '\t[./ACBulk_eta'+str(i_grain)+']\n'+\
                              '\t\ttype = AllenCahn\n'+\
                              '\t\tvariable = eta'+str(i_grain)+'\n'+\
                              '\t\tcoupled_variables = '+coupled_variables+'\n'+\
                              '\t\tmob_name = L\n'+\
                              '\t\tf_name = f_tot\n'+\
                              '\t[../]\n'+\
                              '\t[./ACInterface_eta'+str(i_grain)+']\n'+\
                              '\t\ttype = ACInterface\n'+\
                              '\t\tvariable = eta'+str(i_grain)+'\n'+\
                              '\t\tmob_name = L\n'+\
                              '\t\tkappa_name = kappa_eta\n'+\
                              '\t[../]\n'
        if j == 37:
            # prepare coupled_variables str
            coupled_variables = "'"
            for i_grain in range(len(dict_user['L_M_eta'])):
                coupled_variables = coupled_variables + 'eta'+str(i_grain) + ' '
            coupled_variables = coupled_variables + "'"
            # write line
            line = line[:-1] + ' ' + coupled_variables + '\n'
        if j == 61:
            line = line[:-1] + "'" + str(dict_user['L']) + ' ' + str(dict_user['kappa_eta']) + ' ' +\
                                     str(dict_user['M']) + ' ' + str(dict_user['kappa_c']) + "'\n"
        if j == 67:
            # prepare coupled_variables str
            coupled_variables = "'"
            for i_grain in range(len(dict_user['L_M_eta'])):
                coupled_variables = coupled_variables + 'eta'+str(i_grain) + ' '
            coupled_variables = coupled_variables + "c'"
            # write line
            line = line[:-1] + ' ' + coupled_variables + '\n'
        if j == 69:
            line = line[:-1] + "'" + str(dict_user['A']) + ' ' + str(dict_user['B']) + "'\n"
        if j == 70:
            # prepare expression str
            expression = "'A*(c^2)*((1-c)^2) + B*(c^2 + 6*(1-c)*("
            for i_grain in range(len(dict_user['L_M_eta'])):
                expression = expression + 'eta'+str(i_grain) + '^2+'
            expression = expression[:-1] + ') - 4*(2-c)*('   
            for i_grain in range(len(dict_user['L_M_eta'])):
                expression = expression + 'eta'+str(i_grain) + '^3+' 
            expression = expression[:-1] + ') + 3*(' 
            for i_grain in range(len(dict_user['L_M_eta'])):
                expression = expression + 'eta'+str(i_grain) + '^2+' 
            expression = expression[:-1] + ")^2)'"     
            # write line
            line = line[:-1] + ' ' + expression + '\n'
        if j == 78:
            line = ''
            for i_grain in range(len(dict_user['L_M_eta'])):
                line = line + '\t[eta'+str(i_grain)+'_txt]\n'+\
                              '\t\ttype = PiecewiseMultilinear\n'+\
                              '\t\tdata_file = data/eta'+str(i_grain)+'.txt\n'+\
                              '\t[../]\n'
        if j == 104 or j == 105 or j == 108 or j == 109:
            line = line[:-1] + ' ' + str(dict_user['crit_res']) + '\n'
        if j == 112:
            line = line[:-1] + ' ' + str(dict_user['n_ite_max']) + '\n'
        if j == 116:
            line = line[:-1] + ' ' + str(dict_user['dt_PF']) + '\n'      
        if j == 121:
            line = ''
            for i_grain in range(len(dict_user['L_M_eta'])):
                line = line + '\t[eta'+str(i_grain)+'_pp]\n'+\
                              '\t\ttype = ElementAverageValue\n'+\
                              '\t\tvariable = eta'+str(i_grain)+'\n'+\
                              '\t[../]\n'
        if j == 131:
            # prepare variable str
            variable = "'"
            for i_grain in range(len(dict_user['L_M_eta'])):
                variable = variable + 'eta'+str(i_grain) + ' '
            variable = variable[:-1]+"'"  
            # write line
            line = line[:-1] + ' ' + variable + '\n'
        if j == 150:
            # prepare show str
            show = "'c_pp "
            for i_grain in range(len(dict_user['L_M_eta'])):
                show = show + 'eta'+str(i_grain) + '_pp '
            show = show[:-1]+"'"  
            # write line
            line = line[:-1] + ' ' + show + '\n'

        file_to_write.write(line)
    file_to_write.close()


#-------------------------------------------------------------------------------

def write_input_fem(dict_user):
    '''
    Write an input file from a template.

    This file is used for the loading of the microstructure.
    '''
    pass