import numpy as np
import itertools
import datetime
import sys
from scipy.stats import t
from scipy.stats import f

from omars_analysis.subset_selection import get_soe
from omars_analysis.generate_model_matrix_heredity import create_model_matrix_heredity

# from subset_selection import get_soe
# from generate_model_matrix_heredity import create_model_matrix_heredity

def hat_matrix(smat):
    hat = smat @ np.linalg.pinv(smat.T @ smat) @ smat.T
    return hat

def create_quadratic_interactions(smat):
    no_runs, nf = smat.shape[0], smat.shape[1]
    quad_ix = [i for i,x in enumerate(smat.T) if (abs(x)!=1).any()]
    
    # Quadratic
    if len(quad_ix) > 0:
        quadratic = np.absolute(smat[:,quad_ix])
        soe_mat = np.hstack((np.ones((no_runs,1)),quadratic))
    else:
        soe_mat = np.ones((no_runs,1))

    # Interactions
    for combi in itertools.combinations(range(nf), 2):
        col = (smat[:,combi[0]]*smat[:,combi[1]]).reshape((no_runs,1))
        soe_mat = np.hstack((soe_mat, col))

    #centering (important for later functions)
    soe_mat = soe_mat - soe_mat.mean(axis=0)

    return soe_mat[:,1:], quad_ix

def code(smat):
    a = np.ones((smat.shape[0],1))
    for i in smat.T:
        maxi = max(i)
        mini = min(i)
        avg = (maxi + mini)/2
        diff = (maxi - mini)/2
        temp = (i - avg)/diff
        temp = temp.reshape((smat.shape[0],1))
        a = np.hstack((a,temp))
       
    return a[:,1:]

def get_omars_analysis(smat: np.ndarray, cy: np.ndarray, alpha: list = [0.05, 0.2, 0.2], qheredity: str='n', iheredity: str='n', effects_to_drop: list=[], full: str='y', force_me: list=[], user_limit_for_step_two=None):
    
    web_statements=[]

    i_start = datetime.datetime.now()
    mat = code(smat)
    no_runs, nf = mat.shape[0], mat.shape[1]

    # Centering response in case not centered already
    cy = (cy - cy.mean(axis=0)).reshape(no_runs,1)

    # Calculate error variance
    intercept = np.ones((no_runs, 1))
    mat_qi, quad_ix = create_quadratic_interactions(mat)
    total_mat = np.hstack((np.hstack((intercept, mat)),mat_qi))
    
    new_error_variance_projection = np.identity(no_runs) - (total_mat @ np.linalg.pinv(total_mat))
    denominator = no_runs - np.linalg.matrix_rank(total_mat)
    print('\nInitial error degrees of freedom available - {}'.format(denominator))
    web_statements.append('\nInitial error degrees of freedom available - {}'.format(denominator))

    if denominator == 0:
        print('Analysis cannot continue due to lack of available degrees of freedom to estimate error variance')
        sys.exit()

    new_variance = (cy.T @ new_error_variance_projection @ cy)/denominator
    print('\nInitial estimate of error variance - {}'.format(round(float(new_variance),3)))
    web_statements.append('\nInitial estimate of error variance - {}'.format(round(float(new_variance),3)))
    
    s = np.sqrt(new_variance)
    
    '''
    ANALYSIS OF ME (STEP 1)
    '''
    mebeta = (np.linalg.inv(mat.T @ mat) @ mat.T @ cy).reshape(1,nf)
    vdiags = np.matrix.diagonal(np.linalg.inv(mat.T @ mat))
    svdiags = np.sqrt(vdiags).reshape(1,nf)
    tvalue = np.absolute(mebeta/(s*svdiags))
    ttest = t.ppf(1 - (alpha[0]/2), df=denominator)
    tcomp = tvalue/ttest
    
    indexrme = [i for i,x in enumerate(tcomp[0]) if (abs(x)>=1)] # index of active main effects    
    indexfme = [i for i,x in enumerate(tcomp[0]) if (abs(x)<1)] # index of inactive main effects
    
    new_me = [str(g+1) for g in indexrme]
    print('\nActive main effects - {}'.format(new_me))
    web_statements.append('\nActive main effects - {}'.format(new_me))
    # p_values_ttest = [t.cdf(ttest,denominator)]
    p_value_ttest = 1-t.cdf(tvalue,denominator)
    print('\nThe p-values for double sided t-tests for main effects - {} (threshold alpha/2 = {})'.format(p_value_ttest[0], alpha[0]/2))
    web_statements.append('\nThe p-values for double sided t-tests for main effects - {} (threshold alpha/2 = {})'.format(p_value_ttest[0], alpha[0]/2))

    if full == 'n':
        return web_statements
    
    '''
    END OF STEP 1
    '''

    # Force main effects in step two:
    if len(force_me)>0:
        print('\nForced main effect - {}'.format(force_me))
        web_statements.append('\nForced main effect - {}'.format(force_me))
        ix = [int(x) for x in force_me]
        for x in ix:
            if x-1 in indexrme:
                pass
            elif x-1 in indexfme:
                indexrme.append(x-1)

    indexrme.sort()
    indexfme = [x for x in range(nf) if x not in indexrme]
    inactive = mat[:, indexfme]

    # Update error variance
    if len(indexfme)>0:
        inactive_proj = hat_matrix(inactive)
    else:
        inactive_proj = np.zeros((no_runs, no_runs))
    new_variance_projection_updated = new_error_variance_projection + inactive_proj
    numerator = cy.T @ new_variance_projection_updated @ cy
    new_denominator = denominator + len(indexfme)
    new_variance = numerator/new_denominator
    print('\nUpdated estimate of error variance taking into account inactive main effects- {}'.format(round(float(new_variance),3)))
    web_statements.append('\nUpdated estimate of error variance taking into account inactive main effects- {}'.format(round(float(new_variance),3)))
    '''
    ANALYSIS OF SOE (STEP 2)
    '''
    # Second order proection matrix (make sure to only multiply with centered response since intercept is ignored)
    y_second_projection = np.identity(no_runs) - hat_matrix(mat) - new_error_variance_projection
    cy_second = (y_second_projection @ cy).reshape(no_runs,1)
    TSS = cy_second.T @ cy_second # Total sum of squares (for second order space)
    
    ####        F Test for second order effects
    o_dfleft = np.linalg.matrix_rank(mat_qi) # quadratics need to be centered in mat_qi
    ftest = f.ppf(q= 1-alpha[1], dfn=o_dfleft, dfd=new_denominator)
    F1 = (TSS/o_dfleft)/new_variance
    fcomp = F1/ftest
    p_value_intial = 1-f.cdf(fcomp, o_dfleft, new_denominator)
    print('\nThe p-value for the first F-test is {} (threshold alpha value - {})'.format(round(float(p_value_intial),3),alpha[1]))
    web_statements.append('\nThe p-value for the first F-test is {} (threshold alpha value - {})'.format(round(float(p_value_intial),3),alpha[1]))
    '''
    CREATE SECOND ORDER MATRIX WITH HEREDITY FOR ANALYSIS (WITHOUT INTERCEPT)
    '''
    if fcomp >= 1:
        # standardize columns
        std = 'n'
        soe_mat2, my_dictionary_soe_single  = create_model_matrix_heredity(std, mat, quad_ix, indexrme, indexfme, qheredity, iheredity)
        
        ix_del = [my_dictionary_soe_single[x] for x in effects_to_drop]
        soe_mat2 = np.delete(soe_mat2, ix_del, axis=1)
        new_names = [x for x in my_dictionary_soe_single.keys() if x not in effects_to_drop]
        my_dictionary_soe_single = {new_names[i]:i for i in range(len(new_names))}
        
        # if quadratic_special_coding == 'y':
        #     relevant_q = []
        #     for dumzy in my_dictionary_soe_single:
        #         if dumzy.split('_')[0] == dumzy.split('_')[1]:
        #             relevant_q.append(int(dumzy.split('_')[0])-1)
                    
        #     if len(relevant_q)>0:
        #         temp_mat = mat[:,relevant_q]
        #         temp_mat = temp_mat - temp_mat.mean(axis=0) # center
        #         temp_mat = np.square(temp_mat) # square
        #         temp_mat = code(temp_mat) # normalize
        #         soe_mat2[:,0:len(relevant_q)] = temp_mat.copy()
        active_soe, p_value_omars, actual_soe_df, concluding_statement = get_soe(mat, o_dfleft, soe_mat2, cy_second, new_denominator, new_variance, my_dictionary_soe_single, alpha_val=alpha[2], limit_for_step_two=user_limit_for_step_two)

        print('\nThe p-value for the final F-test is \n {} (threshold alpha value - {})'.format(round(p_value_omars,3), alpha[2]))
        web_statements.append('\nThe p-value for the final F-test is \n {} (threshold alpha value - {})'.format(round(p_value_omars,3), alpha[2]))
    else:
        active_soe = []
        p_value_omars = np.nan
    
    new_qe = []
    new_ie = []
    for d in active_soe:
        x = d.split('_')
        if x[0] == x[1]:
            new_qe.append(d)
        else:
            new_ie.append(d)

    print('\nActive interaction effects - {}'.format(new_ie))
    print('\nActive quadratic effects - {}'.format(new_qe))

    web_statements.append('\nActive interaction effects - {}'.format(new_ie))
    web_statements.append('\nActive quadratic effects - {}'.format(new_qe))
    
    i_finish = datetime.datetime.now()
    timing = (i_finish - i_start).total_seconds()
    print('\nAnalysis performed in {} seconds'.format(round(timing,3)))

    print('\n=============================================================')
    print('\nInformation on the limit of second order terms allowed to enter the model:')
    print(concluding_statement)
    print('Number of second order effects considered - {}'.format(soe_mat2.shape[1]))
    print('Maximum number of second order terms jointly estimable of all the second order effects considered - {}'.format(actual_soe_df))
    print('\nRank of matrix with all {} possible second order terms is {} (in case the rank is less than the number of total second order terms, it is advisable to set a limit that is less than half of the rank)\n'.format(mat_qi.shape[1], o_dfleft))
    
    web_statements.append('\nAnalysis performed in {} seconds'.format(round(timing,3)))
    web_statements.append('\n=============================================================')
    web_statements.append(concluding_statement)
    web_statements.append('Number of second order effects considered - {}'.format(soe_mat2.shape[1]))
    web_statements.append('Maximum number of second order terms jointly estimable of all the second order effects considered - {}'.format(actual_soe_df))
    web_statements.append('\nRank of matrix with all {} possible second order terms is {} (in case the rank is less than the number of total second order terms, it is advisable to set a limit that is less than half of the rank)\n'.format(mat_qi.shape[1], o_dfleft))
    
    return web_statements

# if __name__ == '__main__':
#     file = np.loadtxt('omars_analysis/data/Laser_data.txt')

#     matz = file[:,:-1]
#     resp = file[:,[-1]]

#     # output = get_omars_analysis(smat=matz, cy=resp, qheredity='y', iheredity='n', effects_to_drop=[], force_me=['4'])
#     output = get_omars_analysis(smat=matz, cy=resp, alpha=[0.05, 0.2, 0.2], qheredity='n', iheredity='n', effects_to_drop=[], full='y', force_me=['4'], user_limit_for_step_two=None)